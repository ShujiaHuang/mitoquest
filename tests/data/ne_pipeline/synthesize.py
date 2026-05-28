#!/usr/bin/env python3
"""
Synthesize a small mother-child mtDNA cohort for end-to-end testing of
`mitoquest trans-prep` and `mitoquest ne-estimate`.

The output is fully deterministic (seeded with --seed) and consists of:

  * cohort.vcf                       multi-sample VCF (GT-aligned FORMAT/AD,
                                     Number=., as emitted by `mitoquest caller`)
  * cohort.fam                       PLINK FAM with one trio per mother-child pair
  * cohort.transmission_pairs.tsv    optional pre-baked TSV mimicking the
                                     output of `mitoquest trans-prep`, for
                                     users who want to run ne-estimate alone

Default model:
  * 20 mother-child pairs (40 samples + 20 placeholder fathers in FAM only)
  * 12 mtDNA sites with maternal heteroplasmy ~ Uniform[0.15, 0.85]
  * True bottleneck Ne = 5 (single-generation transmission)
  * Read depth 2000 for both mother and child
  * 5% of (sample, site) cells emit a truncated GT='.' call (GT missing,
    DP positive, AD absent) so that the GT-aware drop logic in
    `mitoquest trans-prep` is exercised end-to-end.

Usage:
    python3 synthesize.py                # writes files in this directory
    python3 synthesize.py --true-ne 10   # tighter / looser bottleneck
    python3 synthesize.py --help

Author: Shujia Huang (hshujia@qq.com)
Date:   2026-05-28
"""

from __future__ import annotations

import argparse
import os
import random
from pathlib import Path
from typing import List, Tuple

# ---------------------------------------------------------------------
# Helpers (no numpy dependency to keep this script standalone).
# ---------------------------------------------------------------------

def binomial(n: int, p: float, rng: random.Random) -> int:
    """Sample k ~ Binomial(n, p).  O(n) but fine for the depths used here."""
    p = max(0.0, min(1.0, p))
    return sum(1 for _ in range(n) if rng.random() < p)


# ---------------------------------------------------------------------
# Synthesizer
# ---------------------------------------------------------------------

def synthesize(out_dir: Path,
               n_pairs: int,
               n_sites: int,
               true_ne: int,
               m_dp: int,
               c_dp: int,
               vaf_low: float,
               vaf_high: float,
               seed: int,
               n_multiallelic: int = 0,
               missing_gt_rate: float = 0.0) -> None:
    rng = random.Random(seed)
    missing_gt_rate = max(0.0, min(1.0, missing_gt_rate))

    # ----- 1. sample names + FAM -------------------------------------------------
    samples: List[str] = []
    fam_lines: List[str] = []
    pairs: List[Tuple[str, str, str]] = []   # (fam_id, mother_id, child_id)
    for i in range(1, n_pairs + 1):
        fam_id   = f"FAM{i:03d}"
        mom_id   = f"mom_{i:03d}"
        child_id = f"child_{i:03d}"
        # PLINK FAM: family child father mother sex phenotype
        fam_lines.append(f"{fam_id} {child_id} 0 {mom_id} 2 -9")
        samples.extend([mom_id, child_id])
        pairs.append((fam_id, mom_id, child_id))

    # ----- 2. variant sites + per-pair read counts ------------------------------
    # Each site is biallelic by default.  The first `n_multiallelic` sites are
    # promoted to tri-allelic SNVs (REF + 2 ALTs), exercising the per-ALT loop
    # in `mitoquest trans-prep`.
    bases = "ACGT"
    sites: List[Tuple[int, str, List[str]]] = []   # (pos, ref, [alt1, alt2, ...])
    n_multiallelic = max(0, min(n_multiallelic, n_sites))
    for s in range(n_sites):
        pos = 200 + s * (16000 // max(1, n_sites))
        ref = bases[rng.randrange(4)]
        # First ALT (always present).
        alt1 = bases[(bases.index(ref) + 1 + rng.randrange(3)) % 4]
        if s < n_multiallelic:
            # Pick a second ALT distinct from REF and ALT1.
            others = [b for b in bases if b != ref and b != alt1]
            alt2 = others[rng.randrange(len(others))]
            sites.append((pos, ref, [alt1, alt2]))
        else:
            sites.append((pos, ref, [alt1]))

    # vcf_rows[i_site] = list of "GT:DP:AD" strings ordered as `samples`
    vcf_rows: List[List[str]] = [[] for _ in sites]
    # Per-(site, pair, alt-index) records for the pre-baked TSV.
    tsv_records: List[Tuple[int, str, str, str, str, str,
                            int, int, int, float,
                            int, int, int, float]] = []

    for i_site, (pos, ref, alts) in enumerate(sites):
        n_alts = len(alts)
        for fam_id, mom_id, child_id in pairs:
            # Maternal allele frequencies for REF and each ALT (sum to 1).
            # For biallelic, p_ref = 1 - p_alt with p_alt drawn from [low, high].
            # For tri-allelic, split [vaf_low, vaf_high] across the two ALTs.
            if n_alts == 1:
                p_alts = [rng.uniform(vaf_low, vaf_high)]
            else:
                # Two ALT frequencies summing to a value in [vaf_low, vaf_high].
                total = rng.uniform(vaf_low, vaf_high)
                share = rng.uniform(0.3, 0.7)
                p_alts = [total * share, total * (1.0 - share)]
            p_ref = max(0.0, 1.0 - sum(p_alts))

            # Bottleneck: sample true_ne copies, distribute by category.
            categories = [p_ref] + p_alts          # length = 1 + n_alts
            counts = [0] * len(categories)
            for _ in range(true_ne):
                u = rng.random()
                acc = 0.0
                idx = 0
                for ci, pc in enumerate(categories):
                    acc += pc
                    if u <= acc:
                        idx = ci
                        break
                else:
                    idx = len(categories) - 1
                counts[idx] += 1
            # child p for each ALT after bottleneck:
            p_c_alts = [counts[1 + a] / true_ne for a in range(n_alts)]

            # Read-level multinomial sampling at depth m_dp / c_dp.
            mom_alt_reads = [0] * n_alts
            for _ in range(m_dp):
                u = rng.random()
                acc = 0.0
                hit = -1
                for ai, pa in enumerate(p_alts):
                    acc += pa
                    if u <= acc:
                        hit = ai
                        break
                if hit >= 0:
                    mom_alt_reads[hit] += 1
            mom_ref_reads = m_dp - sum(mom_alt_reads)

            child_alt_reads = [0] * n_alts
            for _ in range(c_dp):
                u = rng.random()
                acc = 0.0
                hit = -1
                for ai, pa in enumerate(p_c_alts):
                    acc += pa
                    if u <= acc:
                        hit = ai
                        break
                if hit >= 0:
                    child_alt_reads[hit] += 1
            child_ref_reads = c_dp - sum(child_alt_reads)

            # Build the per-sample VCF FORMAT field (GT-aligned AD layout,
            # i.e. AD[i] is the depth of the allele at GT position i).  This
            # mirrors what `mitoquest caller` emits: AD only carries depths
            # for alleles present in that sample's GT, in GT order.  Alleles
            # absent from GT contribute zero reads in the downstream model.
            full_counts_mom   = [mom_ref_reads]   + mom_alt_reads        # len = 1 + n_alts
            full_counts_child = [child_ref_reads] + child_alt_reads

            def _gt_aligned(full_counts: List[int]) -> Tuple[str, str]:
                # Allele indices observed (>0 reads), preserving 0=REF first
                # then ALT indices in ascending order.  Empty support → GT='.'.
                in_gt = [i for i, c in enumerate(full_counts) if c > 0]
                if not in_gt:
                    return (".", "")
                gt_str = "/".join(str(i) for i in in_gt) if len(in_gt) > 1 \
                         else str(in_gt[0])
                ad_str = ",".join(str(full_counts[i]) for i in in_gt)
                return (gt_str, ad_str)

            mom_gt,   mom_ad   = _gt_aligned(full_counts_mom)
            child_gt, child_ad = _gt_aligned(full_counts_child)

            # Inject GT='.' calls at the configured rate.  This mirrors
            # `mitoquest caller`'s truncated emission for samples with no
            # variant evidence, and exercises the GT-aware drop logic in
            # `mitoquest trans-prep`.
            mom_missing   = (rng.random() < missing_gt_rate)
            child_missing = (rng.random() < missing_gt_rate)
            if mom_missing:   mom_gt   = "."
            if child_missing: child_gt = "."

            # Truncate the FORMAT string when GT='.' to mirror the caller's
            # behaviour: under our test FORMAT `GT:DP:AD`, the equivalent of
            # the caller's `.:0:DP` (under its full `GT:GQ:DP:AD:...`) is
            # `.:DP:.`  — GT missing, DP positive, AD absent.
            mom_field   = (f".:{m_dp}:." if mom_gt   == "."
                           else f"{mom_gt}:{m_dp}:{mom_ad}")
            child_field = (f".:{c_dp}:." if child_gt == "."
                           else f"{child_gt}:{c_dp}:{child_ad}")
            vcf_rows[i_site].append(mom_field)
            vcf_rows[i_site].append(child_field)

            # Skip pre-baked TSV rows where either sample has GT='.':
            # `mitoquest trans-prep` would drop the (mother, child) pair at
            # this site (counted as `pair_site_dropped_gt_missing`).
            if mom_gt == "." or child_gt == ".":
                continue

            # One TSV row per ALT (mirroring `trans-prep` per-ALT decomposition).
            for a, alt in enumerate(alts):
                m_alt = mom_alt_reads[a]
                c_alt = child_alt_reads[a]
                tsv_records.append((
                    pos, ref, alt, fam_id, mom_id, child_id,
                    m_dp, m_dp - m_alt, m_alt, m_alt / m_dp if m_dp else 0.0,
                    c_dp, c_dp - c_alt, c_alt, c_alt / c_dp if c_dp else 0.0,
                ))

    # ----- 3. write VCF ---------------------------------------------------------
    vcf_lines: List[str] = [
        "##fileformat=VCFv4.2",
        "##FILTER=<ID=PASS,Description=\"All filters passed\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">",
        # AD is declared as Number=. (per-sample variable) because the layout
        # is GT-aligned: each sample carries one AD value per allele present
        # in that sample's GT, not one per ALT (Number=A) or per REF+ALTs
        # (Number=R).  This matches `mitoquest caller`'s actual emission.
        "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depth "
        "for each allele in the order listed by GT\">",
        "##contig=<ID=chrM,length=16569>",
        "#" + "\t".join(["CHROM", "POS", "ID", "REF", "ALT",
                         "QUAL", "FILTER", "INFO", "FORMAT", *samples]),
    ]
    for (pos, ref, alts), row in zip(sites, vcf_rows):
        vcf_lines.append("\t".join([
            "chrM", str(pos), ".", ref, ",".join(alts),
            "100", "PASS", ".", "GT:DP:AD", *row,
        ]))
    (out_dir / "cohort.vcf").write_text("\n".join(vcf_lines) + "\n")

    # ----- 4. write FAM ---------------------------------------------------------
    (out_dir / "cohort.fam").write_text("\n".join(fam_lines) + "\n")

    # ----- 5. write a pre-baked TSV mimicking `mitoquest trans-prep` output ----
    tsv_header = (
        "CHROM\tPOS\tREF\tALT\tFAM_ID\tMOTHER_ID\tCHILD_ID\t"
        "MOTHER_DP\tMOTHER_AD_REF\tMOTHER_AD_ALT\tMOTHER_VAF\t"
        "CHILD_DP\tCHILD_AD_REF\tCHILD_AD_ALT\tCHILD_VAF\tQC"
    )
    provenance = (
        "#synthesize.py: deterministic mtDNA mother-child cohort\n"
        f"#synthesize.py: n_pairs={n_pairs} n_sites={n_sites} "
        f"true_ne={true_ne} m_dp={m_dp} c_dp={c_dp} "
        f"vaf=[{vaf_low},{vaf_high}] seed={seed} "
        f"n_multiallelic={n_multiallelic} missing_gt_rate={missing_gt_rate}\n"
    )
    tsv_lines = [provenance + tsv_header]
    for (pos, ref, alt, fam_id, mom_id, child_id,
         m_dp_v, m_ad_ref, m_ad_alt, m_vaf,
         c_dp_v, c_ad_ref, c_ad_alt, c_vaf) in tsv_records:
        tsv_lines.append("\t".join([
            "chrM", str(pos), ref, alt, fam_id, mom_id, child_id,
            str(m_dp_v), str(m_ad_ref), str(m_ad_alt), f"{m_vaf:.6f}",
            str(c_dp_v), str(c_ad_ref), str(c_ad_alt), f"{c_vaf:.6f}",
            "PASS",
        ]))
    (out_dir / "cohort.transmission_pairs.tsv").write_text("\n".join(tsv_lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir",   default=os.path.dirname(os.path.abspath(__file__)),
                   help="Where to write cohort.{vcf,fam,transmission_pairs.tsv}.")
    p.add_argument("--n-pairs",   type=int,   default=20)
    p.add_argument("--n-sites",   type=int,   default=12)
    p.add_argument("--true-ne",   type=int,   default=5,
                   help="True bottleneck size used for simulation.")
    p.add_argument("--m-dp",      type=int,   default=2000)
    p.add_argument("--c-dp",      type=int,   default=2000)
    p.add_argument("--vaf-low",   type=float, default=0.15)
    p.add_argument("--vaf-high",  type=float, default=0.85)
    p.add_argument("--n-multiallelic", type=int, default=2,
                   help="Of n_sites, how many should be tri-allelic SNVs "
                        "(REF + 2 ALTs).  Set to 0 for purely biallelic data.")
    p.add_argument("--missing-gt-rate", type=float, default=0.05,
                   help="Probability that any (sample, site) cell emits a "
                        "truncated GT='.' call (GT missing, DP positive, AD "
                        "absent), exercising the GT-aware drop logic in "
                        "`mitoquest trans-prep`.  Set to 0 to disable.")
    p.add_argument("--seed",      type=int,   default=20260528)
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    synthesize(out_dir,
               n_pairs=args.n_pairs,
               n_sites=args.n_sites,
               true_ne=args.true_ne,
               m_dp=args.m_dp,
               c_dp=args.c_dp,
               vaf_low=args.vaf_low,
               vaf_high=args.vaf_high,
               seed=args.seed,
               n_multiallelic=args.n_multiallelic,
               missing_gt_rate=args.missing_gt_rate)
    print(f"Wrote cohort.vcf, cohort.fam, cohort.transmission_pairs.tsv "
          f"to {out_dir}")


if __name__ == "__main__":
    main()
