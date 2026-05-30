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

Bottleneck models:
  * discrete   (default) - sample `round(true_ne)` integer copies via
                            Multinomial; matches `--model discrete` MMLE.
  * continuous           - draw post-bottleneck heteroplasmy from a
                            Beta-diffusion (Wright-Fisher diffusion limit),
                            i.e. `p_child ~ Beta(p_m*(Ne-1), (1-p_m)*(Ne-1))`;
                            matches `--model continuous` MMLE.  Permits a
                            non-integer `--true-ne`.

Outlier injection (--outlier-frac):
  A fraction of (site, pair) cells is forced to extreme drift
  (`p_child` set to 0 or 1 with equal probability) to simulate NUMTs,
  mosaicism, or genotyping errors.  Used to demonstrate the
  outlier sensitivity of the Kimura method-of-moments estimator.

Multi-generation transmission (--n-generations):
  Iterates the bottleneck `g` times before read sampling, simulating a
  pedigree of depth `g` (e.g. grand-mother -> grand-child uses g=2).
  Wright-Fisher predicts that the per-pair accumulated drift is
      F_g = 1 - (1 - 1/Ne)^g,
  so feeding multi-generation pairs into a single-generation estimator
  will yield an apparent `Ne_apparent = 1 / F_g`, biased *downward*
  relative to the true per-generation Ne.  Useful for studying the
  pedigree-depth bias of all three estimators.

Usage:
    python3 synthesize.py                          # default cohort
    python3 synthesize.py --true-ne 10             # tighter / looser bottleneck
    python3 synthesize.py --bottleneck-model continuous --true-ne 7.5
    python3 synthesize.py --outlier-frac 0.05      # 5% high-drift outliers
    python3 synthesize.py --n-generations 5        # 5-generation pedigree
    python3 synthesize.py --help

Author: Shujia Huang (hshujia@qq.com)
Date:   2026-05-28
"""

from __future__ import annotations

import argparse
import math
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


def gamma_sample(shape: float, rng: random.Random) -> float:
    """Sample X ~ Gamma(shape, 1) using Marsaglia-Tsang (shape >= 1) or the
    boost-shape trick for shape < 1 (Stuart's theorem).  No numpy needed.
    """
    if shape <= 0.0:
        return 0.0
    if shape < 1.0:
        # Stuart: if Y ~ Gamma(shape+1, 1) and U ~ Uniform(0,1), then
        # Y * U^(1/shape) ~ Gamma(shape, 1).
        u = rng.random()
        return gamma_sample(shape + 1.0, rng) * (u ** (1.0 / shape))
    d = shape - 1.0 / 3.0
    c = 1.0 / (9.0 * d) ** 0.5
    while True:
        x = rng.gauss(0.0, 1.0)
        v = (1.0 + c * x) ** 3
        if v <= 0.0:
            continue
        u = rng.random()
        x2 = x * x
        if u < 1.0 - 0.0331 * x2 * x2:
            return d * v
        if math.log(u) < 0.5 * x2 + d * (1.0 - v + math.log(v)):
            return d * v


def dirichlet_sample(alphas: List[float], rng: random.Random) -> List[float]:
    """Sample p ~ Dirichlet(alphas) via Gamma decomposition.  Falls back
    to a uniform vector if all alphas are tiny (degenerate)."""
    xs = [gamma_sample(max(a, 1e-12), rng) for a in alphas]
    s = sum(xs)
    if s <= 0.0:
        # Degenerate: return REF=1, ALTs=0
        return [1.0] + [0.0] * (len(alphas) - 1)
    return [x / s for x in xs]


def bottleneck_step(p_in: List[float], true_ne: float,
                    model: str, rng: random.Random) -> List[float]:
    """Apply one Wright-Fisher / Kimura bottleneck step to a category
    probability vector `p_in` (sums to 1, length = 1 + n_alts) and
    return the post-bottleneck category vector `p_out`.

    Discrete:   k ~ Multinomial(round(true_ne), p_in); p_out = k / true_ne.
                  This is the exact Wright-Fisher transition.
    Continuous: p_out ~ Dirichlet(p_in * (true_ne - 1)).
                  Kimura diffusion limit; for biallelic this reduces to
                  Beta(p_alt*(Ne-1), (1-p_alt)*(Ne-1)).

    Both have E[p_out|p_in] = p_in and the marginal variance of any
    category is p_in*(1-p_in)/Ne, matching one-generation drift.
    """
    if model == "discrete":
        ne_i = max(1, int(round(true_ne)))
        counts = [0] * len(p_in)
        for _ in range(ne_i):
            u = rng.random()
            acc = 0.0
            idx = len(p_in) - 1
            for ci, pc in enumerate(p_in):
                acc += pc
                if u <= acc:
                    idx = ci
                    break
            counts[idx] += 1
        return [c / ne_i for c in counts]
    else:
        # Continuous Dirichlet sampler.  For Ne <= 1 this is degenerate;
        # we already gate that in synthesize().
        ne1 = float(true_ne) - 1.0
        alphas = [c * ne1 for c in p_in]
        return dirichlet_sample(alphas, rng)


# ---------------------------------------------------------------------
# Synthesizer
# ---------------------------------------------------------------------

def synthesize(out_dir: Path,
               n_pairs: int,
               n_sites: int,
               true_ne: float,
               m_dp: int,
               c_dp: int,
               vaf_low: float,
               vaf_high: float,
               seed: int,
               n_multiallelic: int = 0,
               missing_gt_rate: float = 0.0,
               bottleneck_model: str = "discrete",
               outlier_frac: float = 0.0,
               n_generations: int = 1) -> None:
    rng = random.Random(seed)
    missing_gt_rate = max(0.0, min(1.0, missing_gt_rate))
    outlier_frac    = max(0.0, min(1.0, outlier_frac))
    n_generations   = max(1, int(n_generations))
    if bottleneck_model not in ("discrete", "continuous"):
        raise ValueError(f"Unknown bottleneck_model: {bottleneck_model}")
    # Continuous Beta-diffusion needs Ne > 1; discrete needs integer >= 1.
    if bottleneck_model == "discrete":
        true_ne_i = max(1, int(round(true_ne)))
    else:
        if true_ne <= 1.0:
            raise ValueError("continuous bottleneck requires true_ne > 1")

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

            # Bottleneck sampling (possibly iterated for multi-generation
            # transmission): produces `p_c_alts` (per-ALT child VAFs).
            categories = [p_ref] + p_alts          # length = 1 + n_alts
            p_curr = list(categories)
            for _g in range(n_generations):
                p_curr = bottleneck_step(p_curr, true_ne, bottleneck_model, rng)
            p_c_alts = p_curr[1:]   # drop REF, keep per-ALT

            # Optional outlier injection: with prob `outlier_frac`, replace
            # the post-bottleneck child VAFs with an extreme value (0 or 1).
            # This simulates NUMTs / mosaicism / errors that drive Kimura's
            # method-of-moments estimator upward.
            if outlier_frac > 0.0 and rng.random() < outlier_frac:
                if rng.random() < 0.5:
                    p_c_alts = [1.0] + [0.0] * (n_alts - 1)
                else:
                    p_c_alts = [0.0] * n_alts

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
        f"n_multiallelic={n_multiallelic} missing_gt_rate={missing_gt_rate} "
        f"bottleneck_model={bottleneck_model} outlier_frac={outlier_frac} "
        f"n_generations={n_generations}\n"
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
    p.add_argument("--true-ne",   type=float, default=5,
                   help="True bottleneck size used for simulation. Float "
                        "is allowed under --bottleneck-model continuous; "
                        "rounded to int under --bottleneck-model discrete.")
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
    p.add_argument("--bottleneck-model",
                   choices=["discrete", "continuous"], default="discrete",
                   help="Bottleneck process: `discrete` (Multinomial over "
                        "true_ne integer copies, matches `--model discrete` "
                        "MMLE) or `continuous` (Beta-diffusion / Dirichlet, "
                        "matches `--model continuous` MMLE; allows fractional "
                        "true_ne).")
    p.add_argument("--outlier-frac", type=float, default=0.0,
                   help="Fraction of (site, pair) cells whose post-bottleneck "
                        "child VAF is forced to 0 or 1, simulating NUMTs / "
                        "mosaicism / genotyping errors. Used to demonstrate "
                        "the outlier sensitivity of the Kimura estimator.")
    p.add_argument("--n-generations", type=int, default=1,
                   help="Number of bottleneck generations between mother and "
                        "child (1 = direct mother-child, 2 = grandmother to "
                        "grandchild, ...). Drift accumulates as "
                        "F_g = 1 - (1 - 1/Ne)^g, so feeding multi-generation "
                        "pairs into a single-generation estimator under-"
                        "estimates the true per-generation Ne.")
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
               missing_gt_rate=args.missing_gt_rate,
               bottleneck_model=args.bottleneck_model,
               outlier_frac=args.outlier_frac,
               n_generations=args.n_generations)
    print(f"Wrote cohort.vcf, cohort.fam, cohort.transmission_pairs.tsv "
          f"to {out_dir} (model={args.bottleneck_model}, true_ne={args.true_ne}, "
          f"n_pairs={args.n_pairs}, n_sites={args.n_sites}, "
          f"depth={args.m_dp}/{args.c_dp}, outliers={args.outlier_frac}, "
          f"n_generations={args.n_generations})")


if __name__ == "__main__":
    main()
