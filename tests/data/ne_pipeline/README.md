# `tests/data/ne_pipeline/` — end-to-end smoke-test data

This directory ships a small, fully synthetic mother-child mtDNA cohort
that lets you verify `mitoquest trans-prep` and `mitoquest ne-estimate`
run smoothly end-to-end without touching real data.

## Files

| File                              | Purpose                                                                                           |
| --------------------------------- | ------------------------------------------------------------------------------------------------- |
| `cohort.vcf`                      | Multi-sample VCF with 40 samples (20 mothers + 20 children) and 12 mtDNA SNV sites.               |
| `cohort.fam`                      | PLINK FAM file describing the 20 mother-child trios (fathers are placeholders).                   |
| `cohort.transmission_pairs.tsv`   | Pre-baked TSV in the exact schema emitted by `mitoquest trans-prep`, for testing `ne-estimate` alone. |
| `synthesize.py`                   | Deterministic regenerator (no numpy dependency). Default seed = `20260528`.                       |
| `run_demo.sh`                     | One-shot end-to-end pipeline: `trans-prep` → `ne-estimate` (with Kimura cross-check).             |

Simulation parameters baked into the default dataset:

* 20 mother-child pairs, 12 mtDNA sites, depth = 2,000 on both sides.
* **2 of the 12 sites are tri-allelic SNVs** (REF + 2 ALTs) to exercise
  the multi-allelic code path; the rest are biallelic.
* Maternal heteroplasmy `p_m ~ Uniform(0.15, 0.85)`.
* **True bottleneck `Ne = 5`** (single-generation transmission).

## Multi-allelic handling

`mitoquest trans-prep` decomposes each multi-allelic VCF record into one
output row **per ALT × per trio**.  For the default dataset this gives:

```
  10 biallelic sites  × 20 trios × 1 ALT  = 200 rows
   2 tri-allelic SNVs × 20 trios × 2 ALTs =  80 rows
   ----------------------------------------------
                                            280 rows
```

`ne-estimate` then treats each row as an independent
`(REF reads, ALT reads)` 2-allele observation — the standard biallelic
decomposition used by gnomAD-style estimators.  This is exact for any
single ALT and is a mild approximation only when two ALTs co-segregate
in the same mother (rare in mtDNA).

Mixed multi-allelic sites such as `A>G,GT` (one SNV + one indel) are
also handled: with `--snv-only` (default), the SNV ALT is kept and the
indel ALT is silently skipped within the same record.

## Quick start

Build `mitoquest` first (`cmake --build build -j8` from the repo root),
then:

```bash
bash tests/data/ne_pipeline/run_demo.sh
```

Expected stderr / stdout (numbers are deterministic):

```
==> trans-prep:  cohort.vcf + cohort.fam -> cohort.trans_prep.run.tsv
[trans-prep] Processed 12 variant records.
[trans-prep] Wrote 280 rows (280 PASS, 0 LOW_DEPTH).

==> ne-estimate (with --cross-check kimura) ...
[ne-estimate] Fitting Ne on 272 pairs (maternal VAF in [0.1, 0.9]).
[ne-estimate] Optimal Ne = 5 (95% CI: 5 - 5), max logL = -1371.16 on 272 pairs.
[ne-estimate] Kimura cross-check: b = 0.772, Ne_kimura = 4.38 on 272 informative pairs
```

`Ne = 5` recovers the simulated truth, and the Wonnapinij/Kimura cross-check
lands within ~15% of it as expected.

## Running the two subcommands by hand

```bash
# Step 1: extract M-C transmission pairs from the VCF + FAM
./bin/mitoquest trans-prep \
    -v tests/data/ne_pipeline/cohort.vcf \
    -f tests/data/ne_pipeline/cohort.fam \
    -d 100 \
    -o /tmp/cohort.trans_prep.tsv

# Step 2: fit Ne by Beta-Binomial MLE, with optional Kimura cross-check
./bin/mitoquest ne-estimate \
    -i /tmp/cohort.trans_prep.tsv \
    --min-ne 1 --max-ne 50 \
    --cross-check kimura \
    -t 4 \
    -o /tmp/cohort.ne.json
```

Or skip step 1 by feeding the pre-baked TSV directly to `ne-estimate`:

```bash
./bin/mitoquest ne-estimate \
    -i tests/data/ne_pipeline/cohort.transmission_pairs.tsv \
    --min-ne 1 --max-ne 50 \
    --cross-check kimura
```

## Regenerating with different settings

```bash
# Tighter bottleneck (true Ne = 2), 50 trios, fewer sites, fixed seed:
python3 tests/data/ne_pipeline/synthesize.py \
    --true-ne 2 --n-pairs 50 --n-sites 8 --seed 42

# Larger Ne (true Ne = 30):
python3 tests/data/ne_pipeline/synthesize.py --true-ne 30 --n-pairs 80

# Pure biallelic dataset (no multi-allelic sites):
python3 tests/data/ne_pipeline/synthesize.py --n-multiallelic 0

# Stress-test multi-allelic handling (every site tri-allelic):
python3 tests/data/ne_pipeline/synthesize.py --n-multiallelic 12
```

The script has no third-party dependencies (uses only the Python
standard library).

## Notes

* This dataset is **synthetic** and is intended only for verifying that
  the binaries run cleanly. It is not a benchmark of statistical
  performance — for that, see the `NeEstFindOptimal.RecoversTrueNeOnSyntheticCohort`
  test in [`tests/test_ne_estimate.cpp`](../../test_ne_estimate.cpp).
* The VCF declares `FORMAT/AD` with `Number=R` (one entry per allele,
  REF first), as required by `mitoquest trans-prep`.
* Positions are spread across `chrM` for plausibility but are not real
  mtDNA variants; do not use this VCF for any biological interpretation.
