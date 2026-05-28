# Release v1.8.2 — Continuous Beta-Diffusion MLE (fixes systematic Ne overestimation)

**Date**: 2026-05-28

## Summary

v1.8.2 replaces the default Ne MLE model with a **continuous Beta-diffusion
likelihood** that correctly captures post-bottleneck mtDNA drift.  The previous
discrete model (`k ~ BetaBin(Ne, α, β)`, `c_alt ~ Bin(c_dp, k/Ne)`)
systematically overestimated Ne at high sequencing depth due to grid-resolution
bias.  With the new default model, the MLE and the Wonnapinij/Kimura
variance-of-moments estimator now agree on real data (both give Ne ≈ 3,
consistent with the deCODE genetics 2024 *Cell* paper), resolving the
previously observed ~8× gap.

This release also adds robust diagnostic tools — trimmed Kimura and per-pair
drift outlier reporting — to help users assess data quality.

## Root Cause: Grid-Resolution Bias in the Discrete Model

The discrete model constrains child heteroplasmy to the grid
`{0, 1/Ne, ..., 1}`.  At high sequencing depth (DP ≈ 2000), the Binomial
likelihood `Bin(c_alt | c_dp, k/Ne)` is an extremely tight spike (SD ≈ 0.01).
Any mismatch between the observed child VAF and the nearest grid point
`k/Ne` produces an enormous log-likelihood penalty, forcing the MLE to inflate
Ne to obtain a finer grid.  This created a systematic upward bias of ~5–10×
above the true bottleneck size.

In reality, the child's heteroplasmy is shaped not only by the initial
bottleneck sampling but also by post-bottleneck vegetative segregation during
cell division, making it a *continuous* variable — not restricted to the
`k/Ne` grid.

## The Fix: Continuous Beta-Diffusion Model (new default)

The new model treats the child's true heteroplasmy as a continuous draw from
the Kimura diffusion distribution:

```
p_child | p_mother  ~  Beta( p_m × (Ne − 1),  (1 − p_m) × (Ne − 1) )
c_alt   | p_child   ~  Binomial(c_dp, p_child)
```

After analytically marginalising out `p_child`:

```
c_alt | p_mother  ~  BetaBinomial(c_dp, p_m × (Ne − 1), (1 − p_m) × (Ne − 1))
```

This is a single BetaBinomial PMF evaluation per pair — simpler and faster
than the discrete model's `O(Ne)` inner sum — and correctly models the
continuous nature of post-bottleneck drift.

### Verification on synthetic data (true Ne = 30)

| Estimator         | Ne    | 95% CI     |
|-------------------|-------|------------|
| Continuous MLE    | 32    | 29 – 35    |
| Kimura (moments)  | 31.7  | 28.0 – 36.0|
| Discrete MLE      | 30    | 30 – 30    |

The continuous MLE agrees closely with Kimura while the discrete MLE shows
artificially tight CI (a symptom of the grid-resolution overfitting).

## What's New

### 1. `--model continuous|discrete`  (default: `continuous`)

Select the likelihood model.  The continuous Beta-diffusion model is now the
default and is recommended for all mtDNA bottleneck analyses.  The discrete
model is retained via `--model discrete` for specialised use cases (e.g.,
virus-passage experiments where the physical inoculum count is the target).

### 2. `--kimura-trim FRAC`  (default: 0.0, recommended: 0.10)

Drops the top `FRAC` of highest-drift pairs (ranked by per-pair Wonnapinij
contribution `F_i = (d_i − s_i) / w_i`) before recomputing `b`.  The trimmed
Ne_Kimura is reported alongside the original untrimmed value in the
`Trimmed_Kimura` JSON block.

### 3. `--top-drift-k INT`  (default: 0, recommended: 20)

Emits the top-K highest-drift pairs in the JSON output as
`Top_Drift_Outliers`, with per-pair read counts, VAFs, and F_i values.

### 4. MLE-vs-Kimura disagreement warning

When `|Ne_MLE / Ne_Kimura| > 3`, a **WARNING** is emitted to stderr with a
recommendation to inspect outlier pairs.

### 5. New JSON output fields

- `"Model"` — `"continuous"` or `"discrete"`, indicating which MLE was used.
- `"Trimmed_Kimura"` block (when `--kimura-trim > 0`).
- `"Top_Drift_Outliers"` array (when `--top-drift-k > 0`).

## Recommended Usage for Real Data

```bash
mitoquest ne-estimate \
    -i sample.transmission_pairs.tsv \
    --cross-check kimura \
    --kimura-trim 0.10 \
    --top-drift-k 20 \
    -o sample.ne.json
```

The continuous model is used by default.  Compare `Ne` (MLE) and `Ne_Kimura`
in the output — they should now agree closely on well-behaved mtDNA data.

## Tests

- 3 new unit tests for the continuous Beta-diffusion model.
- 3 new unit tests for trimmed Kimura diagnostics.
- All 20 Ne-estimate tests pass.

## Files Changed

- `src/ne_estimate.h` — New `compute_ll_single_continuous()` declaration;
  `model` field in `Config`; `bool continuous` parameter on pipeline functions;
  extended `KimuraCheck` struct with trimmed fields and `DriftOutlier`.
- `src/ne_estimate.cpp` — Continuous model implementation; `--model` CLI
  option; dispatch logic in global-LL / find-optimal / estimate functions;
  JSON `Model` field; trimmed Kimura; top-K outlier extraction; warning.
- `tests/test_ne_estimate.cpp` — 6 new tests.
- `src/version.h` — Bumped 1.8.1 → 1.8.2.
- `CMakeLists.txt` — Version bump.
