# Release v1.8.6 — Rename MLE → MMLE (Maximum Marginal Likelihood Estimator)

**Date**: 2026-05-30

## Summary

v1.8.6 is a **terminology / naming refinement** release.  No algorithm,
no statistical model, and no numerical output value has changed — only
the *labels* used in code, console output, JSON output, TSV column
headers, on-disk metadata, and documentation.

The previous label "MLE" / "Likelihood" was technically loose: the
estimator that `mitoquest ne-estimate` reports is **not** a plain
maximum-likelihood fit of the joint observed-data likelihood.  It is a
**Maximum Marginal Likelihood Estimator (MMLE)**:

1. The latent maternal / child true allele frequencies are
   **analytically integrated out** (Beta-Binomial conjugacy / Beta-
   diffusion), yielding a per-pair *marginal* likelihood that depends
   only on Ne and the observed counts.
2. Mother-child pairs are treated as independent (composite /
   pseudo-likelihood), so the global objective is the sum of per-pair
   *marginal* log-likelihoods.
3. The estimator maximises this composite marginal likelihood over Ne.

Calling that "MLE" obscured both assumptions (the marginalisation and
the composite-likelihood independence assumption).  The new label
"MMLE" makes those assumptions explicit and is consistent with the
hierarchical-model literature.

## What Changed

### 1. CLI / `usage()` text

`mitoquest ne-estimate -h` now describes the estimator explicitly:

```
  Estimate the mitochondrial DNA bottleneck size (Ne) from mother-child
  transmission pairs using the Maximum Marginal Likelihood Estimator
  (MMLE).  The latent mother / child true allele frequencies are
  analytically integrated out (Beta-Binomial conjugacy), and pairs are
  treated as independent (composite likelihood), yielding a per-pair
  marginal log-likelihood that is maximised jointly over Ne.
```

`mitoquest -h` likewise advertises `ne-estimate` via "Beta-Binomial
MMLE" instead of "Beta-Binomial MLE".

### 2. JSON output (**breaking**)

| Old key                  | New key                          |
|--------------------------|----------------------------------|
| `Max_LogLik`             | `Max_Marginal_LogLik`            |
| *(none)*                 | `Estimator: "MMLE (composite marginal likelihood)"` (new) |

The numerical value of `Max_Marginal_LogLik` is bit-for-bit identical
to the v1.8.5 `Max_LogLik` value for the same input.

### 3. `--ne-profile` TSV (**breaking**)

Column header rename:

| Old column     | New column      |
|----------------|-----------------|
| `mle_log_lik`  | `mmle_log_lik`  |
| `mle_delta_2ll`| `mmle_delta_2ll`|

Metadata header line rename:

| Old metadata key            | New metadata key               |
|-----------------------------|--------------------------------|
| `#fitted_ne_mle=`           | `#fitted_ne_mmle=`             |
| `#fitted_ne_mle_ci_low=`    | `#fitted_ne_mmle_ci_low=`      |
| `#fitted_ne_mle_ci_high=`   | `#fitted_ne_mmle_ci_high=`     |
| `#max_log_lik=`             | `#max_marginal_log_lik=`       |
| `#best_ne_mle_on_grid=`     | `#best_ne_mmle_on_grid=`       |

### 4. `--bin-simulation` TSV metadata (**breaking**)

`#max_log_lik=` → `#max_marginal_log_lik=`.

### 5. Console output

```
[ne-estimate]   Best Ne on grid: MMLE = 2.6, Kimura SSR = 4.8 ...
```

(was `MLE = 2.6`).

### 6. C++ public API (**source-breaking for downstream code**)

In `src/ne_estimate.h`, struct `NeProfileRow`:

| Old field      | New field        |
|----------------|------------------|
| `mle_log_lik`  | `mmle_log_lik`   |
| `mle_delta_2ll`| `mmle_delta_2ll` |

The internal struct `NeEstimator::Result.max_log_lik` is **kept**
unchanged for source stability; only the public-facing labels move.

### 7. Python tools

- `tools/plot_ne_profile.py`:
  - The function previously named `plot_mle_panel` is now
    `plot_mmle_panel`.
  - `load_ne_profile()` accepts both the new (`mmle_*`) and old
    (`mle_*`) column names, and both new (`#fitted_ne_mmle*`,
    `#max_marginal_log_lik`) and old (`#fitted_ne_mle*`,
    `#max_log_lik`) metadata keys.  v1.8.5 TSVs still plot correctly.
- `tools/plot_bottleneck_simulation.py` — legend label `Ne_MLE` →
  `Ne_MMLE`.
- `tools/plot_mito_bottlenechk_smulation.py` — all in-figure labels and
  log messages renamed.
- `tools/README.md` — TSV column documentation updated.

### 8. Tests

- `tests/test_ne_estimate.cpp` — comments updated; one test renamed
  `TrimmedCloserToMLEWhenOutliersPresent` →
  `TrimmedCloserToMMLEWhenOutliersPresent` (the assertion logic is
  identical).
- `tests/data/ne_pipeline/scenarios/**/*.json` and
  `tests/data/ne_pipeline/cohort.ne.run.json` — committed reference
  fixtures regenerated to use `Max_Marginal_LogLik` (23 files).
- `tests/data/ne_pipeline/validate_estimators.py` — reads
  `Max_Marginal_LogLik`, falling back to `Max_LogLik` for
  backward compatibility with v1.8.5-generated JSONs.
- `tests/data/ne_pipeline/synthesize.py`,
  `tests/data/ne_pipeline/ne_diagnostic.py` — diagnostic prints
  updated.

### 9. Documentation

- `README.md` — section heading renamed to "Bottleneck size (Ne)
  Maximum Marginal Likelihood Estimation (MMLE)"; introductory
  paragraph rewritten to explain the marginal-likelihood and
  composite-likelihood assumptions; "Why two MMLE models?" / "Why
  Ne_MMLE is often smaller than Ne_Kimura" / decision-rule subsections
  consistently relabelled; JSON example updated.
- `tools/README.md` — all MLE references relabelled MMLE.

## Files Changed

- `CMakeLists.txt` — version bump 1.8.5 → 1.8.6.
- `src/main.cpp` — top-level `usage()`.
- `src/ne_estimate.h` — file-header doc + `NeProfileRow` field rename.
- `src/ne_estimate.cpp` — comments, `usage()`, JSON output, TSV
  emission, metadata, console messages, warnings, local variables.
- `tests/test_ne_estimate.cpp` — comments + one test name.
- `tests/data/ne_pipeline/{validate_estimators,synthesize,ne_diagnostic}.py`.
- `tests/data/ne_pipeline/scenarios/**/*.json`,
  `tests/data/ne_pipeline/cohort.ne.run.json` — bulk JSON key rename.
- `tools/plot_ne_profile.py`,
  `tools/plot_bottleneck_simulation.py`,
  `tools/plot_mito_bottlenechk_smulation.py`,
  `tools/README.md`.
- `README.md` — `ne-estimate` section rewritten with MMLE terminology.

## Verification

End-to-end smoke test on
`tests/data/ne_pipeline/cohort.transmission_pairs.tsv` (236 pairs,
`--max-ne 30`, `--ne-profile-step 0.1`):

- `Optimal Ne`, `95% CI`, `Max_Marginal_LogLik`, `Kimura b`,
  `Ne_kimura` all bit-for-bit identical to v1.8.5.
- All 22 `NeEst*` GoogleTest cases pass.
- `tools/plot_ne_profile.py` renders new-format TSV identically to
  v1.8.5 output, and still parses v1.8.5 TSVs via the
  backward-compatibility shim.

## Upgrade Notes (Breaking Changes)

Downstream consumers of `mitoquest ne-estimate` output must update:

- **JSON readers**: `Max_LogLik` → `Max_Marginal_LogLik`.  A safe
  one-liner shim is:
  ```python
  max_loglik = j.get("Max_Marginal_LogLik", j.get("Max_LogLik"))
  ```
- **TSV readers** of `--ne-profile` output: rename column accesses
  `mle_log_lik` → `mmle_log_lik`, `mle_delta_2ll` → `mmle_delta_2ll`.
  The bundled `tools/plot_ne_profile.py` already auto-detects either
  schema.
- **Awk/grep parsers** of TSV metadata header lines: rename keys per
  the table in §3 above.
- **Downstream C++ code** that linked against `libne_estimate` and
  read `NeProfileRow::mle_log_lik` must rename to
  `NeProfileRow::mmle_log_lik`.

There are **no other behaviour changes** in v1.8.6.  Anyone willing to
treat the new keys as new field names (and the old ones as deprecated)
needs no further code changes.
