# Release v1.8.5 — deCODE-style Ne-Profile Diagnostic (MLE vs Kimura)

**Date**: 2026-05-30

## Summary

v1.8.5 reproduces the deCODE 2024 *Cell* paper's "best-fit Ne" exercise
inside `mitoquest ne-estimate`.  The deCODE authors picked Ne ≈ 3 by
scanning candidate Ne values, simulating the Kimura allele-frequency-
change distribution at each one, and choosing the Ne that best fits the
observed distribution across 137 variants × 53,041 mother-child pairs.

This release adds an analogous Ne-profile diagnostic that scores **every**
candidate Ne under both estimators in the program (MLE and Kimura) on the
*same* informative pair set, so the user can see directly whether the
two estimators agree on the location of the best Ne — and, when they
don't, which one each piece of the data actually supports.

No breaking changes — the two new flags are strictly additive.

## What's New

### 1. New CLI flags on `mitoquest ne-estimate`

- `--ne-profile FILE` — when provided, writes a 5-column TSV scoring every
  candidate Ne in `[--min-ne, --max-ne]` under two metrics:
  - `mle_log_lik(Ne)` — global log-likelihood under the configured model
    (continuous Beta-diffusion or discrete Beta-Binomial).  Maximised at
    the fitted `Ne_MLE`.
  - `kimura_ssr(Ne)` = Σᵢ ((dᵢ − sᵢ) − p_mᵢ (1 − p_mᵢ) / Ne)².
    Per-pair least-squares fit of the one-generation Wright-Fisher
    prediction.  Minimised at the closed-form
    `Ne_Kimura_SSR = Σ w² / Σ rw`.
- `--ne-profile-step FLOAT` — grid step on the Ne axis (default 0.1).

The TSV is fully self-describing: a commented header records the command
line, version, model, fitted `Ne_MLE` (with 95% profile CI), Wonnapinij
`b` and `Ne_Kimura` (with bootstrap CI when enabled), the analytic
`Ne_Kimura_SSR`, and grid-best Ne under each metric, so the file can be
re-plotted without re-reading the source pair table.

Example:
```bash
mitoquest ne-estimate \
    -i cohort.transmission_pairs.tsv \
    --cross-check kimura --kimura-bootstrap 200 \
    --ne-profile cohort.ne_profile.tsv --ne-profile-step 0.1 \
    --max-ne 30 \
    -o cohort.ne.json
```

### 2. New Python tool: `tools/plot_ne_profile.py`

Renders the Ne-profile TSV as a two-panel figure:

- **Left panel** — `−2(logL − logL_max)` vs Ne, with the χ²(1, 0.95) =
  3.841 threshold marked, the fitted `Ne_MLE` annotated, and the 95%
  profile CI shaded.
- **Right panel** — Kimura per-pair SSR (normalised by its minimum) vs
  Ne, with the analytic `Ne_Kimura_SSR` and the Wonnapinij/method-of-
  moments `Ne_Kimura = 1 / (1 − b)` annotated, the bootstrap 95% CI
  shaded, and the deCODE Ne = 3 reference marked.

```bash
python tools/plot_ne_profile.py \
    -i cohort.ne_profile.tsv \
    -o cohort.ne_profile.png
```

The tool follows the same coding style as the other `tools/` Python
utilities (argparse with `ArgumentDefaultsHelpFormatter`, Agg backend
for headless use, `parse_args() → argparse.Namespace`,
`main() → None`).

### 3. New public C++ helpers

Two new static methods on `NeEstimator` (declared in `src/ne_estimate.h`,
implemented in `src/ne_estimate.cpp`):

- `compute_ne_profile(data, lf, min_ne, max_ne, step, threads, continuous)`
  — performs the full grid scan and returns a vector of `NeProfileRow`
  with both raw (`mle_log_lik`, `kimura_ssr`) and normalised
  (`mle_delta_2ll`, `kimura_norm_ssr`) metrics.
- `kimura_ssr_best_ne(data)` — closed-form Kimura-SSR best-fit Ne;
  returns `NaN` when `Σ rw ≤ 0`.

## Files Changed

- `src/ne_estimate.h` — adds `Config::ne_profile_file`,
  `Config::ne_profile_step`, `NeProfileRow` struct, and declarations of
  `compute_ne_profile()` / `kimura_ssr_best_ne()`.
- `src/ne_estimate.cpp` — implements both helpers, CLI parsing,
  `usage()` text, and emission block in `run()`.
- `tools/plot_ne_profile.py` — new (227 lines).
- `tools/README.md` — new section documenting purpose, invocation, and
  interpretation of the figure.
- `CMakeLists.txt` — version bump 1.8.4 → 1.8.5.

## Verification

End-to-end smoke test on
`tests/data/ne_pipeline/cohort.transmission_pairs.tsv` (236 pairs,
`--max-ne 30`, `--ne-profile-step 0.1`, 291 grid points):

```
[ne-estimate]   Best Ne on grid: MLE = 2.6, Kimura SSR = 4.8 (analytic = 4.77744)
[ne-estimate] Optimal Ne = 2.60917 (95% CI: 2.35917 - 2.88917)
[ne-estimate] Kimura cross-check: b = 0.787734, Ne_kimura = 4.71108
              [95% CI: Ne_kimura 4.06607 - 5.41786 via 200 bootstraps]
```

This is exactly the diagnostic the user wanted: the MLE and Kimura
estimators disagree by ~2× on this synthetic cohort, the figure shows
the disagreement clearly, and the `--kimura-trim 0.10 --top-drift-k 20`
recipe (already present since v1.8.x) explains what to do about it.

All 22 `NeEst*` GoogleTest cases still pass.

## Upgrade Notes

- No breaking changes.  Existing JSON output, CLI flags, and on-disk
  formats are unchanged.  The two new flags are strictly additive.
- The Ne = 1 grid point is reported as `mle_log_lik = -inf` under the
  continuous model (degenerate point: complete drift to fixation).  The
  plotting tool drops these rows automatically.
