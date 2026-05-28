# Release v1.8.4 — deCODE-style Bottleneck Bin-Simulation Output

**Date**: 2026-05-28

## Summary

v1.8.4 adds a new diagnostic output to `mitoquest ne-estimate` that reproduces
the per-VAF-bin bottleneck-parameter figure from the deCODE 2024 *Cell* paper
("Figure 5. Observed and simulated means for bottleneck parameter, b").  Two
new CLI flags emit a TSV containing observed and theoretically expected drift
per maternal-VAF bin, and a companion Python tool renders it as the familiar
two-panel figure.  This release also bundles the multi-generation Wright-Fisher
validation suite added since v1.8.3.

No breaking changes — all existing invocations and outputs are byte-identical.

## What's New

### 1. New CLI flags on `mitoquest ne-estimate`

- `--bin-simulation FILE` — when provided, writes a per-bin TSV alongside the
  JSON result. Each row contains observed drift (`obs_var`), sampling-corrected
  drift (`obs_var_corr`), per-bin `F_i = 1 - b` (`obs_F`, `obs_F_se`), and the
  theoretical curves at the fitted Ne and 95% CI.
- `--bin-simulation-bins INT` — number of equal-width VAF bins (default 10).

The TSV is fully self-describing: a commented header records the
`mitoquest_ne_estimate_command`, version, model, fitted Ne (with CI), Kimura
cross-check (when computed), and bin configuration, so the file can be
re-plotted without re-reading the source pair table.

Example:
```bash
mitoquest ne-estimate \
    -i cohort.transmission_pairs.tsv \
    --cross-check kimura --kimura-bootstrap 200 \
    --bin-simulation cohort.bin_sim.tsv --bin-simulation-bins 10 \
    -o cohort.ne.json
```

### 2. New Python tool: `tools/plot_bottleneck_simulation.py`

Renders the deCODE-style two-panel figure directly from the bin-simulation TSV:

- **Left panel** — observed `(p_c − p_m)²` per bin vs `p_m`, overlaid with the
  parabolic curve `p_m(1 − p_m) / Ne` at `Ne_MLE` (red, with 95% CI ribbon)
  and at `Ne_Kimura` (blue dashed).  Marker size is proportional to the number
  of pairs per bin.
- **Right panel** — observed per-bin `F_i = 1 − b` vs `p_m`, overlaid with
  horizontal lines at `1/Ne_MLE` (red, with 95% CI band) and `1/Ne_Kimura`
  (blue dashed).

```bash
python tools/plot_bottleneck_simulation.py \
    -i cohort.bin_sim.tsv \
    -o cohort.bottleneck.png
```

The script follows the same coding style as the other `tools/` Python utilities
(argparse with `ArgumentDefaultsHelpFormatter`, Agg backend for headless use,
`parse_args() → argparse.Namespace`, `main() → None`).

### 3. Multi-generation Wright-Fisher validation suite

Added under `tests/data/ne_pipeline/`:

- `synthesize.py` — extended to simulate `g > 1` generations of Wright-Fisher
  drift between mother and child.
- `verify_popgen.py` — population-genetics sanity checks (drift variance,
  fixation/loss rates, segregating sites) across generations.
- `validate_estimators.py` — head-to-head accuracy comparison of Continuous
  MLE vs Kimura vs Discrete MLE across generation counts and Ne values.
- `README.md` — adds a "When discrete MLE is actually the right choice"
  section documenting the regimes where each estimator wins.

## Files Changed

- `src/ne_estimate.h` — adds `Config::bin_simulation_file`,
  `Config::bin_simulation_n_bins`, `BinSimulationRow` struct, and
  `compute_bin_simulation()` static helper declaration.
- `src/ne_estimate.cpp` — implements `compute_bin_simulation()` (Welford-style
  per-bin accumulation reusing `prepare_pair_contributions`), CLI parsing,
  `usage()` text, and emission block in `run()` after the Kimura cross-check.
- `tools/plot_bottleneck_simulation.py` — new (229 lines).
- `tools/README.md` — new section documenting purpose, invocation, and
  interpretation of the figure.
- `tests/data/ne_pipeline/synthesize.py` — multi-generation extension.
- `tests/data/ne_pipeline/verify_popgen.py` — new validation script.
- `tests/data/ne_pipeline/validate_estimators.py` — new comparison script.
- `tests/data/ne_pipeline/README.md` — guidance on estimator selection.
- `CMakeLists.txt` — version bump 1.8.3 → 1.8.4.

## Verification

- All 22 `NeEst*` GoogleTest cases still pass.
- End-to-end smoke test on `tests/data/ne_pipeline/cohort.transmission_pairs.tsv`
  (236 pairs) reports `Ne_MLE = 2.61` (95% CI 2.36 – 2.89),
  `Ne_Kimura = 4.71`, with 8 bins emitted to TSV and the corresponding figure
  rendered correctly.

## Upgrade Notes

- No breaking changes.  Existing JSON output, CLI flags, and on-disk formats
  are unchanged.  The two new flags are strictly additive.
