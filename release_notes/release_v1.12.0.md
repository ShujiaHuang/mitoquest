# mitoquest v1.12.0 — ne-estimate: marginalized maternal frequency by default

## Highlights

This is a **minor release** (one feature commit since v1.11.4). It promotes
**maternal-frequency marginalization to the default continuous estimator** in
`mitoquest ne-estimate`: the mother's latent allele frequency `p_m` is now
integrated out against its Beta read-sampling posterior instead of being
plugged in as the point estimate `p_m = m_alt / m_dp`. It also adds
`--min-depth`, per-family summary statistics, and depth/QC diagnostics to the
JSON / TSV / stderr outputs. The new options and output fields are
**backward-compatible additions**; only the default estimator's numerical
behaviour changes, and that is opt-outable via `--no-maternal-marginalization`.
Version bump `1.11.4 → 1.12.0`.

---

## Behavior Change

- **Default continuous-model estimator**: the maternal frequency `p_m` is now
  marginalized (integrated out) **by default**. The previous default used the
  plug-in point estimate `p_m = m_alt / m_dp`, which folds the mother's own
  read noise into the drift signal and makes the fitted `Ne` depth-dependent.
  Marginalizing `p_m` removes that read-noise term, so the reported `Ne` is far
  less sensitive to maternal/child sequencing depth.
- Opt out with **`--no-maternal-marginalization`** to reproduce the legacy
  plug-in path exactly.
- The legacy `--model discrete` path already integrated a maternal Beta
  posterior and is unaffected.
- The mode actually used is echoed in the output as
  `Maternal_Marginalization: true | false`.

---

## New Functionality

### Command line

- **`--min-depth INT`** — drop mother–child pairs whose **mother or child**
  depth is below the threshold (applied to both sides; a shallow child gives an
  equally noisy estimate). Default `0` disables the filter; values `< 0` are
  rejected. Echoed back as `Min_Depth` and `Min_Depth_Skipped`.

### Output (all additive)

- Top-level JSON depth descriptors: `Mean_Mother_DP`, `Mean_Child_DP`,
  `Harmonic_Mother_DP`, `Harmonic_Child_DP` — descriptive cohort-comparability
  diagnostics (no aggregate of them is a plug-in back-correction handle; see
  `ne_estimate.h`).
- Trio load accounting surfaced in JSON: `Trio_Founder_Mismatch_Skipped`,
  `Trio_Founder_Hom_Skipped`.
- With `--per-family`, a new **`Per_Family_Summary`** JSON object
  (`N_Families_Estimated`, `N_Families_Skipped`, `N_Families_CI_Clipped`,
  `Mean_Ne`, `Median_Ne`) alongside the existing `Per_Family_Estimates` array.
  The per-family TSV gains the header keys `#n_families_estimated`,
  `#n_families_skipped`, `#n_families_ci_clipped`, `#mean_family_ne`,
  `#median_family_ne`, and stderr now reports the mean + median and warns when
  any family's CI is clipped at the search boundary (a shallow-maternal-depth
  artifact of the marginalized fit — see the caveats in `README.md`).

---

## Bug Fixes

- **Per-family depth precision**: `MEAN_MOTHER_DP` / `MEAN_CHILD_DP` and the
  marginal-log-likelihood columns in the per-family TSV were written with
  `setprecision(2)`; under `defaultfloat` this silently dropped digits (e.g.
  `674.75` rendered as `6.7e+02`). Raised to `setprecision(8)`.
- **`-Wl,-no_compact_unwind` removed from the production binary.** The flag
  disabled the unwind tables the C++ runtime needs to propagate an exception
  across translation units, so a user-facing `std::runtime_error` (e.g. an
  invalid `--min-depth`) reached `__cxa_throw` with no reachable handler and
  aborted with **SIGABRT (exit 134)** instead of printing `Error: ...` and
  exiting `1`. This had affected error handling for *every* shipped subcommand;
  scoping the flag to the test-only binary was the inverted fix. Dropping it
  emits no warnings on Apple clang / arm64.
- Corrected the direction of the `--model discrete` plug-in warning text.

---

## Tests

- `ctest` **3/3** green on macOS / arm64 (Apple clang), including
  `TrioLikelihoodMpmathCrossCheck` (runs against mpmath — not silently skipped),
  across **186** GoogleTest cases.

---

## Documentation

- `README.md`: refreshed the `ne-estimate` help and the prior-mismatch /
  depth-bias sections to describe the marginalized default, the `--min-depth`
  lever, and the measured gate-isolation result (the VAF window contributes an
  `Ne`-independent factor and does not move the MLE).
- `docs/INSTALL_FROM_SOURCE.md`: documented the ctest mpmath-skip pitfall —
  `TrioLikelihoodMpmathCrossCheck` exits `77` (rendered `Skipped`, not a
  failure) when `mpmath` is not importable for the discovered interpreter — and
  how to point CMake at an mpmath-enabled Python
  (`-DPython3_EXECUTABLE=/path/to/python-with-mpmath`).

---

## Files changed (main)

| Area | Files |
|---|---|
| ne-estimate | `src/ne_estimate.*`, `tests/test_ne_estimate.cpp` |
| build / version | `CMakeLists.txt` (project VERSION 1.11.4 → 1.12.0; drop `-Wl,-no_compact_unwind`) |
| docs | `README.md`, `docs/INSTALL_FROM_SOURCE.md`, `release_notes/release_v1.12.0.md` |
