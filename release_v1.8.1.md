# MitoQuest v1.8.1

## Highlights

Bug-fix and feature release for `mitoquest ne-estimate`:

- **Fix** the Wonnapinij sampling-noise correction in the Kimura cross-check
  (denominators `m_dp` / `c_dp` were swapped).
- **Add** a non-parametric pair-level bootstrap **95% CI** for the Kimura
  cross-check (`b` and `Ne_Kimura`).
- **Document** the expected real-data divergence between the Beta-Binomial
  MLE and the Kimura cross-check.

## What changed

### Kimura sampling-correction bug fix

The Wonnapinij 2010 sampling-noise correction subtracts the binomial
sampling variance contributed by the mother and the child read counts:

```
s_i  =  p_m (1 - p_m) / m_dp   +   p_c (1 - p_c) / c_dp
```

Each variance term must be divided by **the depth at which that side was
sampled**.  v1.8.0 had the depths swapped (`p_m * (1-p_m) / c_dp` and
`p_c * (1-p_c) / m_dp`).  The bug is essentially invisible when mother
and child are sequenced at similar depth (the typical mtDNA cohort
configuration), but is now corrected to match Wonnapinij 2010.

### Non-parametric bootstrap CI for the Kimura cross-check

Two new CLI options on `mitoquest ne-estimate` (active only when
`--cross-check kimura` is set):

| Flag | Default | Description |
|------|---------|-------------|
| `--kimura-bootstrap N` | `1000` | Pair-level bootstrap iterations.  `0` disables the CI. |
| `--kimura-seed K`       | `42`   | RNG seed for the bootstrap (reproducible). |

When `N > 0` the JSON output now includes:

```json
"Kimura_Cross_Check": {
  "b":                    0.787734,
  "Ne_Kimura":            4.7110793,
  "b_CI_95_Low":          0.75497181,
  "b_CI_95_High":         0.81970617,
  "Ne_Kimura_CI_95_Low":  4.0811631,
  "Ne_Kimura_CI_95_High": 5.5465015,
  "N_Bootstrap":          1000,
  "Bootstrap_Seed":       42,
  ...
}
```

The CI is computed by resampling the informative pairs with replacement,
recomputing `b` per resample, and taking the 2.5 / 97.5 percentiles.
Infinite `Ne_Kimura` values (which arise when a bootstrap resample
produces `V <= 0`) are emitted as the JSON string `"Infinity"` so
downstream consumers can detect them unambiguously.

### MLE vs. Kimura: when they diverge

On synthetic data drawn from the model, both estimators agree (e.g., on
the included demo with true Ne = 5: `Ne_MLE = 5`, `Ne_Kimura = 4.71
[4.08, 5.55]`).

On real cohorts the two can diverge significantly — for example
`Ne_MLE = 23` vs `Ne_Kimura = 2.78`.  This is **not** a software defect;
the two estimators answer different questions:

- The **Beta-Binomial MLE** is the exact discrete-Wright-Fisher
  likelihood for a single transmission.  It is sensitive to the discrete
  grid resolution `k / Ne`: hundreds of low-drift heteroplasmic pairs
  (the typical mtDNA cohort) require `Ne` large enough that `k / Ne` can
  land near the observed concordant VAFs.
- The **Wonnapinij `b`** is a method-of-moments estimator that uses only
  the cohort-level variance of `(p_c - p_m)`.  A small number of
  high-drift outliers (errors / NUMTs / mixed populations) inflate that
  variance and pull `Ne_Kimura` downward.

A wide gap therefore signals heterogeneous transmission or contamination,
not an algorithm defect.  The Beta-Binomial MLE remains the primary
reported estimate; the Kimura cross-check (now with CI) is provided for
comparison with the deCODE 2024 *Cell* paper and as a sanity check.

## Tests

- All 22 `NeEst*` and `TransPrep*` tests still pass.
- Synthetic demo (`tests/data/ne_pipeline/run_demo.sh`) recovers
  `Ne_MLE = 5` and `Ne_Kimura = 4.71 [4.08, 5.55]` (true Ne = 5).

## Compatibility

- No CLI changes to existing flags; new flags are opt-in with backwards-
  compatible defaults.
- The JSON shape under `Kimura_Cross_Check` adds new keys
  (`b_CI_95_Low`, `b_CI_95_High`, `Ne_Kimura_CI_95_Low`,
  `Ne_Kimura_CI_95_High`, `N_Bootstrap`, `Bootstrap_Seed`).  Older
  consumers reading `b`, `Ne_Kimura`, `N_Informative`, `Note`, or
  `Method` are unaffected.
