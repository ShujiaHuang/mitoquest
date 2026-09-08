# mitoquest v1.12.1 — ne-estimate: opt-in deCODE-style binned Kimura estimator

## Highlights

This is a **patch release** (one feature commit since v1.12.0). It adds an
**opt-in, default-off** comparability report to `mitoquest ne-estimate`:
`--kimura-decode-style` emits a binned method-of-moments `Ne` that reproduces
the deCODE pedigree estimator (Árnadóttir et al., *Cell* 2024, Table 2) so the
mitoquest `Ne` can be placed on the same numerical footing as the published
deCODE values. The change is **purely additive to the report layer** — the
primary continuous MMLE `Ne`, its CI, the Kimura cross-check, and every existing
JSON field are untouched. With the flag off (the default) the JSON output is
**byte-for-byte identical** to v1.12.0. Version bump `1.12.0 → 1.12.1`.

Because the estimator is opt-in, non-breaking, and changes no default behaviour,
it ships as a patch increment following the project's convention for
backward-compatible report-layer additions.

---

## New Functionality

### Command line

- **`--kimura-decode-style`** — *also* emit a deCODE-style binned Kimura moment
  estimator. Independent of `--cross-check kimura`; default **off**, never
  implied by any other flag.

### Output (all additive)

- A new top-level **`Decode_Style_Kimura`** JSON object, present only when the
  flag is set: `Ne_Decode_Style`, `g`, `Bin_Width`, `Min_VAF`, `Max_VAF`,
  `N_Bins`, `N_Bins_Valid`, `N_Pairs_Used`, a per-bin `Bins[]` array
  (`Lo`, `Hi`, `N`, `Pbar_M`, `Hbar`, `Var_H`, `Denom`, `b`, `Ne_Bin`,
  `Status`), plus `Note` and `Method`. Non-finite values serialise as the JSON
  strings `"NaN"` / `"Infinity"` / `"-Infinity"`, consistent with the existing
  emitter.
- A one-line stderr summary of the binned `Ne`, valid-bin count, and pairs used.

---

## Method

Mothers are grouped by their plug-in read-frequency `p_hat_M = m_alt / m_dp`
into `5%` VAF bins spanning `[--min-vaf, --max-vaf]` (default `[0.10, 0.90]` →
exactly deCODE's 16 intervals). Within each bin, for a single mother→child
transmission (`g = 1`):

```
b_bin   = 1 - var(h) / [ pbar_M (1 - pbar_M) ]
Ne_bin  = -g / ln(b_bin)                      (diffusion convention)
```

`var(h)` is the **uncorrected** unbiased sample variance (`ddof = 1`) of the
child read-frequencies about the bin's mean child frequency `hbar`. The overall
`Ne` is the **child-count-weighted** mean of `Ne_bin` over all bins that yield a
finite `Ne` (i.e. `b ∈ (0, 1)`).

This is a **comparability-only** report. It deliberately omits the Wonnapinij
read-sampling-noise correction that the default `--cross-check kimura` applies,
and it uses the diffusion inversion `-g/ln b` rather than the discrete
Wright–Fisher `1/(1-b)`; at conventional mtDNA depth it is therefore
statistically **worse** than the corrected cross-check. It exists so the number
can be compared like-for-like with deCODE's published Table-2 values, not to
replace the default estimator. Bins that are empty, single-pair, or yield
`b ≤ 0` / `b ≥ 1` are reported with an explicit `Status` and excluded from the
weighted mean.

---

## Validation

- **Three-way numerical agreement** on a 3 571-pair simulated cohort
  (`true Ne = 3`), all computed independently:

  | Path | `Ne_Decode_Style` |
  |---|---|
  | C++ (JSON, 8-digit) | `2.6149705` |
  | Python reference (full precision) | `2.6149704969056877` |
  | Independent numpy re-implementation | `2.6149704969056886` |

  `|numpy − Python| = 8.9e-16` (machine epsilon); `|numpy − C++| = 3.1e-9`
  (the 8-digit JSON print rounding). All 16 bins matched field-by-field
  (`Lo`, `Hi`, `N`, `Pbar_M`, `Hbar`, `Var_H`, `Denom`, `b`, `Ne_Bin`, `Status`)
  to within the print precision, with zero status mismatches. The Python
  reference was written from the paper's Table-2 description, not ported from
  the C++.
- **Default path is unchanged**: with the flag off, the top-level `Ne` (`2.99`)
  and CI (`[2.92, 3.06]`) are identical to a run with the flag on, and the
  `Decode_Style_Kimura` block is absent. The estimator is a pure report-layer
  addition.
- **JSON robustness**: verified valid JSON for the flag alone (no
  `--cross-check`), for an all-bins-degenerate cohort (`Ne_Decode_Style` →
  `"NaN"`, exit `0`), and for all three optional blocks co-existing
  (`Kimura_Cross_Check` + `Decode_Style_Kimura` + `Per_Family_Estimates`).

---

## Tests

- **3 new GoogleTest cases** (`NeEstDecodeStyle.*`): one exercising every
  bin-status branch (`ok` / `single` / `empty` / `b<=0` / `b>=1`) against
  full-precision expectations from the Python reference, one empty-cohort edge
  case, and one large-cohort case asserting the **diffusion** inversion
  (`-1/ln b ≈ 2.47`, not the discrete `3.0`) is used.
- `ctest` **3/3** green on macOS / arm64 (Apple clang), across **189**
  GoogleTest cases (186 in v1.12.0 + 3 new).

---

## Files changed (main)

| Area | Files |
|---|---|
| ne-estimate | `src/ne_estimate.cpp`, `src/ne_estimate.h` (`DecodeStyleCheck` result struct + `compute_decode_style_check`), `tests/test_ne_estimate.cpp` |
| build / version | `CMakeLists.txt` (project VERSION 1.12.0 → 1.12.1) |
| docs | `release_notes/release_v1.12.1.md` |
