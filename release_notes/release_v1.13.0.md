# mitoquest v1.13.0 — variant-qc: posterior-based QUAL, INFO/KL, reworked H1 training

## Highlights

This **minor release** reworks `mitoquest variant-qc`. It contains two
behaviour changes that can alter existing results, so it is **not** a drop-in
patch:

1. `QUAL` in the output VCF is now a Phred-scaled posterior score instead of a
   rounded KL divergence; the KL value moves to a new `INFO/KL` field.
2. H1 (Beta) parameter training was restructured: new ALT-selection and
   convergence rules, a frozen final model applied once to all records, and an
   optional training trace.

GT, `FORMAT/PP`, `GOOD_CALL`, `FILTER` rules and the TSV report format are
unchanged (values may still shift because of the training changes in section 2).

---

## What changed

### 1. QUAL and INFO/KL

| Field | Before | Now |
|---|---|---|
| `QUAL` | rounded site KL divergence | `-10·log10(1 - p*)`, `p*` = strongest original-ALT posterior |
| `INFO/KL` | — | histogram KL divergence (natural log, not rounded), diagnostic only |

- `QUAL` is computed from the final-model log-posterior odds `z` with the stable
  form `(10/ln10)·[max(z,0) + ln(1+e^{-|z|})]`, so it is not capped by `PP`
  rounding to 1.
- It is computed **before** the PP threshold: a rejected ALT still has a score,
  and a high REF posterior cannot inflate it.
- Only original non-reference GT entries with a declared ALT and valid
  AD/SRF/AQ contribute. Records with no scoreable ALT (REF_ONLY, skipped
  records) have `QUAL=.`.
- `QUAL` is an approximate strongest-ALT evidence score, not a calibrated joint
  site-error probability.
- `INFO/KL` is `.` for skipped records; a processable REF-only record can have
  `KL=0`. `--bins` affects KL only.
- The VCF header gains `##INFO=<ID=KL,...>` and `##variant_qc_qual=...`;
  the `FILTER=LOW_QUALITY` description now reads "No sample posterior exceeds
  the QC threshold" (the logic itself is unchanged).
- Site metrics are stored per input record rather than keyed by `(chrom,pos)`,
  so multiple records at one position no longer overwrite each other.

### 2. H1 training rework

- Refinement uses only evaluable ALT observations (`0 < AD < DP`); each ALT must
  pass its own posterior threshold to enter the next fit.
- Convergence requires all of: fraction of changed ALT calls, max absolute ALT
  posterior change, and max absolute change of `log(alpha_h1)` / `log(beta_h1)`
  to be below `--convergence-eps`.
- Training also stops when fewer than three ALT observations qualify, when no
  ALT is evaluable, or when the fit limit is reached; the log states which.
- After training, the frozen H1 model is applied once to all processable
  records (including `REF_ONLY` sites from `caller --output-ref-sites`).
- `--max-iter` now counts H1 fits **including** the initial fit.
- New option `--training-trace FILE` writes a per-fit trace TSV (`.gz`/`.bgz`
  supported); the path must differ from every input/output file.
- New helper `tools/plot_variant_qc_training.py` plots the trace.

## Compatibility

- **`QUAL` thresholds are not transferable.** Any pipeline that filtered on the
  old KL-based `QUAL` must be re-tuned; use `INFO/KL` for the old quantity
  (unrounded).
- Results can differ slightly from v1.12.x because of the new ALT-selection and
  convergence rules, even without REF-only sites in the input.
- `--max-iter 20` now means 20 fits including the initial one.
- No change to command-line options other than the added `--training-trace`.

## Files changed

| Area | Files |
|---|---|
| variant-qc | `src/variant_qc.cpp`, `src/variant_qc.h` |
| tools | `tools/plot_variant_qc_training.py` |
| tests | `tests/test_variant_qc.cpp`, `tests/test_plot_variant_qc_training.py` |
| docs | `README.md`, `UPDATELOG.md`, `release_notes/release_v1.13.0.md` |
| build / version | `CMakeLists.txt` (project VERSION 1.12.4 → 1.13.0) |

## Verification

- Rebuilt and ran `ctest` (AllTests, TrioLikelihoodMpmathCrossCheck,
  NeEstimateJsonSerialization): 3/3 passed.
- `python -m unittest tests.test_plot_variant_qc_training`: 21 tests OK.
