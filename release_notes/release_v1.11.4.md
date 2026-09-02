# mitoquest v1.11.4 — io hardening, variant-qc parity & ne-estimate performance

## Highlights

This is a **patch release** spanning ten commits since v1.11.3. It hardens the
`src/io` htslib wrapper layer (fixes + extensions across BAM/VCF/CRAM paths and
a ThreadPool deadlock), brings the C++ `variant-qc` output to full VCF parity
with the Python reference implementation, accelerates the continuous trio
likelihood with a quadrature-rule cache, and expands the test suite with new
GoogleTest units plus mpmath/pysam cross-check harnesses. No CLI or output
schema changes; version bump `1.11.3 → 1.11.4`.

---

## Bug Fixes

### `io` layer — wrapper hardening

- **ThreadPool deadlock fix**: a worker previously executed queued tasks while
  still holding the queue mutex, so threads could not observe the stop flag and
  wait forever at shutdown. Tasks are now run outside the lock; the old
  try/catch around task execution was removed because `submit()` uses
  `packaged_task`, whose exceptions are stored and re-thrown by `future::get()`.
- **VCF record / header fixes**: float-missing detection rewritten as an exact
  bit-pattern comparison against htslib's sentinel (`is_float_missing`), plus
  signed/unsigned comparison and bounds fixes across `vcf_record.cpp`.
- **BAM record / header / bgzf extensions** (`iobgzf`, `bam_record`,
  `bam_header`, `fasta`, `hts_utils`, `io/utils`): bgzf read/write wrapper
  additions and correctness fixes feeding the `subsam` and caller paths.
- **htslib submodule**: pointer updated (submodule bump commit).

### `variant-qc` — Python reference parity restored

C++ output now matches `tools/mtDNA_variant_QC.py` record-for-record
(8/8 VCF records identical in the three-way comparison). Reverted semantic
divergences (Plan A):

- SB-missing SRF default `1.0 → 0.5` (neutral prior).
- Removed the SRF resize padding block.
- Pre-filter no longer triggers on ploidy == 0 and no longer range-checks
  AF/AD (restores v1.11.3 semantics: only missing essential fields trigger).
- Background error rate uses unconditional `1 - sum(AF)` (no `sum_af` guard).
- `_iterative_fit_and_call` restores full-allele computation (removes the
  ref-skip), `n_alleles = min(ad, srf, aq)`, and GT updates mirror Python.

Additionally, `variant_qc.cpp` no longer calls the raw htslib API at all:
every `bcf_*` / `BCF_*` use was replaced by the `ngslib::VCFRecord` /
`VCFHeader` wrappers (unpack flags, `INT32_MISSING` / `INT32_VECTOR_END` /
`FLOAT_MISSING` sentinels, `FieldNumber`-based format lookup), and header
includes were trimmed to declaration-only dependencies.

### CI — cross-version gtest fix for `NeEstProfile`

`EXPECT_NEAR(-inf, -inf, 1e-12)` computes `fabs(-inf - -inf) = NaN`, which
fails on gtest releases without the same-sign-infinity special case (the apt
`libgtest-dev` 1.13 on the `ubuntu-latest` job) even though both sides are
genuinely equal. The Ne = 1 discrete profile row is `-inf` by design
(a segregating child has zero probability under complete drift). The assertion
now uses `EXPECT_DOUBLE_EQ` (bit-for-bit, all gtest versions). `ubuntu-latest`
CI is green again.

---

## Refinements

### `ne-estimate`

- **Trio quadrature cache**: Gauss–Jacobi rules for the grandmother–mother
  posterior depend only on (posterior α, β, node count); `TrioQuadratureCache`
  now caches the per-rule eigensolver output within one global-likelihood
  evaluation (double-checked locking, construction outside the lock,
  evaluation-scoped lifetime). The trio smoke benchmark runs roughly **4×**
  faster with output identical field-by-field.
- Input/`PairData` refinement: `locus_id` carried through rows, `LoadStats`
  reporting added to `load_pairs`, and the VAF-window arguments are validated
  and threaded through the API.
- Related docs: `ne_estimate.h` header comment now states the plug-in
  calibration caveat (MMLE biased down by ~3% at 500× and >10% at 100× depth;
  the trio path marginalises `p_mother` and attenuates but does not remove the
  bias).

### Caller / core

- `mt_variant_caller`, `basetype`, `mt_utils`, `algorithm`: base-calling and EM
  genotyping refinements, INFO annotation cleanup, and include hygiene
  (`mt_variant_caller.h` trimmed to declarations, `io/bam.h` moved into the
  `.cpp`).
- `mt_copynum`: region-parsing and GC-counting fixes.
- `trans_prep` / `vcf_subset_samples`: raw htslib calls replaced by io wrappers
  and sentinel-constant updates.

---

## Tests

- New GoogleTest suites: basetype, BAM aligned pairs, FASTA safety, mt utils,
  variant caller, VCF subset samples.
- Python cross-check harnesses registered with CTest: mpmath trio-likelihood
  comparison (`test_trio_likelihood_mpmath.py`) and ne-estimate JSON
  serialization checks.
- Extended ne-estimate / trans-prep / copynum / algorithm / io tests —
  **183 GoogleTest cases** in total; full `ctest` suite green on
  `ubuntu-latest` and `macos-latest`.

---

## Documentation

- `README.md`: `copynum` calibration semantics (autosomal length-normalized
  fragment density as reference, mtDNA ×2), `trans-prep` always-emitted wide
  schema and the three supported `FORMAT/AD` layouts (`Number=.` / `R` / `A`
  with `DP - sum(AD_ALT)` REF inference), and a formal description of the
  default continuous working Beta transition model (with the Ne = 1
  absorbing-state limit).
- Release notes moved under `release_notes/`.

---

## Files changed (main)

| Area | Files |
|---|---|
| io layer | `src/io/bam*`, `src/io/iobgzf.*`, `src/io/vcf*`, `src/io/fasta.*`, `src/io/hts_utils.h`, `src/io/utils.*` |
| threading | `src/external/thread_pool.h` |
| variant-qc | `src/variant_qc.*`, `tests/test_variant_qc.cpp` |
| ne-estimate | `src/ne_estimate.*`, `src/log_factorial.*`, `tests/test_ne_estimate.cpp`, `tests/trio_likelihood_probe.cpp` |
| caller/core | `src/mt_variant_caller.*`, `src/basetype.*`, `src/mt_utils.*`, `src/algorithm.*`, `src/trans_prep.*`, `src/mt_copynum.*`, `src/vcf_subset_samples.*` |
| docs/version | `README.md`, `CMakeLists.txt` (1.11.3 → 1.11.4), `release_notes/*` |
