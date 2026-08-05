
## mitoquest v1.11.3 — `variant-qc` correctness fixes & GT phase encoding

### Highlights

This is a **patch release** that corrects several logic bugs in the C++ `mitoquest variant-qc` subcommand so its results faithfully match the Python reference implementation (`tools/mtDNA_variant_QC.py`), plus a genotype phase-encoding fix in the VCF I/O layer. No CLI, output-schema, or API changes.

---

### Bug Fixes

#### `variant-qc` — background error rate now matches the Python reference

`_fetch_background_error_rates` previously dropped perfect-homoplasmy samples (`AF = 1.0`, error rate `0`) via a `(0, 1)` boundary filter and used only the first allele frequency (`1 - AF[0]`). This biased the per-site `p_error` mean upward, which propagates through the binomial threshold in `_collect_mutation_vafs` and can change which VAFs feed the global Beta fit.

- Now uses `1 - sum(AF)` for **every** qualifying haploid sample, regardless of GT, **without** the `(0, 1)` boundary filter — matching Python `fetch_background_error_rates`.
- Boundary values are legitimate observations for the `p_error` mean; the Beta fit (`fit_beta_mle`) still filters the `(0, 1)` interior on its own, mirroring Python `estimate_background_noise`.

#### VCF I/O — correct unphased GT encoding on write-back

`VCFRecord::update_genotypes` received **already-decoded** allele indices (produced by `get_genotypes` via `bcf_gt_allele`, which strips phase). The old code applied `bcf_gt_is_phased()` to those raw indices, so any **odd** allele index (`1, 3, …`) at position ≥ 2 was misencoded as *phased*, emitting a `|` separator instead of `/` (e.g. `0/1 → 0|1`, `0/1/2 → 0|1/2`).

- Now always encodes with `bcf_gt_unphased()`, since phase cannot round-trip through this interface.
- **Impact:** allele values were never affected (decoding strips the phase bit), and haploid genotypes were unaffected; only the phase separator of multi-allelic / heteroplasmy genotypes was wrong. This fix removes those incorrect phase assertions.

#### `variant-qc` — earlier fixes rolled into this release

- Use total depth **DP** as the Bayesian denominator and write filtered genotypes back to the output VCF.
- Stop-at-vector-end parsing for variable-length FORMAT arrays (AD/AF/AQ/FS/SOR), preventing spurious padding values.
- Strand-ratio penalty direction: `penalty_srf = srf` (biased strands → small SRF → heavy H1 penalty), correcting the previously inverted `1 - srf`.

---

### Validation

- C++ vs Python `variant-qc` on the test VCF: **GT and FILTER are 100% identical** across all sites and samples (the only remaining difference is a cosmetic `:PP:GOOD_CALL` suffix that pysam omits for blacklisted, never-assigned records).
- GT phase-encoding reasoning verified with a standalone program against htslib's own macros (`bcf_gt_unphased` / `bcf_gt_phased` / `bcf_gt_is_phased` / `bcf_gt_allele`).
- Full unit-test suite rebuilt; no regressions in `variant-qc`-related tests.

---

### Files changed

| File | Change |
|---|---|
| `src/variant_qc.cpp` | `_fetch_background_error_rates` uses `1 - sum(AF)` without `(0,1)` filter |
| `src/io/vcf_record.cpp` | `update_genotypes` always encodes unphased GT |
| `CMakeLists.txt` | Version bump 1.11.2 → 1.11.3 |
