# mitoquest v1.12.4 — trans-prep REF-only sites; variant-qc compressed output

## Highlights

This **patch release** is a follow-up to v1.12.3 (which introduced
`caller --output-ref-sites`). It makes the downstream toolchain *aware* of the
REF-only records that switch produces, and adds output-compression control to
`variant-qc`. All changes are non-breaking: default behaviour is preserved, and
every new capability is either opt-in or purely additive.

Three things changed:

1. `mitoquest trans-prep` no longer silently drops REF-only records.
2. `mitoquest variant-qc` decides output compression from the file suffix.
3. `mitoquest caller` tags its command line in the VCF header with the version.

---

## What changed

### 1. trans-prep: faithful reporting of REF-only sites

When the input VCF was produced with `mitoquest caller --output-ref-sites`, it
contains REF-only (non-variant) records (`ALT=.`, `FILTER=REF_ONLY`). Previously
`trans-prep` skipped every record with `n_alt == 0`, so non-variant sites were
invisible in the transmission table. Now each REF-only record yields one row per
mother–child pair (and grandmother column when present) with:

| Column | Value |
|---|---|
| `ALT` | `.` |
| `MOTHER_AD_ALT` / `CHILD_AD_ALT` | `0` |
| `VAF` | `0` |
| `*_AD_REF` / `*_DP` | the real REF depth |

Details:

- With `--require-pass`, a REF-only record is kept when its `FILTER` is
  `REF_ONLY` (raw `caller` output) or `PASS` (after `variant-qc` rewrites
  FILTER); other filters (e.g. `LOW_QUALITY`) are dropped, exactly as for
  variant sites. Without `--require-pass`, REF-only records are kept regardless.
- `--snv-only` does not remove REF-only rows (they carry a single-base REF).
- The usual `GT='.'` and `--min-depth` rules apply unchanged.
- The summary on STDERR now reports the REF-only record count separately
  (`Processed N variant records and M REF-only records.`).
- Interaction with `ne-estimate`: these rows have maternal VAF 0, so the default
  `--min-vaf 0.10` window filters them out; they only enter the fit if you lower
  `--min-vaf` to 0. Filter on `ALT=.` if you do not want that.

### 2. variant-qc: suffix-controlled output compression

`-o/--output-vcf` and `-t/--output-tsv` now choose their compression from the
file suffix instead of always writing plain text:

| Suffix | Written as |
|---|---|
| `.gz` / `.bgz` | BGZF (gzip-compatible, block-gzipped) |
| anything else (incl. no suffix) | plain text |

- `-i/--input-vcf` accepts plain or compressed VCF/BCF.
- No index is created automatically; `.gz` VCF output can be indexed with
  `tabix` and read back by position.
- The TSV writer was migrated from `std::ofstream` to the project's BGZF writer
  so both sinks share the same compression path.

### 3. caller: version-stamped command header

The command line recorded in the VCF header now carries the tool version, so a
VCF can be traced back to the exact mitoquest build that produced it:

```
##mitoquest(1.12.4)_caller_command=mitoquest caller ...
```

## Compatibility

- `trans-prep`: previously REF-only records were dropped; they are now emitted.
  Consumers that assume every row is a variant should filter `ALT=.` if needed.
  This aligns `trans-prep` with the `--output-ref-sites` feature from v1.12.3;
  without that switch the input VCF has no REF-only records, so output is
  unchanged.
- `variant-qc`: passing an output path without a `.gz`/`.bgz` suffix keeps the
  old plain-text behaviour. Existing uncompressed invocations are unaffected.
- `caller`: the only header difference is the version token in the command
  comment line.

## Files changed

| Area | Files |
|---|---|
| trans-prep | `src/trans_prep.cpp` |
| variant-qc | `src/variant_qc.cpp` |
| caller | `src/mt_variant_caller.cpp` |
| tests | `tests/test_trans_prep.cpp`, `tests/test_variant_qc.cpp` |
| docs | `README.md`, `UPDATELOG.md`, `release_notes/release_v1.12.4.md` |
| build / version | `CMakeLists.txt` (project VERSION 1.12.3 → 1.12.4) |

## Verification

- Full `mitoquest_tests` suite run from `tests/` (where the E2E fixtures under
  `tests/data/` resolve): all modules green, including the new
  `TransPrepRun.ReportsRefOnlySitesFaithfully`,
  `TransPrepRun.RefOnlySitesKeptWithoutPassGateOnlyWhenNotFiltered`, and the
  parameterized `Formats/VariantQCOutput.*` suite over `''`/`.gz`/`.bgz`/`.txt`
  (compression mode + byte-for-byte content equality with plain output).
