# MitoQuest v1.8.0

## Highlights

This release introduces two new subcommands that together build a complete
**mtDNA bottleneck-size (Ne) estimation pipeline** for cohorts with mother–child
transmission pairs, and tightens VCF FORMAT metadata to follow the VCFv4.2
spec.

- **`mitoquest trans-prep`** — extract per-site mother→child mtDNA allele
  transmission pairs from a multi-sample VCF + a PLINK-style FAM file.
- **`mitoquest ne-estimate`** — estimate the mtDNA germline bottleneck size
  Ne from those pairs via Beta-Binomial maximum-likelihood, with a Kimura
  cross-check.
- **VCF FORMAT cardinality fix** — six per-sample-variable FORMAT fields
  (`AD`, `AF`, `AQ`, `LAF`, `FS`, `SOR`) emitted by `mitoquest caller` are
  now declared with the spec-correct `Number=.`.

## What's new

### `mitoquest trans-prep`

Extracts informative (mother, child, site) triples for downstream Ne
estimation. Honours `mitoquest caller`'s **GT-aligned per-sample variable**
AD layout: at each sample, `AD[i]` is the depth of the allele at GT position
`i`, and the number of values varies per sample (1 for homoplasmic, 2 for
heteroplasmic, 3+ for tri-allelic). GT='.' calls are correctly treated as
missing and counted in a dedicated drop bucket.

Output is a per-pair, per-site TSV with mother/child VAFs, depths, allele
labels, and quality flags (`PASS` / `LOW_DEPTH`).

### `mitoquest ne-estimate`

Estimates Ne under the Beta-Binomial sampling model
(Wonnapinij et al. 2008/2010), with:

- Bias-corrected MLE on informative heteroplasmic pairs.
- VAF filter rationale: drop `m_vaf == 0` and homoplasmic-mother sites
  (no Fisher information for Ne).
- Kimura distribution cross-check on the same pair set.
- JSON and human-readable reports.

### VCF FORMAT header standardization

Per VCFv4.2, `Number=A` denotes one value per ALT (fixed across samples) and
`Number=R` denotes one per REF+ALT (also fixed). `mitoquest caller` instead
emits **per-sample variable-length** values (one per allele present in that
sample's GT). The correct cardinality for that layout is `Number=.`.

Six FORMAT declarations are corrected accordingly:

| Field | Before        | After         |
|-------|---------------|---------------|
| AD    | `Number=A`    | `Number=.`    |
| AF    | `Number=A`    | `Number=.`    |
| AQ    | `Number=A`    | `Number=.`    |
| LAF   | `Number=A`    | `Number=.`    |
| FS    | `Number=A`    | `Number=.`    |
| SOR   | `Number=A`    | `Number=.`    |

Fixed-cardinality fields (`GT`, `GQ`, `DP`, `CI`, `SB`, `VT`) and INFO-level
per-ALT fields (`VAF_MEAN`, `VAF_MEAN_HET`) are unchanged.

## Tests

- New `TransPrep` GoogleTest suite (8 cases) covering AD/GT alignment,
  multi-allelic sites, GT='.' missing-call drops, low-depth gating, FAM
  parsing, and end-to-end `run()` correctness.
- New `NeEstimate` GoogleTest suite covering MLE convergence, edge cases,
  and JSON output.
- `tests/data/ne_pipeline/synthesize.py` synthesizes cohort VCFs with
  configurable bottleneck size, multi-allelic rate, and a new
  `--missing-gt-rate` knob (default 5%) that injects truncated `GT='.'`
  calls to exercise the missing-genotype drop path end-to-end.
- End-to-end `run_demo.sh` reproduces the published Ne under noisy
  conditions (Ne=5 recovered with ~5% per-sample missing GT).

Full TransPrep test suite: **8/8 passing**. Full project test suite excluding
pre-existing data-dependent BAM/CopyNum tests: **passing**.

## Compatibility

- **Adds** two new subcommands; no existing CLI is changed.
- The FORMAT cardinality change (`Number=A` → `Number=.`) is a metadata-only
  correction. Conforming VCF readers (bcftools, htslib, pysam, GATK) handle
  `Number=.` for variable per-sample values without issue. Downstream tools
  that hard-coded `Number=A` parsing for these fields against `mitoquest
  caller` output should be re-checked, but in practice the value strings
  are identical — only the declared cardinality changed.
- `trans-prep` requires FORMAT/GT and FORMAT/AD to be present on each record.
