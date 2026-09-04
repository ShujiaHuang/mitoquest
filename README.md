<p align="center">
  <a href="https://github.com/ShujiaHuang/mitoquest">
    <img height="230" src="docs/assets/images/mitoquest_logo.svg">
  </a>
</p>

# MitoQuest: Human Mitochondrial sequencing data Analysis Toolkit

***MitoQuest*** is a cross-platform, efficient and practical bioinformatics
toolkit written in C++ (C++17) for analyzing **human mitochondrial DNA
(mtDNA)** from whole-genome / mtDNA-targeted sequencing (WGS) data. It calls
mitochondrial **SNVs and Indels**, quantifies **heteroplasmy and homoplasmy**,
estimates **mtDNA copy number**, and provides downstream utilities for VCF
sub-sampling, annotation, and quality recalibration.

MitoQuest is built on top of [htslib](https://github.com/samtools/htslib)
(vendored as a submodule) and is designed to scale from a single-sample
clinical workflow to population-level cohorts with thousands of samples.

```bash
mitoquest: Human Mitochondrial sequencing data Analysis Toolkit
Version: 1.11.0

Usage: mitoquest <command> [options]
Commands:
  caller       Mitochondrial variants and heteroplasmy/homoplasmy caller.
  subsam       Extract mitochondrial variants for specified samples from VCF files and output a new VCF file.
  copynum      Estimate per-chromosome (incl. mtDNA) relative copy number from a BAM/CRAM file.
  trans-prep   Extract mother-child mtDNA allele transmission pairs from a multi-sample VCF + FAM file.
  ne-estimate  Estimate the mtDNA bottleneck size (Ne) from transmission pairs via Beta-Binomial MMLE (Maximum Marginal Likelihood).

```

In addition to the main `mitoquest` binary, the project ships:

- `tools/` — a suite of Python helper scripts for VCF QC, annotation,
  pipeline assembly, and downstream analysis.
- `data/` — curated reference resources (population databases, in-silico
  predictors, RNA/protein domain annotations, blacklist, Phylotree
  variants, …) used by the annotation and QC tools.

---

## Installation and quick start

Pre-built static binaries are published on the
[GitHub Releases page](https://github.com/ShujiaHuang/mitoquest/releases) — **most user shoild simply download the binary and run**.

| Platform              | Download                                                                                                              | Notes                                  |
| --------------------- | --------------------------------------------------------------------------------------------------------------------- | -------------------------------------- |
| Linux (x86_64)        | [mitoquest-linux-static](https://github.com/ShujiaHuang/mitoquest/releases/latest/download/mitoquest-linux-static)    | Requires **glibc ≥ 2.35** (see below)  |
| macOS (arm64 / Intel) | [mitoquest-macos-static](https://github.com/ShujiaHuang/mitoquest/releases/latest/download/mitoquest-macos-static)    | Requires **macOS 12+**                 |

#### System requirements for `mitoquest-linux-static`

The Linux binary is a **partial-static** build, produced on Ubuntu 22.04
(glibc 2.35) in CI. It bundles `libstdc++`, `libgcc`, `htslib`, `zlib`,
`bzip2`, `xz`, and `openssl` statically — only the system C library
(**glibc**) is linked dynamically. Because glibc symbol versions are
**forward-compatible only**, the binary requires the host `glibc` version
to be **≥ 2.35**.

**Quick check on your machine:**

```bash
# If the printed glibc version is >= 2.35, mitoquest-linux-static will run.
ldd --version | head -1
```

If you see this — or you are on CentOS / RHEL / Rocky / AlmaLinux / older Ubuntu / older Debian — you can compile from source instead (see below).

```bash
# Linux
wget https://github.com/ShujiaHuang/mitoquest/releases/latest/download/mitoquest-linux-static
chmod +x mitoquest-linux-static
mv mitoquest-linux-static motiquest
./mitoquest --help
```

```bash
# macOS
curl -LO https://github.com/ShujiaHuang/mitoquest/releases/latest/download/motiquest-macos-static
chmod +x mitoquest-linux-static motiquest
mv mitoquest-linux-static motiquest
./mitoquest --help
```

> [!IMPORTANT]
> **Rename the downloaded binary** to `mitoquest` for convenience. You may also move it to a directory in your `$PATH` (e.g. `/usr/local/bin`) for system-wide access.

> [!TIP]
> **If the pre-built binary does not work on your system** (e.g. glibc < 2.35 on older Linux), you can compile from source instead. Detailed build instructions (CMake build, static build, and manual g++ fallback) are available in **[docs/INSTALL_FROM_SOURCE.md](docs/INSTALL_FROM_SOURCE.md)**.

## Commands overview

```bash
Usage: mitoquest <command> [options]

Commands:
  caller       Call mitochondrial variants (SNVs/Indels) and quantify
               heteroplasmy / homoplasmy from BAM/CRAM files.
  subsam       Extract a subset of samples from an mtDNA VCF and recompute
               INFO fields.
  copynum      Estimate per-chromosome (incl. mtDNA) relative copy number
               from a BAM/CRAM file.
  trans-prep   Extract mother-child mtDNA allele transmission pairs from a
               multi-sample VCF + a PLINK FAM file.
  ne-estimate  Estimate the mtDNA bottleneck size (Ne) from transmission
               pairs via the Beta-Binomial Maximum Marginal Likelihood Estimator (MMLE).
```

---

## `mitoquest caller` — Variant calling and heteroplasmy/homoplasmy detection

`mitoquest caller` reads aligned reads from one or more BAM/CRAM files and
emits a multi-sample VCF containing per-sample heteroplasmy fractions
(HF) and per-site INFO fields suitable for downstream filtering and
annotation.

### Full parameter reference of `caller`

```bash
Usage: mitoquest caller [options] -f ref.fa -o output.vcf.gz in1.bam [in2.bam ...]

Required options:
  -f, --reference FILE       Reference FASTA file
  -o, --output    FILE       Output VCF file (use .vcf.gz for bgzipped output)

Optional options:
  -b, --bam-list FILE        List of input BAM/CRAM filenames, one per line.
  -Q, --min-BQ INT           Skip bases with base quality smaller than INT
                             (default: 20).
  -q, --min-MQ INT           Skip alignments with mapQ smaller than INT
                             (default: 20).
  -r, --regions REG[,...]    Comma-separated regions to process
                             (default: entire mtDNA contig).
                             Format: chr | chr:start | chr:start-end
                             Example: chrM or chrM:1-1000,chrM:8000-8200
  -p, --pairs-map-only       Only use paired reads where mate maps to the
                             same chromosome.
  -P, --proper-pairs-only    Only use properly paired reads (SAM flag 0x2).
  --filename-has-samplename  When BAM/CRAM filenames are 'SampleID.xxxx.bam',
                             skip reading sample names from the BAM header.
                             Saves significant time on large cohorts.
  -j, --het-threshold FLOAT  Heteroplasmy fraction threshold below which an
                             allele is treated as reference (default: 0.01).
  -c, --chunk INT            Chunk size (bp) for parallel region processing
                             (default: 1000).
  -t, --threads INT          Number of threads (default: all available CPUs).
  -h, --help                 Show this help message and exit.
```

### Usage examples

**Single-sample variant calling on the entire mtDNA contig:**

```bash
mitoquest caller \
    -f reference.fasta \
    -o sample.mt.vcf.gz \
    sample.bam
```

**Multi-sample call from a list of BAM/CRAM files (with quality filters):**

```bash
mitoquest caller \
    -f reference.fasta \
    -o cohort.mt.vcf.gz \
    -b bamfile.list \
    -Q 30 -q 30 -t 16
```

**Recommended call with sample-name optimisation (large cohorts):**

```bash
# When BAMs are named `SampleID.bam` / `SampleID.cram`, --filename-has-samplename
# avoids opening every BAM just to read @RG SM tags.
mitoquest caller \
    -f reference.fasta \
    -o cohort.mt.vcf.gz \
    -Q 30 -q 30 -t 24 \
    --filename-has-samplename \
    -b bamfile.list
```

**Restrict calling to a specific region (e.g., the D-loop / control region):**

```bash
mitoquest caller \
    -f reference.fasta \
    -o cohort.dloop.vcf.gz \
    -r chrM:1-576,chrM:16024-16569 \
    -Q 30 -q 30 -t 16 \
    --filename-has-samplename \
    -b bamfile.list
```

**Stricter heteroplasmy detection (lower the HF threshold):**

```bash
# Default heteroplasmy threshold is 0.01 (1%).  Lower it to 0.005 to capture
# very low-frequency heteroplasmy (use with caution; requires high coverage).
mitoquest caller \
    -f reference.fasta \
    -o cohort.lowhet.vcf.gz \
    -j 0.005 \
    -Q 30 -q 30 -t 24 \
    --filename-has-samplename \
    -b bamfile.list
```

**Use only properly paired reads (more stringent for NUMT filtering):**

```bash
mitoquest caller \
    -f reference.fasta \
    -o cohort.proper.vcf.gz \
    -P --filename-has-samplename \
    -Q 30 -q 30 -t 16 \
    -b bamfile.list
```

**Mix a file list with extra BAMs on the command line:**

```bash
mitoquest caller \
    -f reference.fasta \
    -o cohort.vcf.gz \
    -b bamfile.list \
    -Q 30 -q 30 -t 16 \
    extra1.cram extra2.bam
```

---

## `mitoquest subsam` — Extract samples from an mtDNA VCF

Extract a subset of samples from a multi-sample VCF/BCF file, optionally
recomputing INFO fields (AC/AN/AF/NS, …) so they reflect only the kept
samples.

### Full parameter reference of `subsam`

```bash
Usage: mitoquest subsam [options] -i <input.vcf> -o <output.vcf> [-s <samplelist>] [<sample1> ...]

Options:
  -i, --input FILE     Input VCF/BCF file (required).
  -o, --output FILE    Output VCF/BCF file (required).
  -s, --sample FILE    File listing sample names to keep, one per line.
  -O, --output-type    Output file type: v|z|b|u
                       (v: VCF, z: bgzipped VCF, b: BCF, u: uncompressed BCF).
                       Default: inferred from output filename extension.
  --no-update-info     Do not recalculate AC/AN/AF/NS in the INFO column.
  --keep-all-site      Retain sites that become reference-only after
                       subsetting (default: drop them).
  -h, --help           Show this help message and exit.
```

### Usage examples for extracting samples

**Extract samples listed in a file (output bgzipped VCF):**

```bash
mitoquest subsam \
    -i cohort.mt.vcf.gz \
    -o subset.mt.vcf.gz \
    -s sample_names.txt
```

**Extract two specific samples directly from the command line:**

```bash
mitoquest subsam \
    -i cohort.mt.vcf.gz \
    -o pair.mt.vcf \
    -O v \
    SampleA SampleB
```

**Keep all sites (including those that become ref-only) and skip INFO updates:**

```bash
# Useful when you want to preserve a fixed multi-cohort site list, e.g. for
# downstream join-call merging.
mitoquest subsam \
    -i cohort.mt.vcf.gz \
    -o subset.mt.vcf.gz \
    -s sample_names.txt \
    --keep-all-site --no-update-info
```

---

## `mitoquest copynum` — mtDNA copy-number estimation

`mitoquest copynum` estimates **per-chromosome relative copy number** from a
sorted/indexed BAM or CRAM file. The mean autosomal length-normalized
fragment density is the reference value: autosomal `CopyNum` values are
reported relative to that reference (near 1), whereas mitochondrial values
are multiplied by 2 so that equal autosomal and mtDNA molecular densities
map to 2 mtDNA copies per diploid autosomal baseline. This is a relative
library-depth calibration, not an absolute molecule count or a GC-corrected
quantification. The output TSV reports fragment counts, GC content as a
diagnostic, length-normalized fragment ratio, and copy-number mean + 95%
confidence interval for every contig in the BAM header.

### Full parameter reference

```bash
Usage: mitoquest copynum [options] <input.bam/cram>

Options:
  -r, --reference FILE   Reference genome FASTA file (required; needed for
                         CRAM decoding and GC content calculation).
  -o, --output    FILE   Output TSV file (default: stdout).
  -q, --mapq      INT    Minimum mapping quality score [0].
  -t, --threads   INT    Number of worker threads [hardware_concurrency].
  -s, --seqtype   STR    Sequencing type: auto | pe | se [auto].
                         pe = paired-end, se = single-end.
  -L, --regions   STR    Restrict counting / GC / length normalization to
                         the listed regions, given either as a comma-
                         separated list (e.g. 'chrM:1-300,chrM:16000-16569')
                         or as a path to a file (one region per line, in
                         either 'chr:start-end' samtools form or BED-style
                         whitespace-delimited `chr start end` triples; '#'
                         starts a
                         comment). The intervals REPLACE the whole-
                         chromosome window for the chromosomes they cover;
                         chromosomes not mentioned keep their full-length
                         behaviour. Typical use: exclude NUMT-affected zones
                         from the mtDNA copy-number estimate.
  -h, --help             Show this help message and exit.
```

### Usage examples for estimating copynum

**Estimate copy numbers from a BAM file (paired-end auto-detected):**

```bash
mitoquest copynum -r reference.fasta sample.bam > sample.cn.tsv
```

**Estimate from a CRAM file with a stricter MAPQ filter and 8 threads:**

```bash
mitoquest copynum \
    -r reference.fasta \
    -q 30 -t 8 \
    -o sample.cn.tsv \
    sample.cram
```

**Force single-end counting (useful for legacy unpaired data):**

```bash
mitoquest copynum \
    -r reference.fasta \
    -s se -q 30 -t 4 \
    sample.bam > sample.cn.tsv
```

**Exclude NUMT-affected zones via inline regions (recommended for mtCN):**

```bash
# Restrict chrM measurement to the two stretches least affected by NUMTs;
# autosomes are still measured over their full length, so the diploid
# baseline is unchanged.
mitoquest copynum \
    -r reference.fasta \
    -L 'chrM:1-300,chrM:16000-16569' \
    -q 30 -t 8 \
    sample.bam > sample.cn.tsv
```

**Exclude NUMT zones via a BED file (one region per line):**

```bash
# my_chrM_regions.bed (0-based half-open, standard BED):
#   chrM    0       300
#   chrM    15999   16569
mitoquest copynum \
    -r reference.fasta \
    -L my_chrM_regions.bed \
    -q 30 -t 8 \
    sample.bam > sample.cn.tsv
```

The TSV output has one row per contig in the BAM header, with the
following columns:

```
#Chromosome  Fragments  Chrom_Length  GC_Content  Fragment_Normalized_Ratio  CopyNum  CopyNum-CI95-Lower  CopyNum-CI95-Upper  Effective_Length  Regions_Used
```

- `Chrom_Length` always reports the full contig length from the BAM header.
- `Effective_Length` is the number of bases that were actually counted
  (equals `Chrom_Length` when `-L` did not target that contig; otherwise it
  is the sum of the merged user-supplied intervals on that chromosome).
- `Regions_Used` is `.` for unrestricted chromosomes, or a comma-separated
  list of merged `start-end` intervals (1-based inclusive) when `-L` was
  applied to that chromosome.
- When `-L` is supplied, the output also includes a `#Regions argument: ...`
  header comment recording the original CLI value for reproducibility.

---

<!-- 

## `mitoquest trans-prep` — Extract mother-child transmission pairs

`mitoquest trans-prep` walks a multi-sample VCF (typically the output of
`mitoquest caller`) together with a PLINK-format FAM file describing the
trios, and emits a per-allele TSV of mother-child mtDNA transmission
pairs. The TSV is the input format expected by `mitoquest ne-estimate`.

Only mother-child relationships are retained; fathers do not transmit
mtDNA in mammals and are deliberately ignored.

### Full parameter reference

```bash
Usage: mitoquest trans-prep [options] -v <vcf> -f <fam> -o <pairs.tsv>

Required options:
  -v, --vcf         FILE   Input multi-sample VCF/BCF (must declare
                           FORMAT/GT, FORMAT/DP and FORMAT/AD; see
                           "AD interpretation" below).
  -f, --fam         FILE   PLINK 6-column FAM file: family_id, child_id,
                           father_id, mother_id, sex, phenotype
                           (whitespace-delimited).
  -o, --output      FILE   Output TSV file (use `-` for stdout).

Optional options:
  -d, --min-depth   INT    Minimum DP at both mother and child for the
                           pair to be tagged QC=PASS [100].
      --require-pass       Only keep VCF records whose FILTER == PASS
                           (default: enabled).
      --no-require-pass    Disable the FILTER=PASS gate.
      --snv-only           Restrict to SNV records (default: enabled).
      --no-snv-only        Disable the SNV-only gate (also keeps indels).
  -g, --gm-fam      FILE   Optional GM FAM file describing grandmother–
                           mother relationships (same 6-column PLINK
                           format).  When set, each MC trio is looked up
                           as a child in the GM FAM; if the MC mother is
                           found as a child with a defined mother there,
                           the row is tagged HAS_G=1 and grandmother
                           genotype columns are emitted.
  -h, --help               Print this help message.
```

### FAM format

The FAM file is the standard PLINK 6-column whitespace-delimited format:

```
FAM_ID  CHILD_ID  FATHER_ID  MOTHER_ID  SEX  PHENOTYPE
```

- `MOTHER_ID == 0` denotes "unknown"; that line is skipped (counted under
  `ignored_no_mother_in_fam` in the matching report).
- Sample IDs are matched **case-sensitively** against the VCF sample names.
- A FAM line whose mother or child is not present in the VCF is dropped
  (counted under `ignored_missing_mother` / `ignored_missing_child`).
- The FATHER_ID column is preserved in the report but otherwise ignored.

### Output TSV columns

One row is emitted per (transmission-pair, ALT allele) combination.
Multi-allelic sites produce one row per ALT. The output always uses the
wide schema with grandmother genotype fields and a `HAS_G` flag: without a
matched `--gm-fam` relationship, the grandmother fields are `NA`/zero and
`HAS_G=0`.

| Column         | Type   | Description                                          |
| -------------- | ------ | ---------------------------------------------------- |
| `CHROM`        | string | Reference contig (e.g. `chrM`).                      |
| `POS`          | int    | 1-based position.                                    |
| `REF`          | string | Reference allele.                                    |
| `ALT`          | string | One ALT allele.                                      |
| `FAMILY_ID`    | string | FAM family ID.                                       |
| `GRANDMOTHER_ID` | string | Grandmother sample ID (`NA` when `HAS_G=0`).  |
| `MOTHER_ID`    | string | Mother sample ID (must match VCF sample column).     |
| `CHILD_ID`     | string | Child sample ID (must match VCF sample column).      |
| `GRANDMOTHER_DP` | int  | Grandmother total depth (0 when `HAS_G=0`).         |
| `GRANDMOTHER_AD_REF` | int | Grandmother REF reads (0 when `HAS_G=0`).      |
| `GRANDMOTHER_AD_ALT` | int | Grandmother ALT reads (0 when `HAS_G=0`).      |
| `GRANDMOTHER_VAF` | float | `GRANDMOTHER_AD_ALT / GRANDMOTHER_DP` (0 when `HAS_G=0`). |
| `HAS_G`        | int    | 1 = G-M-C trio (3-gen); 0 = MC pair (2-gen).        |
| `MOTHER_DP`    | int    | Mother total depth at this site.                     |
| `MOTHER_AD_REF`| int    | Mother reads supporting REF.                         |
| `MOTHER_AD_ALT`| int    | Mother reads supporting ALT.                         |
| `MOTHER_VAF`   | float  | `MOTHER_AD_ALT / MOTHER_DP`.                         |
| `CHILD_DP`     | int    | Child total depth at this site.                      |
| `CHILD_AD_REF` | int    | Child reads supporting REF.                          |
| `CHILD_AD_ALT` | int    | Child reads supporting ALT.                          |
| `CHILD_VAF`    | float  | `CHILD_AD_ALT / CHILD_DP`.                           |
| `QC`           | string | `PASS`, `LOW_DEPTH`, or other failure reason.        |

Legacy 16-column TSVs without the grandmother fields are still accepted by
`ne-estimate`; a TSV that contains any trio column must contain all of
`HAS_G`, `GRANDMOTHER_DP`, and `GRANDMOTHER_AD_ALT`. Declared trio rows with
invalid grandmother counts are excluded rather than silently treated as
ordinary mother-child rows.
A leading provenance comment line (`#mitoquest_trans_prep_command=...`)
records the exact CLI invocation for reproducibility.

### Usage example

```bash
# Step 1: variant calling on a cohort that includes mother-child trios
mitoquest caller \
    -f rCRS.fasta \
    -o cohort.vcf.gz \
    -Q 30 -q 30 -t 24 \
    -b bamfile.list

# Step 2: extract per-allele transmission pairs
mitoquest trans-prep \
    -v cohort.vcf.gz \
    -f cohort.fam \
    -d 500 \
    -o cohort.transmission_pairs.tsv
```

`mitoquest trans-prep` writes a small "matching report" to STDERR
summarising how many FAM lines were retained / dropped and why; review
it before downstream analysis.

### Caveats

- **Garbage in, garbage out**: the QC of the upstream variant calls
  fundamentally determines the quality of the Ne estimate. Always run
  this on a high-quality, filtered VCF (e.g. after `mitoquest variant-qc`).
- **AD interpretation is GT-aligned**: see "FORMAT/AD layout" below.
  Any `Number=` declaration is accepted; AD is always decoded
  per-sample using `FORMAT/GT`.
- **Sample IDs are case-sensitive** and must match exactly between the
  FAM file and the VCF sample column.

### Multi-allelic sites

`mitoquest trans-prep` decomposes each multi-allelic VCF record into
one output row **per ALT × per trio**.  AD is decoded per-sample
using `FORMAT/GT` (see "FORMAT/AD layout" below), so a tri-allelic
SNV like `chrM 1000 . A G,T` produces two rows per trio, one for the
`A>G` allele and one for the `A>T` allele.

* **Pure multi-allelic SNVs** (e.g. `A>G,T,C`) — every ALT is kept.
* **Mixed multi-allelic** (e.g. `A>G,GT`, one SNV + one indel) — with
  the default `--snv-only` flag, the SNV ALT is kept and the indel ALT
  is silently skipped within the same record. Disabling `--snv-only`
  keeps both.

Downstream, `mitoquest ne-estimate` treats each row as an independent
`(REF reads, ALT reads)` 2-allele observation — the standard biallelic
decomposition. This is exact for any single ALT and is a mild
approximation only when two heteroplasmic ALTs co-segregate in the
same mother (rare in mtDNA).

The `tests/data/ne_pipeline/` directory ships a synthetic dataset that
includes 2 tri-allelic SNV sites, so you can verify multi-allelic
handling end-to-end with one command:

```bash
bash tests/data/ne_pipeline/run_demo.sh
```

### FORMAT/AD layout (GT-aligned, per sample)

`mitoquest caller` emits `FORMAT/AD` with a **GT-aligned per-sample
layout** rather than the standard `Number=R` layout: for each sample,
`AD[i]` is the read depth of the allele at GT position *i*, and only
alleles present in that sample's `GT` are listed.  Examples:

| Sample call           | Meaning                                          |
| --------------------- | ------------------------------------------------ |
| `0/1:DP:r,a`          | heteroplasmic; `r`=REF reads, `a`=ALT reads      |
| `0:DP:r`              | homoplasmic REF; only REF depth recorded         |
| `1:DP:a`              | homoplasmic ALT; only ALT depth recorded         |
| `0/1/2:DP:r,a1,a2`    | tri-allelic heteroplasmy; depths in GT order     |
| `.:GQ:DP` (truncated) | GT missing / no variant evidence; AD absent      |

`mitoquest trans-prep` decodes AD by walking each sample's `GT`
and mapping AD positions back to allele indices.  An allele *not
present* in a sample's GT yields **0 supporting reads** for that
sample (this is the canonical "transmission-loss" case where a
heteroplasmic mother transmits homoplasmic REF or ALT to the child).

`mitoquest trans-prep` supports three explicit `FORMAT/AD` layouts:

| Header declaration | Interpretation |
| --- | --- |
| `Number=.` | MitoQuest GT-aligned variable-length values |
| `Number=R` | Standard fixed-width REF plus all ALT depths by allele index |
| `Number=A` | Standard fixed-width ALT depths by allele index; REF is inferred as `DP - sum(AD_ALT)` |

The `Number=A` inference requires non-negative ALT depths whose sum does not
exceed `DP`; otherwise the affected transmission row is dropped as malformed.

### Missing genotypes (`GT='.'`)

When a sample's `FORMAT/GT` is `'.'`, AD positions cannot be mapped
to allele indices.  In that case `mitoquest trans-prep` drops the
`(mother, child)` pair at that site.  The count of dropped
`(pair, site)` observations is reported on STDERR after the run as:

```
[trans-prep] (pair, site) dropped due to GT='.' or AD/GT mismatch: <N>
```

The same drop also fires when the lengths of `GT` and `AD` for a
sample disagree (which prevents an unambiguous GT-aligned AD lookup).

---

## `mitoquest ne-estimate` — Bottleneck size (Ne) Maximum Marginal Likelihood Estimation (MMLE)

`mitoquest ne-estimate` fits the mitochondrial **transmission bottleneck
size** Ne from the per-allele TSV produced by `mitoquest trans-prep`,
using a count-based **Maximum Marginal Likelihood Estimator (MMLE)**. The
default continuous model plugs the maternal count frequency into a working
Beta transition model and analytically marginalizes the child's latent
frequency; the alternative discrete model also integrates a maternal Beta
posterior. Rows are treated as independent working observations, so the
global objective is a composite log likelihood. We label the resulting
estimator **MMLE** to keep those assumptions explicit.

### Statistical model

#### Default continuous working Beta transition model

For each mother-child (M-C) row, the default `--model continuous` uses the
count-derived maternal plug-in frequency

$$
\hat p_M=\frac{k_M}{d_M},
$$

then models the child as

```txt
  p_child | p_hat_mother ~ Beta(p_hat_mother × (Ne − 1),
                                 (1 − p_hat_mother) × (Ne − 1))
  c_alt   | p_child      ~ Binomial(c_dp, p_child)
```

Marginalizing $p_C$ gives the per-row Beta-Binomial likelihood

$$
P(k_C\mid\hat p_M,N_e)=
\binom{d_C}{k_C}
\frac{B\{\hat p_M(N_e-1)+k_C,\ (1-\hat p_M)(N_e-1)+d_C-k_C\}}
{B\{\hat p_M(N_e-1),\ (1-\hat p_M)(N_e-1)\}}.
$$

This is a mean/variance-matched working transition approximation, not the
stationary distribution or exact transient density of a Wright-Fisher
diffusion. At $N_e=1$ it has its exact absorbing-state limit:

$$
P(k_C=0)=1-\hat p_M,\qquad P(k_C=d_C)=\hat p_M,
$$

and any interior child count has probability zero. If no candidate in the
configured range has finite likelihood, `ne-estimate` exits with an error
rather than reporting a fictitious optimum.

The cohort objective is the sum of these row log likelihoods. The optimizer
scans all integer candidates, preserves an exact constrained optimum at a
search boundary, and otherwise refines the neighboring interval with a
golden-section search. Its nominal 95% profile interval is the connected
region satisfying $\log L(N_e)\ge\log L(\widehat N_e)-1.92$; because this
is a composite likelihood, it is model-based rather than automatically
calibrated for linkage or repeated-family dependence.

#### Integrating the maternal frequency out (the default; `--no-maternal-marginalization` opts out)

The plug-in $\hat p_M=k_M/d_M$ treats the maternal *estimate* as if it were
the maternal *true* frequency, so the mother's own read-sampling noise is
mistaken for drift. The resulting fitted value is a pseudo-true $N_e^{\text{app}}$
obeying

$$
\frac{1}{N_e^{\text{app}}}\;=\;\frac{1}{N_e}+\frac{1}{d_M}
\qquad\Longleftrightarrow\qquad
\frac{N_e^{\text{Kimura}}}{\widehat N_e^{\text{MMLE}}}\approx 1+\frac{N_e}{d_M}
$$

at a **single** maternal depth $d_M$. Measured on synthetic data
($N_e^{\text{true}}=30$, $d_C=2000$, 240 pairs per cohort, 20 replicate
cohorts per depth): $-45.0\%$ at $d_M=30$, $-21.2\%$ at $d_M=100$,
$-7.4\%$ at $d_M=500$, $-0.6\%$ at $d_M=2000$. That depth dependence is why
plug-in estimates from cohorts sequenced at different coverage are not
directly comparable, and why the plug-in is no longer the default.

The default removes that term by modelling it instead of deleting rows. It
replaces the plug-in child factor with an integral over the maternal
read-sampling posterior (uniform Beta(1,1) prior, the same prior
`--model discrete` uses), evaluated by adaptive Gauss-Jacobi quadrature:

$$
\log L_i(N_e)=\log\binom{d_M}{k_M}+\log B(k_M{+}1,\,d_M{-}k_M{+}1)
+\log E_{p_M\sim\text{Beta}(k_M+1,\,d_M-k_M+1)}
\big[P(k_C\mid p_M,N_e)\big].
$$

It is **on by default** — `--no-maternal-marginalization` restores the legacy
plug-in — applies to the continuous model only, and reaches the per-family and
`--ne-profile` paths as well as the cohort fit. Four things to know:

1. **It is slower, but no longer by orders of magnitude in wall-clock terms.**
   Measured on the shipped 236-pair demo cohort (`tests/data/ne_pipeline`,
   uniform 2000× depth, default `--min-ne 1 --max-ne 100`, single thread,
   macOS arm64, best of 3): 0.005 s plug-in vs. 0.382 s marginalised for the
   cohort fit (80×), 0.007 s vs. 0.825 s with `--per-family` (126×), and
   0.016 s vs. 2.723 s with `--ne-profile` at the default 0.1 step, i.e. 992
   grid points (170×). Hoisting the maternal Gauss-Jacobi rule to the lifetime
   of one $N_e$ scan — the rule depends on the row and the quadrature order but
   not on $N_e$ — cut those three from 11.99 s / 27.45 s / 101.5 s, a 31–37×
   saving, and that saving is what makes marginalisation affordable as the
   default. The plug-in wall clock sits at the ~5 ms process-startup floor, so
   the ratio overstates the compute gap; the absolute cost is what matters and
   it is sub-second.
2. **`Max_Marginal_LogLik` is not comparable across the flag.** Turning it on
   replaces the whole child term and adds the maternal read and posterior
   normalisation, so the likelihood moves by hundreds of log units and the
   fitted $N_e$ can move too. Only the normalisation
   $\log\binom{d_M}{k_M}+\log B(k_M{+}1,d_M{-}k_M{+}1)=-\log(d_M+1)$ is
   $N_e$-independent. The JSON `Maternal_Marginalization` field records which
   likelihood the reported value belongs to.
3. **It does not repair a mis-specified bottleneck.** If the true bottleneck
   draws an integer number of genomes, the continuous transition is the wrong
   model *form* and both the plug-in and the marginalised fit stay biased; use
   `--model discrete` in that regime. It also leaves a residual
   prior-mismatch bias (the Beta(1,1) prior does not match a cohort whose true
   maternal frequencies are truncated away from 0 and 1 by the default VAF
   window), so it is **not exactly unbiased**. Over the same sweep as above it
   measures $+14.1\%$ at $d_M=30$, $\le 2.4\%$ in magnitude at every
   $d_M\ge 60$, and $-2.1\ldots+2.4\%$ for $d_M\ge 500$ against a Monte Carlo
   se of $1.3\ldots2.5\%$ — at deep coverage the residual is not separable from
   zero and its **sign flips with $d_M$** ($-2.1$ at 500, $+2.4$ at 1000,
   $+0.9$ at 2000, $+1.6$ at 5000) rather than staying positive. Redrawing the
   maternal frequencies near the VAF-window edges instead of the interior gives
   $+13.0\%$ / $\le 5.0\%$ / $-0.3\ldots+1.8\%$ for the same three brackets.
   At $d_M=100$ it turns the plug-in's
   $-21.2\%$ into $+1.7\%$. A `NOTE:` is printed on stderr when the harmonic
   maternal depth of the two-generation rows falls below 60, and `--min-depth`
   is the lever for dropping those rows.
4. **It is a no-op on rows that already integrate $p_M$.** Trio rows
   (`HAS_G=1`) integrate $p_M$ against the grandmother-informed posterior, so
   the flag does not touch them; and `--model discrete` already integrates
   $p_M$, so combining the two prints a `WARNING:` and ignores the flag.

Under the flag a maternal `0 / d_M` row is **no longer** uninformative: it
still leaves a proper $\text{Beta}(1,d_M{+}1)$ posterior against which the
child's count carries information. Those rows are excluded by the default
`--min-vaf 0.10` window anyway and only enter under `--min-vaf 0`.

`--min-depth INT` is the complementary **QC** lever rather than a model fix:
it drops rows whose mother *or* child depth falls below the threshold (both
sides are gated, 0 disables). The default already *models* the maternal read
noise; `--min-depth` *removes* the rows whose residual prior-mismatch bias is
largest, so the two are complementary rather than alternatives. Both change
`Pairs_Used` and the retained depth composition, so read them together
with the `Min_Depth` / `Mean_*_DP` / `Harmonic_*_DP` fields in the JSON.
Note that those depth fields are **descriptive only** — the
$1/N_e^{\text{app}}$ identity is not invertible on a mixed-depth cohort, so no
arithmetic, harmonic or depth-weighted summary of them reproduces the plug-in
fit. Do not try to back-correct a reported $N_e$ by hand.

#### Alternative discrete model

`--model discrete` retains the older hard-bottleneck model for specialised
fixed-inoculum experiments:

```txt
  Maternal true VAF p0 ~ Beta(alpha = m_alt + 1,
                              beta  = m_ref + 1)        # uniform prior
  Bottleneck count   k  ~ Binomial(Ne, p0)
  Child true VAF     p1 = k / Ne
  Child read counts  c_alt ~ Binomial(c_dp, p1)
```

Integrating `p0` analytically yields `k ~ BetaBinomial(Ne, alpha, beta)`,
so the per-pair likelihood collapses to a finite sum over `k = 0..Ne`:

```txt
  L(Ne) = Sigma_{k=0..Ne} BetaBinom(k | Ne, alpha, beta)
                          * Binom(c_alt | c_dp, k/Ne)
```

The discrete likelihood is scanned on integers only because its grid can
produce secondary local maxima. It does not support `HAS_G=1` trio rows,
and `--per-family` currently requires the continuous model.

### Three-generation G-M-C trios (`--gm-fam`)

When `mitoquest trans-prep` is run with `--gm-fam`, the output TSV
carries a `HAS_G` column: `HAS_G = 1` rows are grandmother–mother–
child (G-M-C) trios and `HAS_G = 0` rows are standard mother–child
(MC) pairs. `mitoquest ne-estimate` auto-detects the `HAS_G`,
`GRANDMOTHER_DP`, and `GRANDMOTHER_AD_ALT` columns by name and
switches each trio row to a three-generation marginal likelihood.

For a G-M-C trio the grandmother’s observed VAF
`p̂_G = g_ad_alt / g_dp` serves as the founder allele frequency.
The mother’s latent heteroplasmy `p_M` follows the same working Beta
transition approximation, `Beta(p̂_G·(Ne−1), (1−p̂_G)·(Ne−1))`. Given `p_M`, the
mother’s read count `k_M ~ Bin(d_M, p_M)` and the child’s read
count `k_C ~ BetaBin(d_C, p_M·(Ne−1), (1−p_M)·(Ne−1))` are
conditionally independent. After analytically marginalising the
child frequency, the remaining mother-frequency marginal is:

```txt
I_trio(Ne) = C(d_M, k_M) · C(d_C, k_C) / B(α_G, β_G)
           · integral_0^1 [
               p_M^(α_G+k_M-1) · (1-p_M)^(β_G+d_M-k_M-1)
               · B((Ne-1)p_M+k_C, (Ne-1)(1-p_M)+d_C-k_C)
               / B((Ne-1)p_M, (Ne-1)(1-p_M))
             ] dp_M

where  α_G = p̂_G · (Ne − 1),   β_G = (1 − p̂_G) · (Ne − 1)
```

This is not a single conjugate Beta-function ratio: the child
Beta-Binomial factor depends on `p_M` through rising factorials.
`ne-estimate` evaluates the integral with adaptive Gauss-Jacobi
quadrature, using the grandmother-plus-maternal-read posterior Beta density
as its weight.
The implementation begins at 16 nodes, refines until stable, and caps
the order at the polynomial-exact value `ceil((d_C + 1) / 2)`.
Low-depth finite-sum identities and 80-decimal `mpmath` calculations
are used in the test suite as independent references.

The composite log-likelihood sums these trio contributions with the
standard 2-gen Beta-Binomial row likelihoods. For `HAS_G = 0`, the
code uses exactly the 2-gen path. Homoplasmic grandmothers and `Ne = 1`
are evaluated as the corresponding degenerate G-M-C probabilities,
rather than discarded as uninformative. The discrete model rejects
declared trio rows. As with other composite
likelihood analyses, independence and neutral-drift assumptions should
be assessed before interpreting profile intervals as calibrated
frequentist confidence intervals.

Example:

```bash
# MC FAM: mother-child pairs
mitoquest trans-prep \
    -v cohort.vcf.gz \
    -f mc_pairs.fam \
  -g gm_pairs.fam \
    -o trio_pairs.tsv

# ne-estimate auto-detects HAS_G and dispatches accordingly
mitoquest ne-estimate \
    -i trio_pairs.tsv \
    --min-ne 2 --max-ne 100 \
    --model continuous \
    -o ne_result.json
```

### Why DP/AD instead of the VCF AF field?

A reasonable instinct is to take the per-sample `AF` (or `HF`) field
that `mitoquest caller` already writes for every site and feed those
values into the bottleneck model. We deliberately do not, for two
reasons:

1. `AF` is a **point estimate** that throws away depth-dependent
   uncertainty. A homoplasmic site at depth 50 (`AD = 50/50`,
   `AF = 1.0`) and at depth 5,000 (`AD = 5000/5000`, `AF = 1.0`) carry
   very different amounts of information about the maternal
   heteroplasmy. The Beta-Binomial step propagates that uncertainty
   into the per-pair likelihood.
2. The bottleneck model is exquisitely sensitive at the boundary
   (small `p0`). Plugging a noisy point estimate biases `Ne` downward
   for low-VAF pairs and underestimates the CI width.

The continuous default still uses the count-derived $\hat p_M$ as a plug-in
rather than treating the VCF `AF` text as a truth parameter, and it retains
the Beta-Binomial child-read likelihood. The discrete alternative additionally
marginalizes a maternal Beta posterior. Neither branch lets a stale or
non-finite `MOTHER_VAF` text value change the fitted rows: the VAF gate is
recomputed from validated `MOTHER_AD_ALT/MOTHER_DP` counts.

### Comparison with the deCODE 2024 *Cell* paper

[Helgason et al. (2024), *Cell*](https://doi.org/10.1016/j.cell.2024.05.022)
estimated the human mtDNA bottleneck on 3,457 mother-child pairs from
Iceland. Their headline number for direct M-C transmissions is
**Ne ≈ 2.29 (95% CI 1.95–2.62)**; pooling deeper-pedigree relatives
places it at **Ne ≈ 3** (best fit in the simulated range 2–30).

Methodologically, deCODE uses the Wonnapinij bottleneck parameter `b`
(based on the Kimura diffusion; Wonnapinij 2008/2010) fitted to the
*variance* of frequency changes between relatives.

Since v1.8.2, `mitoquest ne-estimate` uses a **continuous working Beta
transition MMLE** as the default model:

```txt
p_child | p_mother  ~  Beta(p_m × (Ne − 1), (1 − p_m) × (Ne − 1))
c_alt   | p_child   ~  BetaBinomial(c_dp, p_m × (Ne − 1), (1 − p_m) × (Ne − 1))
```

The continuous MMLE and the Wonnapinij/Kimura cross-check share the
one-generation drift-variance target $p(1-p)/N_e$, but they are not the
same likelihood or interchangeable estimators. Agreement is a useful
diagnostic, not proof that all model assumptions hold. The current Kimura
calculation uses the mother-child contrast on every row; it does not fit an
explicit two-transition G-M-C likelihood or a general generation-count
parameter. For declared G-M-C trios, use the continuous scorer described
above.

### Why two MMLE models?

The previous default (v1.8.0–v1.8.1) was the **discrete model**:

```txt
k       ~ BetaBinomial(Ne, α, β)
c_alt   ~ Binomial(c_dp, k/Ne)
```

This restricts the child's heteroplasmy to the grid {0, 1/Ne, …, 1}.
At high sequencing depths (DP ≥ 2000), the Binomial likelihood is an
extremely tight spike, and any mismatch between the child's true VAF
and the nearest grid point creates an enormous penalty — forcing the
MMLE to inflate Ne upward (systematic ~5–10× bias; see
`release_v1.8.2.md`).  In contrast, the continuous model allows the
child's heteroplasmy to be any value in [0, 1], correctly capturing
post-bottleneck vegetative segregation during cell division.

The discrete model is still available via `--model discrete` for
specialised use cases (e.g. virus-passage experiments where the
physical inoculum count is the target).

### Choosing the right estimator

| Estimator | Role | Strengths | Limitations |
|-----------|------|-----------|-------------|
| **Continuous MMLE** (default) | Primary | Real-valued parameter and full child count likelihood | Maternal frequency integrated out by default, which costs ~80× the plug-in wall clock (`--no-maternal-marginalization`); working Beta transition, composite-likelihood interval |
| **Kimura cross-check** (`--cross-check kimura`) | Secondary diagnostic | Sampling-error-corrected ratio of sums; fast and interpretable | Uses only corrected M-C shifts; sensitive to outliers and does not fit explicit trio transitions |
| **Discrete MMLE** (`--model discrete`) | Specialised | Finite hard-bottleneck model for fixed-inoculum experiments | Integer grid, can inflate on deep mtDNA data, and rejects declared trio rows |

**Decision rule:**

1. Use the **continuous MMLE** as the reported Ne (default behaviour).
2. Run `--cross-check kimura --kimura-trim 0.10` as a sanity check.
3. Treat agreement or CI overlap as a diagnostic consistency check, not an
  independent validation of the shared assumptions.
4. When they disagree materially, inspect outlier pairs with
   `--top-drift-k 20` — a handful of high-drift pairs (NUMTs /
  sequencing errors) can collapse the variance-based Kimura estimate.
5. For G-M-C data, use `--model continuous`; `--model discrete` rejects
  `HAS_G=1`, and the Kimura cross-check remains an M-C diagnostic.
6. If maternal or child depth is below ~500×, the default marginalised fit
   already removes the dominant plug-in bias; add `--min-depth INT` to gate the
   shallowest rows as well. That is a QC lever rather than a model fix — it
   drops rows — so the two are complementary, not interchangeable. Check
   `Mean_Mother_DP` / `Harmonic_Mother_DP` in the JSON to see which regime you
   are in. Only after opting out with `--no-maternal-marginalization` do you
   get the depth-dependent plug-in, whose estimates must not be compared across
   cohorts sequenced at different coverage.

#### Why Ne_MMLE and Ne_Kimura can differ

There is no fixed ordering between the continuous MMLE and the Kimura
moment estimate. Their difference can arise from finite-sample behavior,
outliers, input quality, or failure of the working one-generation model:

1. **Jensen's inequality bias in the Kimura ratio estimator.**
  The Kimura formula computes `Ne = 1/V` where `V = Σ(d_i − s_i) / Σw_i`.
   Because `1/x` is a convex function, the estimator `1/V̂` is biased
   *upward* by Jensen's inequality: `E[1/V̂] > 1/E[V̂]`.  With large
   cohorts (hundreds of pairs) the CIs are tight enough to expose this
   systematic upward bias in Ne_Kimura.

2. **Full distribution vs. second moment only.**
  The MMLE fits the entire Beta-Binomial count shape, while the
   Kimura estimator uses only the second central moment (variance).  At
   small Ne the Beta distribution is highly non-Gaussian (often U-shaped),
   so the higher-order information captured by the MMLE carries real signal
   that the variance-only Kimura discards.

3. **Different use of depth and count information.**
  The likelihood responds to each row's count distribution, while the
  Kimura estimate aggregates a corrected squared shift divided by the
  maternal variance scale.

4. **Plug-in sampling correction noise.**
  The implemented finite-depth correction is
  $s_i=\hat p_M(1-\hat p_M)/(d_M-1)+\hat p_C(1-\hat p_C)/(d_C-1)$.
  It excludes rows with maternal or child depth at most one, and still
  uses plug-in frequencies.

5. **The maternal plug-in — the one term that is predictable, not random.**
  The Kimura side subtracts the maternal sampling variance $s_i$, so it is
  unbiased for the drift; the MMLE side cannot identify that term, because the
  Beta-Binomial already marginalises the *child's* read noise but nothing
  marginalises the *mother's*. To first order
  $1/\widehat N_e^{\text{MMLE}}\approx 1/N_e+1/d_M$, i.e.
  $N_e^{\text{Kimura}}/\widehat N_e^{\text{MMLE}}\approx 1+N_e/d_M$.
  The offset is set entirely by the ratio $N_e/d_M$: ~0.6% at
  $N_e=3,\ d_M=500$ (invisible) but ~30% at $N_e=30,\ d_M=100$. **Before
  deciding that a disagreement is anomalous, subtract this expected offset.**
  It is a property of the plug-in, not a bug; it disappears under the default
  marginalised fit (measured: the ratio returns to 0.98–1.00 at every depth),
  so you only see it after opting out with `--no-maternal-marginalization`.
  Note it does *not* survive small samples — with only a handful
  of sites both estimators' sampling noise dwarfs the systematic term and the
  direction is random (the 5-site worked example in
  `handoff/HANDOFF_ne_estimate_methods.md` §10 has Kimura *below* MMLE).

**Is it always MMLE < Kimura?**  No — when high-drift outliers (NUMTs,
sequencing errors) dominate, they inflate V and collapse Ne_Kimura *below*
Ne_MMLE.  The direction depends on the data:

| Scenario | Typical relationship |
| :----------: | :---------------------: |
| Clean, adequately powered simulation under the working model | Often similar, but no fixed ordering is guaranteed |
| High-drift outliers present | Kimura may be pulled downward by inflated corrected drift |
| Sparse or shallow data | Either estimate can be unstable or boundary-clipped |

**Bottom line:** Use the continuous MMLE as the configured primary estimate
and treat Ne_Kimura as a qualitative diagnostic. Material disagreement is a
reason to inspect data and assumptions; it is neither automatic evidence of
a data problem nor a reason to choose an estimator solely by direction.

### Full parameter reference of `ne-estimate`

```bash
Usage: mitoquest ne-estimate [options] -i <pairs.tsv>

Required options:
  -i, --input     FILE   Input transmission pairs TSV produced by
                         `mitoquest trans-prep`.

Optional options:
  -o, --output    FILE   JSON output file (default: stdout).
      --model     NAME   Marginal-likelihood model: `continuous` (default,
                         recommended for mtDNA) or `discrete`.
      --no-maternal-marginalization
                         Fall back to the legacy plug-in p_m = m_alt /
                         m_dp for the continuous model's two-generation
                         rows.  By default p_m is integrated out against
                         its Beta(m_alt+1, m_ref+1) read-sampling
                         posterior by adaptive Gauss-Jacobi quadrature;
                         the plug-in is faster but biased low by
                         ~d_M/(Ne+d_M) and is not comparable across
                         coverage, so it is available only by explicit
                         request.  Continuous model only: --model discrete
                         always integrates p_m and cannot honour this flag
                         (warned on stderr); trio rows also already
                         integrate p_m given the grandmother, so it is a
                         no-op for them.
      --min-vaf   FLOAT  Lower maternal VAF gate, inclusive [0.10].
      --max-vaf   FLOAT  Upper maternal VAF gate, inclusive [0.90].
      --min-depth INT    Drop pairs whose MOTHER or CHILD depth is below
                         this value (both sides are gated: a shallow child
                         makes k_C uninformative just as a shallow mother
                         makes p_m unreliable).  0 disables [0].
      --min-ne    INT    Smallest Ne value to consider [1].
      --max-ne    INT    Largest Ne value to consider  [100].
  -t, --threads   INT    Worker threads for the inner sum [1].
      --cross-check NAME   Optional secondary estimator alongside the
                           MMLE. Supported value: `kimura`, which
                           computes the Wonnapinij b and the implied
                           Ne (single-generation approximation).
      --kimura-bootstrap  INT  Non-parametric bootstrap iterations for the
                              Kimura cross-check 95% CI.  0 disables [1000].
      --kimura-seed      INT  RNG seed for the Kimura bootstrap [42].
      --kimura-trim    FLOAT  Fraction of highest-drift pairs to drop before
                              recomputing the Kimura b (robust trimmed mean,
                              0.0 disables, recommended 0.10) [0.0].
      --top-drift-k     INT   Emit the top-K highest-drift pairs in the JSON
                              output for outlier inspection (NUMTs / errors).
                              0 disables [0].
      --bin-simulation FILE   Emit a per-bin drift summary TSV: per
                              maternal-VAF bin observed mean drift vs
                              analytical Kimura prediction p_m(1 - p_m) / Ne
                              under the fitted Ne (and its 95% CI).  Bin
                              range = [--min-vaf, --max-vaf].
      --bin-simulation-bins INT  Number of equal-width maternal-VAF bins
                                 for --bin-simulation [10].
      --ne-profile     FILE   Emit a TSV that scores every candidate Ne
                              under both the MMLE marginal log-likelihood
                              and the Kimura per-pair SSR metric.  Useful to
                              visually compare which Ne each estimator
                              prefers (dual-objective Ne scan).
                              Grid range = [--min-ne, --max-ne].
      --ne-profile-step FLOAT Grid step on the Ne axis for --ne-profile
                              [0.1].
      --per-family              Enable per-family Ne estimation.
                                Groups pairs by FAM_ID + MOTHER_ID and
                                estimates Ne independently for each family.
      --min-family-sites INT    Minimum informative sites per family [3].
                                Families with fewer sites are skipped.
      --per-family-output FILE  Write per-family Ne results as TSV (one
                                row per family).  If not set, per-family
                                results are embedded in the JSON output.
  -h, --help               Print this help message.
```

### Output JSON

```json
{
  "Ne":              2.78,
  "CI_95_Low":       2.12,
  "CI_95_High":      3.54,
  "CI_Low_Clipped":  false,
  "CI_High_Clipped": false,
  "Pairs_Used":      147,
  "Trio_Founder_Mismatch_Skipped": 0,
  "Trio_Founder_Hom_Skipped":      0,
  "Min_Depth_Skipped":             0,
  "Estimator":       "MMLE (composite marginal likelihood)",
  "Max_Marginal_LogLik": -812.34561,
  "Model":           "continuous",
  "Min_VAF":         0.10,
  "Max_VAF":         0.90,
  "Search_Min_Ne":   1,
  "Search_Max_Ne":   100,
  "Min_Depth":       0,
  "Maternal_Marginalization": false,
  "Mean_Mother_DP":     1834.2,
  "Mean_Child_DP":      1791.6,
  "Harmonic_Mother_DP": 1602.8,
  "Harmonic_Child_DP":  1577.3,
  "Kimura_Cross_Check": {
    "b":             0.732,
    "Ne_Kimura":     3.731,
    "b_CI_95_Low":   0.614,
    "b_CI_95_High":  0.819,
    "Ne_Kimura_CI_95_Low":  2.59,
    "Ne_Kimura_CI_95_High": 5.52,
    "N_Bootstrap":   1000,
    "Bootstrap_Seed": 42,
    "N_Informative": 145,
    "Trimmed_Kimura": {
      "Trim_Frac":        0.1,
      "N_After_Trim":     130,
      "b_Trimmed":        0.753,
      "Ne_Kimura_Trimmed": 4.05
    },
    "Top_Drift_Outliers": [
      { "Pair_Index": 42, "M_DP": 2000, "M_AD_ALT": 1000,
        "C_DP": 2000, "C_AD_ALT": 0, "M_VAF": 0.5,
        "C_VAF": 0.0, "F_i": 0.99 }
    ],
    "Note":          "",
    "Method":        "Wonnapinij 2008/2010 with sampling-error correction; single-generation Ne = 1 / (1 - b); 95% CI by non-parametric pair-level bootstrap"
  }
}
```

The `Kimura_Cross_Check` block is only emitted when
`--cross-check kimura` is passed.  The `Trimmed_Kimura` and
`Top_Drift_Outliers` sub-blocks appear only when `--kimura-trim > 0`
and `--top-drift-k > 0`, respectively.

- `CI_Low_Clipped` / `CI_High_Clipped` flag confidence-interval bounds
  that hit the search boundary (`--min-ne` / `--max-ne`); when either is
  `true` you should re-run with a wider search range.
- `Trio_Founder_*_Skipped` are non-zero only on input carrying three-generation
  rows (`HAS_G=1`); `Min_Depth_Skipped` counts rows dropped by `--min-depth`.
- `Min_Depth` and `Maternal_Marginalization` echo the model-form switches
  **actually applied**, so a reader can tell which likelihood the reported
  `Max_Marginal_LogLik` belongs to. `Maternal_Marginalization` is `true` for a
  default continuous run, `false` when you opted out with
  `--no-maternal-marginalization`, and also `false` under `--model discrete` —
  that model always integrates $p_M$ and therefore cannot honour the opt-out,
  which it reports with a `WARNING:` rather than silently.
  `Max_Marginal_LogLik` is **not comparable** across that flag.
- `Mean_*_DP` and `Harmonic_*_DP` describe the retained rows' coverage and are
  **descriptive only**. The plug-in pseudo-true identity
  $1/N_e^{\text{app}}=1/N_e+1/d_M$ holds at a single depth and is not
  invertible on a mixed-depth cohort, so no arithmetic, harmonic or
  depth-weighted summary of these fields reproduces the plug-in fit. The
  default marginalised fit removes that bias at the model level, so neither
  back-correct a reported $N_e$ by hand nor try to repair a plug-in run
  arithmetically — drop the opt-out flag instead.
- Per-family records (`--per-family`) carry their own `Mean_*_DP` /
  `Harmonic_*_DP`. The three TSV outputs (`--per-family-output`,
  `--bin-simulation`, `--ne-profile`) keep their previous columns; the new
  fields appear in JSON only.
- Every invocation writes a leading provenance comment
  (`#mitoquest_ne_estimate_command=...`) before the JSON object. After that
  one line is removed, the object is strict JSON: dynamic strings are
  escaped and any non-finite optional diagnostic is represented as a quoted
  string such as `"Infinity"` or `"NaN"`, never as a bare JSON token.

### Usage example for `trans-prep`

End-to-end pipeline chained with `trans-prep`:

```bash
mitoquest trans-prep \
    -v cohort.vcf.gz \
    -f cohort.fam \
    -d 500 \
    -o cohort.transmission_pairs.tsv

mitoquest ne-estimate \
    -i cohort.transmission_pairs.tsv \
    --cross-check kimura \
    --kimura-trim 0.10 \
    --top-drift-k 20 \
    --min-vaf 0.10 --max-vaf 0.90 \
    --min-ne 1 --max-ne 100 \
    -t 8 \
    -o cohort.ne.json
```

### Diagnostic TSV outputs

Two optional TSV outputs let you visually assess and compare the MMLE
and Kimura fits:

| Flag                        | Purpose                                                                                       |
| --------------------------- | --------------------------------------------------------------------------------------------- |
| `--bin-simulation FILE`     | Per maternal-VAF bin observed drift vs. the analytical Kimura prediction `p_m(1 − p_m) / Ne`. |
| `--bin-simulation-bins INT` | Number of equal-width bins across `[--min-vaf, --max-vaf]` (default: 10).                      |
| `--ne-profile FILE`         | Per-candidate-Ne grid of MMLE marginal log-likelihood and Kimura SSR — lets you plot both objectives on the same axis and see which Ne each prefers. |
| `--ne-profile-step FLOAT`   | Grid step on the Ne axis (default: 0.1).                                                       |

```bash
# Emit both diagnostic TSVs alongside the JSON estimate.
mitoquest ne-estimate \
    -i cohort.transmission_pairs.tsv \
    --cross-check kimura \
    --min-ne 1 --max-ne 30 \
    --bin-simulation cohort.bin_sim.tsv \
    --bin-simulation-bins 10 \
    --ne-profile cohort.ne_profile.tsv \
    --ne-profile-step 0.1 \
    -o cohort.ne.json
```

Both TSVs start with `#`-prefixed metadata comments that record the
fitted `Ne`, its 95% CI, and the exact CLI invocation, so downstream
plotting scripts (e.g., `tools/plot_bottleneck_simulation.py`) can
reproduce the figure without re-parsing the JSON.

### Per-family Ne estimation (`--per-family`)

In addition to the population-level Ne estimate, `mitoquest ne-estimate`
can estimate the bottleneck size **independently for each mother-child
family** in the cohort. This is enabled by the `--per-family` flag.

Per-family mode groups the transmission pairs by `(FAM_ID, MOTHER_ID)`
(columns emitted by `mitoquest trans-prep`) and runs the continuous
working-Beta MMLE on each group separately. It requires `--model continuous`.
The results are reported as:

- A **per-family TSV** (`--per-family-output FILE`) with one row per
  family, including: family identifiers, number of children, number of
  pairs, number of informative sites, Ne estimate with 95% CI, mean
  depth, optional Kimura cross-check columns, and, when bootstrapping is
  enabled, Kimura percentile-CI columns.
- A `Per_Family_Estimates` array and `Per_Family_Summary` block in the
  main JSON output.

Families with fewer than `--min-family-sites` (default: 3) informative
sites are **skipped** and flagged with `SKIPPED=TRUE` and a warning
message. This prevents unreliable estimates from very small families.

> **Note:** Per-family Ne estimation on a single family with only a
> handful of variant sites will typically produce a very wide CI. This
> is expected — the working-Beta MMLE is statistically consistent but
> needs sufficient heteroplasmic sites to resolve Ne precisely.

The point estimate remains a row-summed composite likelihood. For its
family-specific Kimura CI, observations sharing the same `CHROM`, `POS`, and
`ALT` are resampled together as a locus cluster; legacy inputs without all
three fields fall back to row clusters and are annotated in the Kimura note.

#### Per-family usage example

```bash
# Population-level + per-family Ne estimation with Kimura cross-check.
mitoquest ne-estimate \
    -i cohort.transmission_pairs.tsv \
    --per-family \
    --min-family-sites 5 \
    --per-family-output cohort.per_family_ne.tsv \
    --cross-check kimura \
    --min-vaf 0.05 --max-vaf 0.95 \
    --min-ne 1 --max-ne 200 \
    -t 8 \
    -o cohort.ne.json
```

#### Per-family TSV format

The `--per-family-output` TSV has one row per family:

| Column             | Description                                                    |
| ------------------ | -------------------------------------------------------------- |
| `FAM_ID`           | Family ID from the input TSV                                   |
| `MOTHER_ID`        | Mother sample ID                                               |
| `N_CHILDREN`       | Number of distinct children in this family                     |
| `N_PAIRS`          | Total variant sites (rows) for this family                     |
| `N_INFORMATIVE`    | Sites with 0 < mother VAF < 1                                  |
| `NE_MMLE`          | Per-family Ne estimate (continuous MMLE); `NA` if skipped     |
| `CI_95_LOW`        | 95% CI lower bound; `NA` if skipped                            |
| `CI_95_HIGH`       | 95% CI upper bound; `NA` if skipped                            |
| `CI_LOW_CLIPPED`   | `TRUE` if CI lower bound hit the search boundary               |
| `CI_HIGH_CLIPPED`  | `TRUE` if CI upper bound hit the search boundary               |
| `MAX_LOG_LIK`      | Maximum marginal log-likelihood                                |
| `MEAN_MOTHER_DP`   | Mean mother sequencing depth across sites                      |
| `MEAN_CHILD_DP`    | Mean child sequencing depth across sites                       |
| `KIMURA_b`         | Kimura bottleneck parameter *b* (only when `--cross-check kimura`) |
| `KIMURA_NE`        | Kimura-implied Ne (only when `--cross-check kimura`)          |
| `KIMURA_N_INFORMATIVE` | Kimura informative pair count                             |
| `KIMURA_B_CI_95_LOW`, `KIMURA_B_CI_95_HIGH` | Bootstrap percentile interval for $b$ (when computed) |
| `KIMURA_NE_CI_95_LOW`, `KIMURA_NE_CI_95_HIGH` | Bootstrap percentile interval for Kimura Ne (when computed) |
| `KIMURA_N_BOOTSTRAP`, `KIMURA_BOOTSTRAP_SEED` | Bootstrap settings used for the family CI |
| `SKIPPED`          | `TRUE` when family has too few informative sites               |
| `WARNING`          | Free-text warning (e.g. "small sample", "CI clipped")          |

The TSV file starts with `#`-prefixed metadata comments recording the
CLI parameters, the population-level Ne estimate, and the number of
families estimated / skipped.

### Validating `ne-estimate` without real data

The `tests/data/ne_pipeline/` directory ships a complete synthetic
cohort that lets you validate the entire `trans-prep → ne-estimate`
pipeline end-to-end **without any real sequencing data**.

#### Quick validation (population-level Ne)

```bash
# Run the bundled smoke-test script (true Ne = 5, 20 trios, seed 20260528):
bash tests/data/ne_pipeline/run_demo.sh

# Expected output:
#   Ne ≈ 3–8  (true Ne = 5; MMLE is consistent but stochastic at N=20)
#   Ne_kimura in the same ballpark
#   Both CI_Low_Clipped and CI_High_Clipped = false
```

The synthetic cohort was generated by `tests/data/ne_pipeline/synthesize.py`
under the continuous Beta-diffusion model with known `true_ne = 5`. The
MMLE should recover a value within ~50% of the truth on this small
cohort; larger cohorts (hundreds of pairs) yield tighter estimates.

#### Quick validation (per-family Ne)

```bash
# Per-family test with 3 synthetic families:
#   F001 (2 children, 15 sites), F002 (1 child, 20 sites), F003 (2 sites, should be skipped)
./bin/mitoquest ne-estimate \
    -i tests/data/ne_pipeline/per_family_validation.tsv \
    --per-family \
    --per-family-output /tmp/per_family_ne.tsv \
    --cross-check kimura \
    --min-vaf 0.05 --max-vaf 0.95 \
    -o /tmp/pop_ne.json

# Expected:
#   F001: Ne estimated (15 informative sites, 2 children)
#   F002: Ne estimated (20 informative sites, 1 child)
#   F003: SKIPPED (only 2 sites < min_family_sites=3)
cat /tmp/per_family_ne.tsv
```

#### Running the full unit-test suite

```bash
# Build with testing enabled, then run all tests:
cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel
./bin/mitoquest_tests --gtest_filter="NeEst*:NeEstFamily*"

# Expected: 56 tests, 0 failures
# (NeEstFamily suite: 23 tests covering grouping, estimation,
#  serial/parallel consistency, Kimura, TSV output, and end-to-end)
```

#### Additional validation scenarios

The `tests/data/ne_pipeline/scenarios/` directory contains 11
pre-generated scenarios (`S1`–`S11`) with known Ne values ranging from
3 to 200, covering: baseline, fractional Ne, small cohorts, low depth,
outliers, large Ne, and multi-generation pedigrees. Each scenario
includes both continuous and discrete model results for comparison.

```bash
# Spot-check one scenario (S8: true Ne = 10, 2 generations):
./bin/mitoquest ne-estimate \
    -i tests/data/ne_pipeline/scenarios/S8_g2_ne10/ne.continuous.json 2>/dev/null || \
cat tests/data/ne_pipeline/scenarios/S8_g2_ne10/ne.continuous.json
```

### Recommendations

- **Sequencing depth ≥ 500×** on chrM for both mother and child. The
  dominant risk at low depth is **bias, not just a wider CI**, and it is the
  reason marginalisation is the default: the legacy plug-in maternal frequency
  makes the fitted value a pseudo-true $N_e^{\text{app}}$ with
  $1/N_e^{\text{app}}=1/N_e+1/d_M$, measured at $-13.5\%$ for $d=100$ and
  $-2.9\%$ for $d=500$ ($N_e^{\text{true}}=20$). The default fit removes that
  term and leaves $+14.1\%$ at $d_M=30$, $\le 2.4\%$ in magnitude at every
  $d_M\ge 60$, and $-2.1\ldots+2.4\%$ for $d_M\ge 500$ — within the
  $1.3\ldots2.5\%$ Monte Carlo se of zero, with the sign varying rather than
  staying positive. Check `Mean_Mother_DP` /
  `Harmonic_Mother_DP` in the JSON output; below ~500× also add
  `--min-depth INT` to drop the shallowest rows. If you do opt out with
  `--no-maternal-marginalization`, its estimates are depth-dependent and must
  not be compared across cohorts sequenced at different coverage.
- **Restrict to maternal VAF in [0.10, 0.90]** (the default). Sites very
  close to 0 or 1 carry virtually no information about Ne, while still
  contributing numerical noise.
- **At least ~30 informative pairs** are needed for a stable CI; with
  fewer pairs the CI tends to be wide and to clip the search boundaries.
- **Verify both `_Clipped` flags are `false`**; if either is `true`,
  widen `--min-ne` / `--max-ne` and re-run.
- **Treat the 95% interval as model-based.** Pairs are summed as a composite
  likelihood and the $\chi^2_1$ threshold 1.92 is applied without a sandwich
  correction, so the interval can be too narrow when many siblings share one
  mother (more so under the default marginalised fit). Use `--per-family` for
  family-resolved inference in that regime.

---

## `mitoquest variant-qc` — Bayesian quality control for mtDNA variants

`mitoquest variant-qc` performs per-site, per-sample quality control on a
multi-sample mtDNA VCF (typically produced by `mitoquest caller`) to
distinguish true heteroplasmic mutations from sequencing artefacts.

### Statistical model

The module implements **Bayesian hypothesis testing** with
**Beta-Binomial** likelihood models:

```txt
  H0 (background noise):  alt_count ~ BetaBinomial(D, q_alpha, q_beta)
  H1 (true mutation):     alt_count ~ BetaBinomial(D, alpha_h1, beta_h1)
```

- Per-site background noise parameters `(q_alpha, q_beta)` are estimated
  via Beta-MLE from homozygous-reference samples at each site.
- Global true-mutation parameters `(alpha_h1, beta_h1)` are estimated via
  Beta-Binomial MLE from high-confidence variant observations.
- The posterior probability `P(H1 | data)` is computed via Bayes' theorem
  with a prior `pi` and quality penalties (strand-bias SRF and allele
  quality HQ sigmoid).
- An EM-style iterative refinement loop (fit → call → refit) repeats
  until mutation calls converge.

The output VCF gains `FORMAT/PP` (posterior probability) and
`FORMAT/GOOD_CALL` (True/False) per sample, plus site-level `FILTER`
flags (`PASS`, `LOW_QUALITY`, `BLACKLISTED_SITE`).

### Full parameter reference of `variant-qc`

```bash
Usage: mitoquest variant-qc [options]

Required:
  -i, --input-vcf FILE       Input indexed VCF/BCF file
  -o, --output-vcf FILE      Output filtered VCF file
  -t, --output-tsv FILE      Output tabular report (TSV)

Optional:
  --max-alt-alleles INT      Max ALT alleles per site [2]
  --dp-threshold INT         Min DP for pre-filter [100]
  --hq-threshold INT         Min allele quality [20]
  --bins INT                 Histogram bins for KL div [100]
  --pi FLOAT                 Prior prob of mutation [5e-8*16569]
  --threshold FLOAT          Posterior threshold [0.9]
  --max-iter INT             Max EM iterations [20]
  --convergence-eps FLOAT    Convergence threshold [0.001]
  -h, --help                 Show this help
```

### Usage examples for `variant-qc`

**Basic quality control on a cohort VCF:**

```bash
mitoquest variant-qc \
    -i cohort.raw.vcf.gz \
    -o cohort.qc.vcf.gz \
    -t cohort.qc.report.tsv
```

**Stricter QC (higher depth and quality thresholds):**

```bash
mitoquest variant-qc \
    -i cohort.raw.vcf.gz \
    -o cohort.qc.vcf.gz \
    -t cohort.qc.report.tsv \
    --dp-threshold 200 --hq-threshold 30 \
    --threshold 0.95
```

### Recommendations

- **Input VCF must be indexed** (`.tbi` or `.csi`).
- **Higher depth is better**: the Beta-Binomial model benefits from
  high coverage (≥ 500× on chrM) for accurate background estimation.
- The TSV report contains one row per (sample, site) flagged as a
  mutation, with posterior probability and fitted Beta parameters.

--- 
-->


## Auxiliary Python tools (`tools/`)

The `tools/` directory ships a set of Python helper scripts that complement
the C++ binaries:

| Script                          | Purpose                                                                                  |
| ------------------------------- | ---------------------------------------------------------------------------------------- |
| `tools/mito_annotate.py`        | Annotate a `mitoquest caller` VCF with population, in-silico, and disease databases.     |
| `tools/parse_annotatedVCF.py`   | Convert an annotated VCF into a tidy table for downstream statistics.                    |
| `tools/mtDNA_variant_QC.py`     | Python QC prototype (superseded by `mitoquest variant-qc`).                              |
| `tools/mtDNA_vcf_to_table.py`   | Flatten a multi-sample VCF into long-format TSV (one row per sample × site).             |
| `tools/parse_vcf.py`            | Generic VCF parsing helper (used by other tools).                                        |
| `tools/filter_mergedVCF.py`     | Apply hard filters on a merged mtDNA VCF.                                                |
| `tools/merge_cr_ncr_vcf.py`     | Merge coding-region and non-coding-region VCFs into a single file.                       |
| `tools/rewrite_vcf.py`          | Rewrite a VCF (header normalisation, CHROM renaming, etc.).                              |
| `tools/shift_fasta.py`          | Produce a circularly shifted FASTA (used by the join-region pipeline).                   |
| `tools/create_join_seq.py`      | Build the joined coding-region / non-coding-region reference for re-alignment.           |
| `tools/detect_NUMT_by_mtCN.py`  | Flag potential NUMT contamination using mtCN ratios per sample.                          |
| `tools/plot_bottleneck_simulation.py` | Plot per-bin observed drift vs. analytical Kimura prediction `p_m(1−p_m)/Ne` from `ne-estimate --bin-simulation`. |
| `tools/plot_ne_profile.py`      | Plot the MMLE and Kimura objective curves over Ne from `ne-estimate --ne-profile`.        |
| `tools/vcf_format_validator.py` | Sanity-check a VCF for downstream compatibility.                                         |

Each script supports `-h / --help`.

### A typical end-to-end workflow

```bash
# 1. Variant calling
mitoquest caller \
    -f rCRS.fasta \
    -o cohort.raw.vcf.gz \
    -Q 30 -q 30 -t 24 \
    --filename-has-samplename \
    -b bamfile.list

# 2. Annotation against MITOMAP / HelixMTdb / gnomAD / in-silico predictors
python tools/mito_annotate.py \
    -i cohort.raw.vcf.gz \
    -o cohort.annotated.vcf.gz \
    --resource-dir data

# 3. Bayesian quality control (filter true mutations from artefacts)
mitoquest variant-qc \
    -i cohort.raw.vcf.gz \
    -o cohort.qc.vcf.gz \
    -t cohort.qc.report.tsv \
    --dp-threshold 200 --hq-threshold 30

# 4. mtDNA copy-number estimation per sample
while IFS= read -r bam; do
    mitoquest copynum -r reference.fasta -q 30 -t 4 "$bam"
done < bamfile.list > cohort.mtCN.tsv

# 5. Convert to a long-format table for analysis in R/Python
python tools/mtDNA_vcf_to_table.py \
    -i cohort.qc.vcf.gz \
    -o cohort.long.tsv
```

---

## Reference resources (`data/`)

The `data/` directory contains curated reference resources used by the
annotation and QC tools. See [data/README.md](data/README.md) and
[data/UPDATE_LOG.md](data/UPDATE_LOG.md) for sources and update history.

Highlights:

- **Population databases**: HelixMTdb, gnomAD v3.1 (chrM), dbSNP (chrM),
  MITOMAP polymorphisms / disease, ClinVar (chrM SNVs).
- **In-silico predictors**: MitoTIP, MitImpact, HmtVar, t-APOGEE,
  phyloP100way.
- **Functional annotations**: tRNA secondary-structure positions, rRNA
  bridge bases, RNA modifications, UniProt protein domains, Complex-I
  proton-pump residues.
- **Phylogenetic context**: Phylotree variants, human-shifted chimp
  alignment, mtDNA gene loci.
- **Blacklist**: known systematic-error sites (`data/blacklist.txt`).

---

## Tips and best practices

- **Reference FASTA**: Use the **rCRS** sequence (`NC_012920.1`) as the
  reference whenever possible — most downstream databases and tools
  (MITOMAP, HelixMTdb, MitoTIP) are coordinated to it.
- **`--filename-has-samplename`**: If your BAM files are named
  `{SampleID}.bam` or `{SampleID}.cram`, always set this flag — it avoids
  reading every BAM header and can save hours on large cohorts.
- **`-j / --het-threshold`**: The default `0.01` (1 %) is a sensible balance
  between sensitivity and false-positive rate at typical WGS depths (≥30×
  on chrM). Lower it to `0.005` only when coverage is very high (≥1000×)
  and you are explicitly targeting low-frequency heteroplasmy.
- **`-r / --regions`**: For large multi-sample analyses, you can split
  chrM into chunks and run several `mitoquest caller` instances in
  parallel; each instance is multithreaded internally via `-t`.
- **Output compression**: Always use `.vcf.gz` as the output filename —
  `mitoquest` automatically writes bgzipped, tabix-ready output when the
  extension matches.
- **NUMT filtering**: For samples with elevated mtCN ratios (`mitoquest copynum`
  output much higher than expected), consider either:
    * restricting the mtDNA measurement to NUMT-poor stretches with
      `mitoquest copynum -L 'chrM:1-300,chrM:16000-16569' ...`
      (chrM is measured only on the listed intervals; autosomes stay full-length),
    * and/or running the caller with `-P/--proper-pairs-only`
      and inspecting outputs with `tools/detect_NUMT_by_mtCN.py`.

---

## Development

MitoQuest is under active development. To update to the latest version:

```bash
git pull
git submodule update --recursive
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

To rebuild from scratch (recommended after a major htslib update):

```bash
rm -rf build
(cd htslib && make distclean) || true
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

Pull requests, bug reports and feature requests are very welcome at
<https://github.com/ShujiaHuang/mitoquest>.

> *"What I cannot create, I do not understand." — Richard Feynman*
