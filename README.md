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
Version: 1.9.1

Usage: mitoquest <command> [options]
Commands:
  caller       Mitochondrial variants and heteroplasmy/homoplasmy caller.
  subsam       Extract mitochondrial variants for specified samples from VCF files and output a new VCF file.
  copynum      Estimate per-chromosome (incl. mtDNA) relative copy number from a BAM/CRAM file.
  trans-prep   Extract mother-child mtDNA allele transmission pairs from a multi-sample VCF + FAM file.
  ne-estimate  Estimate the mtDNA bottleneck size (Ne) from transmission pairs via Beta-Binomial MMLE (Maximum Marginal Likelihood).
  variant-qc   Bayesian quality control for mtDNA variants from VCF files.

```

In addition to the main `mitoquest` binary, the project ships:

- `tools/` — a suite of Python helper scripts for VCF QC, annotation,
  pipeline assembly, and downstream analysis.
- `data/` — curated reference resources (population databases, in-silico
  predictors, RNA/protein domain annotations, blacklist, Phylotree
  variants, …) used by the annotation and QC tools.

---

## Installation

### Option 1 — Download pre-built binary (Recommended, no compilation needed)

Pre-built static binaries are published on the
[GitHub Releases page](https://github.com/ShujiaHuang/mitoquest/releases).

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

**Confirmed compatible distributions** (glibc ≥ 2.35):

| Distribution            | glibc      | `mitoquest-linux-static` |
| ----------------------- | ---------- | :----------------------: |
| Ubuntu 22.04 LTS        | 2.35       | ✅                       |
| Ubuntu 24.04 LTS        | 2.39       | ✅                       |
| Debian 12 (bookworm)    | 2.36       | ✅                       |
| Fedora 36+              | 2.35+      | ✅                       |
| openSUSE Tumbleweed     | rolling    | ✅                       |

**Distributions where `mitoquest-linux-static` will NOT run** (glibc too old —
please [compile from source](#option-2--compile-from-source) instead):

| Distribution                                  | glibc          | `mitoquest-linux-static` |
| --------------------------------------------- | -------------- | :----------------------: |
| CentOS 7 / RHEL 7                             | 2.17           | ❌                       |
| CentOS 8 / RHEL 8 / Rocky 8 / Alma 8          | 2.28           | ❌                       |
| CentOS 9 / RHEL 9 / Rocky 9 / Alma 9          | 2.34           | ❌                       |
| Ubuntu 18.04 / 20.04                          | 2.27 / 2.31    | ❌                       |
| Debian 10 / 11                                | 2.28 / 2.31    | ❌                       |

**Quick check on your machine:**

```bash
# If the printed glibc version is >= 2.35, mitoquest-linux-static will run.
ldd --version | head -1
```

A typical incompatibility error looks like:

```
./mitoquest-linux-static: /lib64/libc.so.6: version `GLIBC_2.35' not found
(required by ./mitoquest-linux-static)
```

If you see this — or you are on CentOS / RHEL / Rocky / AlmaLinux / older
Ubuntu / older Debian — please use [Option 2: compile from source](#option-2--compile-from-source).
The build is straightforward and takes only a few minutes.

```bash
# Linux
wget https://github.com/ShujiaHuang/mitoquest/releases/latest/download/mitoquest-linux-static
chmod +x mitoquest-linux-static
./mitoquest-linux-static --help
```

```bash
# macOS
curl -LO https://github.com/ShujiaHuang/mitoquest/releases/latest/download/mitoquest-macos-static
chmod +x mitoquest-macos-static
./mitoquest-macos-static --help
```

#### System requirements for `mitoquest-macos-static`

The macOS binary is a **best-effort static** build. Apple does not support
fully-static executables, so only Apple system frameworks (`libSystem`,
`libc++abi`, `libcurl`, …) remain dynamic; everything else (htslib, zlib,
bzip2, xz) is statically linked.

- Tested on **macOS 13+** on Apple Silicon (arm64) — should also run on
  Intel macOS 12+.
- If you hit a "code signature invalid" error after `chmod +x`, run
  `xattr -d com.apple.quarantine ./mitoquest-macos-static` once to clear
  Gatekeeper's quarantine flag.

---

### Option 2 — Compile from source

*Requirements: C++17 compiler (GCC 7+ or Apple Clang 10+), CMake ≥ 3.12, and
the system libraries: zlib, bzip2, xz-utils, libcurl, openssl (Linux only).*

#### Step 1 — Clone the repository (including the htslib submodule)

```bash
git clone --recursive https://github.com/ShujiaHuang/mitoquest.git
cd mitoquest
```

> If you forgot `--recursive`, run: `git submodule update --init --recursive`

#### Step 2 — Build with CMake (standard dynamic build)

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

The executable `bin/mitoquest` will be produced. Verify with:

```bash
./bin/mitoquest --help
./bin/mitoquest copynum --help
```

#### Step 3 (Optional) — Run the unit tests

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel
cd build && ctest --output-on-failure
```

> Requires a system GoogleTest installation (e.g., `brew install googletest`
> on macOS or `sudo apt-get install libgtest-dev` on Ubuntu 22.04+).

#### Step 4 (Optional) — Build a static binary locally

**Linux** (portable static via Ubuntu/glibc — same approach used in CI):

```bash
sudo apt-get install -y build-essential cmake autoconf automake \
    zlib1g-dev libbz2-dev liblzma-dev libssl-dev
cmake -B build-static -DSTATIC_BUILD=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-static --parallel
```

This bundles `libstdc++`, `libgcc`, `htslib`, `openssl`, and the compression
libraries statically; glibc remains dynamic. The resulting binary runs on
the build host and on any Linux with the same-or-newer glibc.

**macOS** (Homebrew):

```bash
brew install autoconf automake zlib bzip2 xz curl
cmake -B build-static -DSTATIC_BUILD=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-static --parallel
```

---

### Option 3 — Manual g++ compilation (fallback)

First, build htslib:

```bash
cd htslib && autoreconf -i && ./configure && make -j && cd ..
```

Then compile manually:

**Linux:**

```bash
mkdir -p bin && cd bin/
g++ -O3 -fPIC -std=c++17 \
    ../src/*.cpp ../src/io/*.cpp ../htslib/libhts.a \
    -I ../htslib -I ../src \
    -lz -lbz2 -lm -llzma -lpthread -lcurl -lssl -lcrypto \
    -o mitoquest
```

**macOS:**

```bash
mkdir -p bin && cd bin/
g++ -O3 -fPIC -std=c++17 -Wl,-no_compact_unwind \
    ../src/*.cpp ../src/io/*.cpp ../htslib/libhts.a \
    -I ../htslib -I ../src \
    -lz -lbz2 -lm -llzma -lpthread -lcurl \
    -o mitoquest
```

> **Note:** If you encounter a `test/test_khash.c` compilation error during
> `make` in htslib, you can safely ignore it — the required `libhts.a`
> archive is still produced correctly.

---

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
  variant-qc   Bayesian quality control for mtDNA variants from VCF files.
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
sorted/indexed BAM or CRAM file. The autosomal chromosomes serve as the
diploid baseline (CN = 2); the mitochondrial chromosome is reported on the
same scale, i.e. the expected number of mtDNA molecules per diploid cell.
The output is a TSV with fragment counts, GC content, length-normalized
fragment ratio, and the copy-number mean + 95% confidence interval for
every contig in the BAM header.

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
                         'chr<TAB>start<TAB>end' triples; '#' starts a
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
Multi-allelic sites produce one row per ALT.  When `--gm-fam` is
specified the output switches to a wide 22-column format that includes
grandmother genotype fields and a `HAS_G` flag.

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

The `GRANDMOTHER_*` columns and `HAS_G` are present only when `--gm-fam`
is used; legacy 16-column TSVs (without `--gm-fam`) are still read
correctly by `ne-estimate`.
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

Any `Number=` declaration on `FORMAT/AD` is accepted (`R`, `A`, `.`)
because AD is always decoded via GT.  Standard `Number=R` AD
(REF + all ALTs, fixed width per record) is also handled correctly
as a special case where every sample's GT lists every allele.

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
using a Beta-Binomial **Maximum Marginal Likelihood Estimator (MMLE)**.
The latent maternal / child true allele frequencies are analytically
integrated out (Beta-Binomial conjugacy), and pairs are treated as
independent (composite / pseudo-likelihood), so the global objective is
the sum of per-pair *marginal* log-likelihoods that depends only on Ne
and the observed counts.  We label the resulting estimator **MMLE** in
the code, output, and documentation to keep the assumptions explicit
(rather than the looser "MLE" used in v1.8.5 and earlier).

### Statistical model

For each transmission pair the model assumes:

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

The global log-likelihood is the sum across independent transmission
pairs. The optimum is found by a **brute-force integer scan** over
`[--min-ne, --max-ne]` (the discrete LL is *not* unimodal in `Ne` —
fitted `Ne = k * true_Ne` aligns with the observed VAF grid and
produces secondary local maxima, so golden-section search is unsafe);
the 95% confidence interval is the contiguous range of `Ne` where
`logL >= logL_max - 1.92` (profile likelihood, `chi2_{1, 0.95} / 2`).

### Three-generation G-M-C trios (`--gm-fam`)

When `mitoquest trans-prep` is run with `--gm-fam`, the output TSV
carries a `HAS_G` column: `HAS_G = 1` rows are grandmother–mother–
child (G-M-C) trios and `HAS_G = 0` rows are standard mother–child
(MC) pairs. `mitoquest ne-estimate` auto-detects the `HAS_G`,
`GRANDMOTHER_DP`, and `GRANDMOTHER_AD_ALT` columns by name and
switches each trio row to a three-generation marginal likelihood.

For a G-M-C trio the grandmother’s observed VAF
`p̂_G = g_ad_alt / g_dp` serves as the founder allele frequency.
The mother’s latent heteroplasmy `p_M` is drawn from the Kimura
diffusion `Beta(p̂_G·(Ne−1), (1−p̂_G)·(Ne−1))`. Given `p_M`, the
mother’s read count `k_M ~ Bin(d_M, p_M)` and the child’s read
count `k_C ~ BetaBin(d_C, p_M·(Ne−1), (1−p_M)·(Ne−1))` are
conditionally independent. The mother’s latent `p_M` is
analytically marginalised, giving the closed-form trio marginal
likelihood:

```txt
I_trio(Ne) = C(d_M, k_M) · C(d_C, k_C)
           · B(α_G + k_M + k_C, β_G + (d_M − k_M) + (d_C − k_C))
           / B(α_G, β_G)

where  α_G = p̂_G · (Ne − 1),   β_G = (1 − p̂_G) · (Ne − 1)
```

The composite log-likelihood sums the per-row contributions:
trio rows use `log I_trio`, pair rows use the standard 2-gen
Beta-Binomial log-likelihood. Homoplasmic grandmothers
(`p̂_G ∈ {0, 1}`) contribute 0 (Ne-independent constant). The
mixed-pedigree MMLE is consistent for Ne under the standard
Kimura assumptions and reduces exactly to the 2-gen MMLE when
no trio rows are present.

Example:

```bash
# MC FAM: mother-child pairs
mitoquest trans-prep \
    -v cohort.vcf.gz \
    -f mc_pairs.fam \
    -g gm_pairs.fam \        # grandmother-mother FAM
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

Using raw `(DP, AD)` and marginalising the maternal posterior is the
**exact discrete-Wright-Fisher likelihood** for a single transmission;
any method that uses `AF` as if it were the truth is an approximation
to this.

### Comparison with the deCODE 2024 *Cell* paper

[Helgason et al. (2024), *Cell*](https://doi.org/10.1016/j.cell.2024.05.022)
estimated the human mtDNA bottleneck on 3,457 mother-child pairs from
Iceland. Their headline number for direct M-C transmissions is
**Ne ≈ 2.29 (95% CI 1.95–2.62)**; pooling deeper-pedigree relatives
places it at **Ne ≈ 3** (best fit in the simulated range 2–30).

Methodologically, deCODE uses the Wonnapinij bottleneck parameter `b`
(based on the Kimura diffusion; Wonnapinij 2008/2010) fitted to the
*variance* of frequency changes between relatives.

Since v1.8.2, `mitoquest ne-estimate` uses a **continuous Beta-diffusion
MMLE** as the default model.  This models the child's true heteroplasmy
as a Kimura-diffusion draw:

```txt
p_child | p_mother  ~  Beta(p_m × (Ne − 1), (1 − p_m) × (Ne − 1))
c_alt   | p_child   ~  BetaBinomial(c_dp, p_m × (Ne − 1), (1 − p_m) × (Ne − 1))
```

The continuous MMLE and the Wonnapinij/Kimura cross-check are
**theoretically consistent**: both estimate the same Ne from the same
Kimura diffusion (Var = p(1−p)/Ne), differing only in statistical
method — full marginal likelihood (MMLE) vs method-of-moments (Kimura).  On
well-behaved mtDNA data they agree closely (both yielding Ne ≈ 1–10
in human cohorts, consistent with deCODE).

For multi-generation pedigrees (e.g. cousins, deeper relatives)
Wonnapinij/Kimura is the better framework because it natively handles
`g > 1` generations of drift.

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
| **Continuous MMLE** (default) | Primary | Real-valued Ne, Cramér-Rao efficient, tighter CI, handles per-read sampling uncertainty, robust to moderate outliers | Assumes single-generation bottleneck (`g = 1`) |
| **Kimura cross-check** (`--cross-check kimura`) | Secondary validator | Method-of-moments (fast), independent confirmation, natively supports multi-generational pedigrees (`g > 1`) | Wider CI, sensitive to outliers without `--kimura-trim` |
| **Discrete MMLE** (`--model discrete`) | Specialised | Exact for virus-passage / fixed-inoculum experiments | Systematic upward bias on mtDNA at high DP; integer-only grid |

**Decision rule:**

1. Use the **continuous MMLE** as the reported Ne (default behaviour).
2. Run `--cross-check kimura --kimura-trim 0.10` as a sanity check.
3. When the two agree (CI overlap), confidence is high.
4. When they disagree (MMLE >> Kimura), inspect outlier pairs with
   `--top-drift-k 20` — a handful of high-drift pairs (NUMTs /
   sequencing errors) can collapse the variance-based Kimura estimate
   while the marginal-likelihood-based MMLE is robust.
5. For multi-generation pedigree data (`g > 1`), prefer the Kimura
   framework which supports generation-adjusted Ne.

#### Why Ne_MMLE is often smaller than Ne_Kimura

On real mtDNA data it is common to observe `Ne_MMLE < Ne_Kimura` (e.g.,
Ne_MMLE ≈ 2.7 vs Ne_Kimura ≈ 3.7).  This is **expected behaviour**, not a
bug, and the MMLE estimate is the more accurate of the two.  The gap arises
from several statistical effects:

1. **Jensen's inequality bias in the Kimura ratio estimator.**
   The Kimura formula computes `Ne = 1/V` where `V = Σ(d_i − s_i) / Σw_i`.
   Because `1/x` is a convex function, the estimator `1/V̂` is biased
   *upward* by Jensen's inequality: `E[1/V̂] > 1/E[V̂]`.  With large
   cohorts (hundreds of pairs) the CIs are tight enough to expose this
   systematic upward bias in Ne_Kimura.

2. **Full distribution vs. second moment only.**
   The MMLE fits the entire BetaBinomial shape (all moments), while the
   Kimura estimator uses only the second central moment (variance).  At
   small Ne the Beta distribution is highly non-Gaussian (often U-shaped),
   so the higher-order information captured by the MMLE carries real signal
   that the variance-only Kimura discards.

3. **Information-optimal weighting.**
   The MMLE weights each pair by its Fisher information (how much that
   pair's specific depths and counts reveal about Ne).  The Kimura method
   weights pairs only by `p_m(1 − p_m)`, ignoring how informative the
   child observation actually is.

4. **Plug-in sampling correction noise.**
   The Kimura correction term `s_i = p̂_m(1−p̂_m)/m_dp + p̂_c(1−p̂_c)/c_dp`
   uses noisy point estimates; the MMLE avoids this by marginalising the
   read-sampling process analytically via the BetaBinomial.

**Is it always MMLE < Kimura?**  No — when high-drift outliers (NUMTs,
sequencing errors) dominate, they inflate V and collapse Ne_Kimura *below*
Ne_MMLE.  The direction depends on the data:

| Scenario | Typical relationship |
| :----------: | :---------------------: |
| Clean cohort, large N, small Ne | Ne_MMLE ≲ Ne_Kimura (Jensen's bias dominates) |
| High-drift outliers present | Ne_MMLE > Ne_Kimura (outliers collapse Kimura) |
| Perfect synthetic data | Ne_MMLE ≈ Ne_Kimura within ~10–20% |

**Bottom line:** Trust the continuous MMLE as the primary Ne estimate; treat
Ne_Kimura as a qualitative confirmation that Ne is in the same order of
magnitude.  Non-overlapping CIs between the two estimators (with the MMLE
lower) is a well-understood property of method-of-moments vs.
marginal-likelihood estimators, not a sign of data problems.

### Full parameter reference of `ne-estimate`

```bash
Usage: mitoquest ne-estimate [options] -i <pairs.tsv>

Required options:
  -i, --input     FILE   Input transmission pairs TSV produced by
                         `mitoquest trans-prep`.

Optional options:
  -o, --output    FILE   JSON output file (default: stdout).
      --model     NAME   Likelihood model: `continuous` (default,
                         recommended for mtDNA) or `discrete`.
      --min-vaf   FLOAT  Lower maternal VAF gate, inclusive [0.10].
      --max-vaf   FLOAT  Upper maternal VAF gate, inclusive [0.90].
      --min-ne    INT    Smallest Ne value to consider [1].
      --max-ne    INT    Largest Ne value to consider  [200].
  -t, --threads   INT    Worker threads for the inner sum [1].
      --cross-check NAME   Optional secondary estimator alongside the
                           MMLE. Supported value: `kimura`, which computes
                           the Wonnapinij b and the implied single-
                           generation Ne (Helgason 2024).
      --kimura-bootstrap INT  Non-parametric bootstrap iterations for the
                              Kimura cross-check 95% CI [1000].
      --kimura-seed      INT  RNG seed for the Kimura bootstrap [42].
      --kimura-trim   FLOAT   Fraction of highest-drift pairs to drop
                              before recomputing b (0.0 disables,
                              recommended 0.10) [0.0].
      --top-drift-k     INT   Emit the top-K highest-drift pairs in JSON
                              for outlier inspection [0].
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
  "Estimator":       "MMLE (composite marginal likelihood)",
  "Max_Marginal_LogLik": -812.34561,
  "Model":           "continuous",
  "Min_VAF":         0.10,
  "Max_VAF":         0.90,
  "Search_Min_Ne":   1,
  "Search_Max_Ne":   200,
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
- A leading provenance comment (`#mitoquest_ne_estimate_command=...`) is
  written before the JSON when the output file is written explicitly.

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

### Recommendations

- **Sequencing depth ≥ 500×** on chrM for both mother and child; lower
  depth widens the per-pair likelihood envelope and inflates the CI.
- **Restrict to maternal VAF in [0.10, 0.90]** (the default). Sites very
  close to 0 or 1 carry virtually no information about Ne, while still
  contributing numerical noise.
- **At least ~30 informative pairs** are needed for a stable CI; with
  fewer pairs the CI tends to be wide and to clip the search boundaries.
- **Verify both `_Clipped` flags are `false`**; if either is `true`,
  widen `--min-ne` / `--max-ne` and re-run.

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
