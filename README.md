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
Version: 1.7.1

Usage: mitoquest <command> [options]
Commands:
  caller    Mitochondrial variants and heteroplasmy/homoplasmy caller.
  subsam    Extract mitochondrial variants for specified samples from VCF
            files and output a new VCF file.
  copynum   Estimate per-chromosome (incl. mtDNA) relative copy number
            from a BAM/CRAM file.
```

In addition to the main `mitoquest` binary, the project ships:

- `tools/` — a suite of Python helper scripts for VCF QC, annotation,
  pipeline assembly, and a built-in **VQSR**-style variant recalibrator
  (`tools/vqsr/`).
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
  caller    Call mitochondrial variants (SNVs/Indels) and quantify
            heteroplasmy / homoplasmy from BAM/CRAM files.
  subsam    Extract a subset of samples from an mtDNA VCF and recompute
            INFO fields.
  copynum   Estimate per-chromosome (incl. mtDNA) relative copy number
            from a BAM/CRAM file.
```

---

## `mitoquest caller` — Variant calling and heteroplasmy/homoplasmy detection

`mitoquest caller` reads aligned reads from one or more BAM/CRAM files and
emits a multi-sample VCF containing per-sample heteroplasmy fractions
(HF) and per-site INFO fields suitable for downstream filtering and
annotation.

### Full parameter reference

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

### Full parameter reference

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

### Usage examples

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

### Usage examples

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

## Auxiliary Python tools (`tools/`)

The `tools/` directory ships a set of Python helper scripts that complement
the C++ binaries:

| Script                          | Purpose                                                                                  |
| ------------------------------- | ---------------------------------------------------------------------------------------- |
| `tools/mito_annotate.py`        | Annotate a `mitoquest caller` VCF with population, in-silico, and disease databases.     |
| `tools/parse_annotatedVCF.py`   | Convert an annotated VCF into a tidy table for downstream statistics.                    |
| `tools/mtDNA_variant_QC.py`     | Per-site / per-sample QC metrics (HF distribution, AD filters, …).                       |
| `tools/mtDNA_vcf_to_table.py`   | Flatten a multi-sample VCF into long-format TSV (one row per sample × site).             |
| `tools/parse_vcf.py`            | Generic VCF parsing helper (used by other tools).                                        |
| `tools/filter_mergedVCF.py`     | Apply hard filters on a merged mtDNA VCF.                                                |
| `tools/merge_cr_ncr_vcf.py`     | Merge coding-region and non-coding-region VCFs into a single file.                       |
| `tools/rewrite_vcf.py`          | Rewrite a VCF (header normalisation, CHROM renaming, etc.).                              |
| `tools/shift_fasta.py`          | Produce a circularly shifted FASTA (used by the join-region pipeline).                   |
| `tools/create_join_seq.py`      | Build the joined coding-region / non-coding-region reference for re-alignment.           |
| `tools/detect_NUMT_by_mtCN.py`  | Flag potential NUMT contamination using mtCN ratios per sample.                          |
| `tools/vcf_format_validator.py` | Sanity-check a VCF for downstream compatibility.                                         |
| `tools/vqsr/`                   | A VQSR-style variant recalibrator (GMM-based) for mtDNA VCFs. See `tools/README.md`.     |

Each script supports `-h / --help`. See [tools/README.md](tools/README.md) for
the original Chinese documentation of `tools/vqsr/` (formerly
`mito_classifier.py`).

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

# 3. (Optional) GMM-based variant recalibration
python tools/vqsr/vqsr.py \
    -i cohort.annotated.vcf.gz \
    -o cohort.recal.vcf.gz

# 4. mtDNA copy-number estimation per sample
while IFS= read -r bam; do
    mitoquest copynum -r reference.fasta -q 30 -t 4 "$bam"
done < bamfile.list > cohort.mtCN.tsv

# 5. Convert to a long-format table for analysis in R/Python
python tools/mtDNA_vcf_to_table.py \
    -i cohort.recal.vcf.gz \
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
