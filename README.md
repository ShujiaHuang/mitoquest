# MitoQuest: Human Mitochondrial sequencing data Analysis Toolkit

A cross-platform, efficient and practical bioinformatic analysis toolkit written in C++, which is
designed to call mitochondrial variants (SNPs and Indels), heteroplasmy/homoplasmy, MT copy number,
etc., from whole-genome sequencing(WGS) data.

Scripts to analyze mitochondrial sequencing data.

Seeking information like heteroplasmy, SNV, etc. on human  Mitochondrial genome from NGS data.

## Build and run testing

```bash
$ rm -rf build
$ mkdir build && cd build
$ cmake -DBUILD_TESTING=ON ..
$ make
$ ctest --output-on-failure
```

> *"What I cannot create, I do not understand" — Richard Feynman*
