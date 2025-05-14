# MitoQuest: Human Mitochondrial sequencing data Analysis Toolkit

A cross-platform, efficient and practical human mitochondrial sequencing data analysis toolkit in C++.

Scripts to analyze mitochondrial sequencing data.

Seeking information like heteroplasmy, SNV, etc. on human  Mitochondrial genome from NGS data.

## Build and run testing

```
$ rm -rf build
$ mkdir build && cd build
$ cmake -DBUILD_TESTING=ON ..
$ make
$ ctest --output-on-failure
```

> *"What I cannot create, I do not understand" — Richard Feynman*
