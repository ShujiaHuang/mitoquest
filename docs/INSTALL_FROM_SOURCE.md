# Building MitoQuest from Source

This document describes how to compile `mitoquest` from source code. Most users should use the [pre-built binary](https://github.com/ShujiaHuang/mitoquest/releases) instead.

*Requires: C++17 compiler (GCC 7+ or Apple Clang 10+), CMake ≥ 3.12, and system libraries: zlib, bzip2, xz-utils, libcurl.*

---

## Option 1 — Build with CMake (standard dynamic build)

### Step 1 — Clone the repository (including htslib submodule)

```bash
git clone --recursive https://github.com/ShujiaHuang/mitoquest.git
cd mitoquest
```

> If you forgot `--recursive`, run: `git submodule update --init --recursive`

### Step 2 — Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

The executable `bin/mitoquest` will be produced. Verify with:

```bash
./bin/mitoquest --help
./bin/mitoquest copynum --help
```

### Step 3 (Optional) — Run the unit tests

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel
cd build && ctest --output-on-failure
```

> Requires a system GoogleTest installation (e.g., `brew install googletest`
> on macOS or `sudo apt-get install libgtest-dev` on Ubuntu 22.04+).

### Step 4 (Optional) — Build a static binary locally

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

## Option 2 — Manual g++ compilation


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

> [!NOTE]
> **Note:** If you encounter a `test/test_khash.c` compilation error during `make` in htslib, you can safely ignore it — the required `libhts.a` archive is still produced correctly.
