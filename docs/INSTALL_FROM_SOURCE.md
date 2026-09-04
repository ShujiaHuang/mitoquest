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

> [!WARNING]
> `ctest` can print **`100% tests passed`** while the high-precision reference
> check never actually ran. `TrioLikelihoodMpmathCrossCheck` executes
> `tests/test_trio_likelihood_mpmath.py` with whichever interpreter
> `find_package(Python3)` happened to discover, and that script exits with code
> 77 when `mpmath` is not importable *for that interpreter*
> (`SKIP_RETURN_CODE 77` in `tests/CMakeLists.txt`). ctest renders 77 as
> `Skipped`, which does not count as a failure. On a machine with several
> Pythons installed this silently downgrades the suite from 3 real tests to 2.
>
> To make it run, point CMake at an Python interpreter that has `mpmath`:
>
> ```bash
> cmake -B build -DPython3_EXECUTABLE=/path/to/python-with-mpmath
> cmake --build build --parallel
> cd build && ctest --output-on-failure    # test #2 should now report Passed (~5 s)
> ```
>
> Read the per-test lines, not just the summary line. A `Skipped` entry means
> "not verified", never "verified".

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
g++ -O3 -fPIC -std=c++17 \
    ../src/*.cpp ../src/io/*.cpp ../htslib/libhts.a \
    -I ../htslib -I ../src \
    -lz -lbz2 -lm -llzma -lpthread -lcurl \
    -o mitoquest
```

> [!WARNING]
> Do **not** add `-Wl,-no_compact_unwind` here. It silences the linker warnings that htslib's C-only 
> objects emit on macOS, but it also strips the unwind tables the C++ runtime needs to throw across translation units. 
> Every `mitoquest` subcommand reports user-facing errors by throwing `std::runtime_error` into the `try`/`catch` 
> in `src/main.cpp`, so a binary linked with that flag turns each of them into `libc++abi: terminating due to uncaught exception` 
> plus SIGABRT (exit 134) instead of `Error: ...` plus exit 1. The recipe above emits no linker warnings at all on 
> Apple clang (arm64 macOS); if your toolchain does emit them, they are harmless — ignore them.

> [!NOTE]
> **Note:** If you encounter a `test/test_khash.c` compilation error during `make` in htslib, you can safely ignore it — the required `libhts.a` archive is still produced correctly.
