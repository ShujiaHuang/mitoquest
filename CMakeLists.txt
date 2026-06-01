# Set the minimum version of CMake that can be used
cmake_minimum_required (VERSION 3.12.0)

# Set the C++ standard required for the project
set(CMAKE_CXX_STANDARD 17)  # C++ 17
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the project name
project (MITOQUEST VERSION 1.9.1 LANGUAGES CXX)

# Add project metadata
set(PROJECT_AUTHOR "Shujia Huang")
set(PROJECT_EMAIL "hshujia@qq.com")
set(PROJECT_DESCRIPTION "mitoquest: Human Mitochondrial sequencing data Analysis Toolkit")

# Configure version header
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/version.h.in
    ${CMAKE_CURRENT_SOURCE_DIR}/src/version.h
)

# Set compiler flags
# -O3: optimization level 3
# -fPIC: generate position-independent code
add_compile_options(-O3 -fPIC)

# -Wl,-no_compact_unwind silences ld warnings about missing compact unwind
# information in C-only libraries (e.g. htslib) on macOS. We apply it ONLY to
# the production `mitoquest` executable below (via target_link_options) and
# NOT directory-wide, because doing so breaks gtest's EXPECT_THROW under
# Apple's linker (cross-TU C++ exceptions get std::terminate'd instead of
# being caught by the test framework).

# -------------------------------------------------------
# Static linking configuration
# -------------------------------------------------------
# Linux  : "Portable" partial-static executable.  We build on a glibc
#          host (Ubuntu in CI) and statically link:
#            * libstdc++ and libgcc   (-static-libstdc++ -static-libgcc)
#            * htslib + zlib + bzip2 + lzma + openssl  (.a archives)
#          glibc itself remains dynamic (the only sane choice on glibc
#          systems).  The result runs unmodified on any reasonably
#          modern Linux distro (Ubuntu 20.04+, Debian 11+, RHEL 8+,
#          CentOS 8+, Fedora 32+, ...).
#
# macOS  : Fully static executables are NOT supported by Apple.  We
#          statically link the third-party .a archives; only system
#          frameworks (libSystem, libc++abi, ...) remain dynamic.
# -------------------------------------------------------
if(STATIC_BUILD)
    if(APPLE)
        # Apple Clang does NOT support -static-libgcc / -static-libstdc++.
        # Strategy: pass all third-party archives explicitly as .a files and
        # let the linker prefer static over dynamic for those.  Only Apple
        # system frameworks (/usr/lib/libSystem.B.dylib, etc.) will remain
        # dynamic — this is the best achievable portability on macOS.
        #
        # We intentionally leave CMAKE_EXE_LINKER_FLAGS unchanged here;
        # the static .a paths are handled below in target_link_libraries
        # by replacing the -l<name> flags with explicit .a paths.
        message(STATUS "STATIC_BUILD is ON (macOS) — linking third-party libs statically via explicit .a paths")
    else()
        # Linux portable static build (glibc-based).  Detailed flags are set
        # in the target_link_libraries section below.
        message(STATUS "STATIC_BUILD is ON (Linux/glibc) - producing a portable basevar (libstdc++/libgcc/htslib/compression libs all bundled, glibc dynamic)")
    endif()
endif()

# Ensure htslib is built before the main executable.
# For static builds the configure flags below are used:
# --disable-libcurl  : we only read local files, no remote URLs.
# --without-libdeflate: htslib's ./configure auto-detects libdeflate when
#                       present on the host (e.g. ubuntu-22.04 GH runners,
#                       or developer machines with `brew install libdeflate`),
#                       and then emits libhts.a with libdeflate_* references.
#                       Since we do not link -ldeflate (and want to keep the
#                       binary dependency-light), we explicitly disable it.
#                       htslib falls back to zlib transparently with no
#                       functional loss.  Apply on BOTH Linux and macOS
#                       static builds for consistency.
if(STATIC_BUILD)
    set(HTSLIB_CONFIGURE_FLAGS --disable-libcurl --without-libdeflate)
else()
    set(HTSLIB_CONFIGURE_FLAGS "")
endif()

# Build htslib incrementally and in parallel.
# - add_custom_command with OUTPUT makes the build incremental (rebuild only
#   when libhts.a is missing or its dependency stamp changed).
# - The stamp file encodes HTSLIB_CONFIGURE_FLAGS so a flag change forces a
#   full reconfigure (avoids stale config.status when toggling STATIC_BUILD).
# - $(MAKE) -j${NPROC} parallelizes the build.
include(ProcessorCount)
ProcessorCount(NPROC)
if(NPROC EQUAL 0)
    set(NPROC 2)
endif()

set(HTSLIB_A     ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a)
set(HTSLIB_STAMP ${CMAKE_CURRENT_BINARY_DIR}/htslib_flags.stamp)
string(REPLACE ";" " " _HTSLIB_FLAGS_STR "${HTSLIB_CONFIGURE_FLAGS}")
file(WRITE ${HTSLIB_STAMP}.in "${_HTSLIB_FLAGS_STR}\n")
configure_file(${HTSLIB_STAMP}.in ${HTSLIB_STAMP} COPYONLY)

add_custom_command(
    OUTPUT ${HTSLIB_A}
    COMMAND ${CMAKE_COMMAND} -E echo "Configuring htslib with: ${HTSLIB_CONFIGURE_FLAGS}"
    COMMAND $(MAKE) distclean || true
    COMMAND autoreconf -i
    COMMAND ./configure ${HTSLIB_CONFIGURE_FLAGS}
    COMMAND $(MAKE) -j${NPROC} libhts.a
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/htslib
    DEPENDS ${HTSLIB_STAMP}
    COMMENT "Building htslib (parallel, j=${NPROC})"
    VERBATIM
)
add_custom_target(Buildhts DEPENDS ${HTSLIB_A})

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/htslib)  # include the header of htslib
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)     # include the header of mitoquest

# Define common system libraries
set(COMMON_SYS_LIB pthread z bz2 m lzma curl)
if(NOT APPLE)
    list(APPEND COMMON_SYS_LIB ssl crypto)
endif()

# Define htslib library
set(HTSLIB ${HTSLIB_A})

# Compile mitoquest - main executable.
# NOTE: only one GLOB_RECURSE pattern is needed; the previous version listed
# both "src/*.cpp" and "src/io/*.cpp" which caused every io/*.cpp file to be
# included twice (recursive glob already traverses subdirectories).
# CONFIGURE_DEPENDS makes CMake re-glob when files are added/removed.
file(GLOB_RECURSE SRC_MT CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
add_executable(mitoquest ${SRC_MT})
if (APPLE)
    target_link_options(mitoquest PRIVATE -Wl,-no_compact_unwind)
endif()
add_dependencies(mitoquest Buildhts)
# Link libraries
if(STATIC_BUILD AND APPLE)
    # ---- macOS best-effort static build ----
    # Apple does not support fully static executables.
    # Strategy: statically link zlib, bz2, lzma via Homebrew .a archives.
    # Homebrew's libcurl.a requires ngtcp2/libicucore (QUIC/HTTP3+IDN) which
    # are complex to resolve as .a on macOS — use system curl dylib instead.
    #
    # Homebrew on Apple Silicon: /opt/homebrew/opt/<pkg>
    # Homebrew on Intel Macs:    /usr/local/opt/<pkg>
    find_library(ZLIB_STATIC   libz.a
        HINTS /opt/homebrew/opt/zlib/lib /usr/local/opt/zlib/lib
        REQUIRED NO_DEFAULT_PATH)
    find_library(BZ2_STATIC    libbz2.a
        HINTS /opt/homebrew/opt/bzip2/lib /usr/local/opt/bzip2/lib
        REQUIRED NO_DEFAULT_PATH)
    find_library(LZMA_STATIC   liblzma.a
        HINTS /opt/homebrew/opt/xz/lib /usr/local/opt/xz/lib
        REQUIRED NO_DEFAULT_PATH)

    message(STATUS "Static libs: zlib=${ZLIB_STATIC} bz2=${BZ2_STATIC} lzma=${LZMA_STATIC}")
    message(STATUS "Dynamic (system): curl, ssl, crypto")

    target_link_libraries(mitoquest
        ${HTSLIB}
        curl            # system /usr/lib/libcurl.4.dylib
        ${LZMA_STATIC}
        ${BZ2_STATIC}
        ${ZLIB_STATIC}
        pthread
        m
    )
elseif(STATIC_BUILD AND NOT APPLE)
    # ---- Linux portable partial-static build (glibc-based) ----
    # Build host: Ubuntu (CI uses ubuntu-22.04, glibc 2.35).
    # Strategy:
    #   * libstdc++ + libgcc are bundled via -static-libstdc++ -static-libgcc
    #   * htslib + zlib + bzip2 + lzma + openssl are linked as .a archives
    #   * glibc remains DYNAMIC (the binary still depends on the host's
    #     /lib64/ld-linux-x86-64.so.2 + libc.so.6, but those exist on every
    #     glibc-based Linux distribution).
    # This produces a binary that runs unmodified on any Linux with
    # glibc >= 2.31 (Ubuntu 20.04+, Debian 11+, RHEL/CentOS 8+, Fedora 32+).
    find_library(ZLIB_STATIC   NAMES libz.a     REQUIRED)
    find_library(BZ2_STATIC    NAMES libbz2.a   REQUIRED)
    find_library(LZMA_STATIC   NAMES liblzma.a  REQUIRED)
    find_library(SSL_STATIC    NAMES libssl.a   REQUIRED)
    find_library(CRYPTO_STATIC NAMES libcrypto.a REQUIRED)

    message(STATUS "Linux static libs:")
    message(STATUS "  zlib=${ZLIB_STATIC} bz2=${BZ2_STATIC} lzma=${LZMA_STATIC}")
    message(STATUS "  ssl=${SSL_STATIC} crypto=${CRYPTO_STATIC}")

    # Bundle the C++ runtime statically for THIS target only; glibc stays
    # dynamic.  Using target_link_options instead of mutating the global
    # CMAKE_EXE_LINKER_FLAGS avoids leaking these flags into other targets
    # (e.g. mtcopynum, gtest binaries).
    target_link_options(mitoquest PRIVATE -static-libstdc++ -static-libgcc)

    target_link_libraries(mitoquest
        ${HTSLIB}
        ${SSL_STATIC}
        ${CRYPTO_STATIC}
        ${LZMA_STATIC}
        ${BZ2_STATIC}
        ${ZLIB_STATIC}
        pthread
        m
        dl
    )
else()
    # ---- Standard dynamic build ----
    target_link_libraries(mitoquest ${HTSLIB} ${COMMON_SYS_LIB})
endif()



####################################################################################################
####### Compile the test and other executables scripts
####################################################################################################
# NOTE: The standalone `mtcopynum` executable has been folded into `mitoquest` as the
# `mitoquest copynum` subcommand (src/mt_copynum.{h,cpp}). The recursive GLOB above
# picks the new sources up automatically; no extra add_executable target is needed.

####################################################################################################
####### Testing Configuration
####################################################################################################
option(BUILD_TESTING "Build the testing tree" OFF)  # 默认关闭测试构建

if(BUILD_TESTING)
    # Enable testing support
    enable_testing()

    # Define common GTest installation paths
    set(GTEST_SEARCH_PATHS
        "/usr/local"
        "/opt/homebrew"  # Homebrew default path on Apple Silicon
        "/usr/local/Cellar"  # Homebrew default path on Intel Mac
        "/opt/googletest"
        "/usr"
    )

    # Try to find GTest headers and libraries
    find_path(GTEST_INCLUDE_DIR
        NAMES gtest/gtest.h
        PATHS ${GTEST_SEARCH_PATHS}
        PATH_SUFFIXES include
        NO_DEFAULT_PATH
    )

    find_library(GTEST_LIBRARY
        NAMES gtest
        PATHS ${GTEST_SEARCH_PATHS}
        PATH_SUFFIXES lib lib64
        NO_DEFAULT_PATH
    )

    # Print debug information
    message(STATUS "Looking for GTest headers in: ${GTEST_INCLUDE_DIR}")
    message(STATUS "Looking for GTest library in: ${GTEST_LIBRARY}")

    if(NOT GTEST_INCLUDE_DIR OR NOT GTEST_LIBRARY)
        message(WARNING "GTest not found in standard locations. You may need to install it:")
        message(WARNING "  Method 1: brew install googletest")
        message(WARNING "  Method 2: Download and build from source:")
        message(WARNING "    git clone https://github.com/google/googletest.git")
        message(WARNING "    cd googletest && mkdir build && cd build")
        message(WARNING "    cmake .. && make && sudo make install")
    endif()

    # Find GTest package
    find_package(GTest REQUIRED)

    if(APPLE)
        # Set RPATH options for macOS
        set(CMAKE_MACOSX_RPATH ON)                  # Enable macOS RPATH support
        set(CMAKE_SKIP_BUILD_RPATH FALSE)           # Do not skip build RPATH
        set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)    # Use install RPATH during build
        set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE) # Use link paths for RPATH
        
        # Set RPATH paths
        list(APPEND CMAKE_INSTALL_RPATH 
            "@rpath"
            "@loader_path"
            "${GTEST_INCLUDE_DIR}"
            "${GTEST_INCLUDE_DIR}/../lib"  # Adjusted for GTest library location
        )
    endif()
    
    # Add test directory
    add_subdirectory(tests)
endif()