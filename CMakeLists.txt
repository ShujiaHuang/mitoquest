# Set the minimum version of CMake that can be used
cmake_minimum_required (VERSION 3.12.0)

# Set the C++ standard required for the project
set(CMAKE_CXX_STANDARD 17)  # C++ 17
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the project name
project (MITOQUEST VERSION 1.5.5 LANGUAGES CXX)

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
# -Wl,-no_compact_unwind: disable compact unwind information for MacOS only
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fPIC")
if (APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-no_compact_unwind")
endif()

# Ensure htslib is built before the main executable
add_custom_target(
    Buildhts ALL
    COMMAND autoreconf -i && ./configure && make
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/htslib
)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/htslib)  # include the header of htslib
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)     # include the header of mitoquest

# Define common system libraries
set(COMMON_SYS_LIB pthread z bz2 m lzma curl)
if(NOT APPLE)
    list(APPEND COMMON_SYS_LIB ssl crypto)
endif()

# Define htslib library
set(HTSLIB ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a)

# Compile mitoquest - main executable
file(GLOB_RECURSE SRC_MT 
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/io/*.cpp"
)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
add_executable(mitoquest ${SRC_MT})
add_dependencies(mitoquest Buildhts)
target_link_libraries(mitoquest ${HTSLIB} ${COMMON_SYS_LIB})

####################################################################################################
####### Compile the test and other executables scripts
####################################################################################################
# Compile count_align_fragments.cpp
# add_executable(count_align_fragments ${CMAKE_CURRENT_SOURCE_DIR}/scripts/count_align_fragments.cpp)
# add_dependencies(count_align_fragments Buildhts)
# target_link_libraries(count_align_fragments ${HTSLIB} ${COMMON_SYS_LIB})

add_executable(mtcopynum ${CMAKE_CURRENT_SOURCE_DIR}/scripts/mt_copynum.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/algorithm.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/mt_utils.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/io/iobgzf.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/io/fasta.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/io/utils.cpp
)
add_dependencies(mtcopynum Buildhts)
target_link_libraries(mtcopynum ${HTSLIB} ${COMMON_SYS_LIB})

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