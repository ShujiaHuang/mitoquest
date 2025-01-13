cmake_minimum_required (VERSION 3.12.0)

# Project information
project (HMTK VERSION 1.2.0 LANGUAGES CXX)

# Set the C++ standard required for the project
set(CMAKE_CXX_STANDARD 17)  # C++ 17
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fPIC -Wl,-no_compact_unwind")
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# Compile htslib
add_custom_target(
    Buildhts ALL
    COMMAND autoreconf -i && ./configure --disable-shared && make
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/htslib
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/htslib)

# Find all source codes
file(GLOB SRC_CAF ${CMAKE_CURRENT_SOURCE_DIR}/scripts/count_align_fragments.cpp)

# Create executable
add_executable(count_align_fragments ${SRC_CAF})

# Ensure htslib is built before the main executable
add_dependencies(count_align_fragments Buildhts)

# Link htslib static library
target_link_libraries(count_align_fragments ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a)

set(COMMON_SYS_LIB pthread z bz2 m lzma curl)
if (APPLE)
    target_link_libraries(count_align_fragments ${COMMON_SYS_LIB})
else ()
    target_link_libraries(count_align_fragments ${COMMON_SYS_LIB} ssl crypto)
endif()

# Compile test_thread_pool.cpp
add_executable(test_thread_pool ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_thread_pool.cpp)

# Link htslib static library to the test executable
target_link_libraries(test_thread_pool ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a)

# Link common system libraries
if (APPLE)
    target_link_libraries(test_thread_pool pthread z bz2 m lzma curl)
else ()
    target_link_libraries(test_thread_pool pthread z bz2 m lzma curl ssl crypto)
endif()
