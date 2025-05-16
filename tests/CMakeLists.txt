# Create test library
add_library(mitoquest_lib STATIC
    ${CMAKE_SOURCE_DIR}/src/algorithm.cpp
    
    ${CMAKE_SOURCE_DIR}/src/io/iobgzf.cpp
    ${CMAKE_SOURCE_DIR}/src/io/vcf.cpp
    ${CMAKE_SOURCE_DIR}/src/io/vcf_header.cpp
    ${CMAKE_SOURCE_DIR}/src/io/vcf_record.cpp
    # ... adding more source files for testing
)

# Set include directories for the library
target_include_directories(mitoquest_lib PUBLIC
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/htslib
)

# Link against htslib
target_link_libraries(mitoquest_lib PUBLIC
    ${HTSLIB}
    ${COMMON_SYS_LIB}
)

# Create test executable
add_executable(mitoquest_tests
    algorithm_test.cpp

    io/iobgzf_test.cpp
    io/vcf_test.cpp
    io/vcf_record_test.cpp
    # Add other test files here as they are created
)

# Set include directories for the test executable
target_include_directories(mitoquest_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

# Link test executable with Google Test and our library
target_link_libraries(mitoquest_tests PRIVATE
    mitoquest_lib
    GTest::GTest
    GTest::Main
)

# Register tests
add_test(NAME AllTests COMMAND mitoquest_tests)