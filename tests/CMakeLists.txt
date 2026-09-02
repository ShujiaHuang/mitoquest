# Copy test data files
file(COPY ${CMAKE_SOURCE_DIR}/tests/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Create test library
add_library(mitoquest_lib STATIC
    ${CMAKE_SOURCE_DIR}/src/algorithm.cpp
    ${CMAKE_SOURCE_DIR}/src/basetype.cpp
    ${CMAKE_SOURCE_DIR}/src/log_factorial.cpp
    ${CMAKE_SOURCE_DIR}/src/mt_copynum.cpp
    ${CMAKE_SOURCE_DIR}/src/mt_variant_caller.cpp
    ${CMAKE_SOURCE_DIR}/src/trans_prep.cpp
    ${CMAKE_SOURCE_DIR}/src/ne_estimate.cpp
    ${CMAKE_SOURCE_DIR}/src/mt_utils.cpp
    ${CMAKE_SOURCE_DIR}/src/variant_qc.cpp
    ${CMAKE_SOURCE_DIR}/src/vcf_subset_samples.cpp

    ${CMAKE_SOURCE_DIR}/src/io/utils.cpp
    ${CMAKE_SOURCE_DIR}/src/io/bam.cpp
    ${CMAKE_SOURCE_DIR}/src/io/bam_header.cpp
    ${CMAKE_SOURCE_DIR}/src/io/bam_record.cpp
    ${CMAKE_SOURCE_DIR}/src/io/fasta.cpp
    ${CMAKE_SOURCE_DIR}/src/io/iobgzf.cpp
    ${CMAKE_SOURCE_DIR}/src/io/vcf.cpp
    ${CMAKE_SOURCE_DIR}/src/io/vcf_header.cpp
    ${CMAKE_SOURCE_DIR}/src/io/vcf_record.cpp
)

# Make sure htslib (headers + libhts.a) is built before this library compiles.
# Without this, parallel builds can start compiling mitoquest_lib's TUs before
# htslib's ./configure has generated config.h / htslib/version.h.
add_dependencies(mitoquest_lib Buildhts)

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
    test_basetype.cpp
    test_bam_aligned_pairs.cpp
    test_fasta_safety.cpp
    test_mt_variant_caller.cpp
    mt_copynum_test.cpp
    test_trans_prep.cpp
    test_vcf_subset_samples.cpp
    test_ne_estimate.cpp
    test_mt_utils.cpp
    test_variant_qc.cpp

    io/bam_test.cpp
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

# Tiny command-line adapter used by the optional high-precision Python
# cross-check. It links the same library code as the GTest suite.
add_executable(mitoquest_trio_probe trio_likelihood_probe.cpp)
target_include_directories(mitoquest_trio_probe PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(mitoquest_trio_probe PRIVATE mitoquest_lib)
add_dependencies(mitoquest_trio_probe mitoquest)

# Register tests
add_test(NAME AllTests COMMAND mitoquest_tests)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    add_test(
        NAME TrioLikelihoodMpmathCrossCheck
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/test_trio_likelihood_mpmath.py
                $<TARGET_FILE:mitoquest_trio_probe>
            $<TARGET_FILE:mitoquest>
    )
    set_tests_properties(TrioLikelihoodMpmathCrossCheck PROPERTIES
        SKIP_RETURN_CODE 77
    )
    add_test(
        NAME NeEstimateJsonSerialization
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/test_ne_estimate_json.py
                $<TARGET_FILE:mitoquest>
    )
endif()