// Header file for the VCF Subsetting tool
// Author: Shujia Huang
// Date: 2025-04-27
#ifndef __INCLUDE_VCF_SUBSET_SAMPLES_H__
#define __INCLUDE_VCF_SUBSET_SAMPLES_H__

#include <getopt.h> // For getopt_long
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <cmath>     // For isnan, isinf
#include <limits>    // For numeric_limits
#include <stdexcept> // For std::runtime_error

#include "version.h"   // For version info in usage
#include "algorithm.h" // For ngslib::mean, ngslib::median
#include "mt_utils.h"
#include "io/utils.h"  // For get_firstcolumn_from_file, suffix_name
#include "io/vcf.h"

/**
 * @brief Class to handle the logic for the 'subsam' command.
 *
 * Parses arguments, reads an input VCF, subsets it based on provided samples,
 * updates relevant INFO fields, and writes the output VCF.
 */
class VCFSubsetSamples {
private:
    // Command line arguments
    std::string _input_vcf_path;
    std::string _output_vcf_path;
    std::vector<std::string> _samples_to_keep;
    std::string _output_mode;    // e.g., "w", "wb", "wz"
    bool _update_info   = true;  // Flag to update INFO fields
    bool _keep_all_site = false; // Flag to update INFO fields
    std::string _cmdline_string; // Command line string for VCF header

    /**
     * @brief Parses command line arguments specific to the 'subsam' command.
     * @param argc Argument count (relative to the command, e.g., argc-1 from main).
     * @param argv Argument vector (relative to the command, e.g., argv+1 from main).
     */
    void parse_args(int argc, char* argv[]);

    /**
     * @brief Prints the usage information for the 'subsam' command.
     */
    static void print_usage();

    /**
     * @brief Recalculates INFO fields like AC, AN, AF based on the subset samples.
     * @param hdr The VCFHeader object (needed for tag definitions).
     * @param rec The VCFRecord object to update (must have FORMAT fields unpacked).
     */
    //  * @param sample_indices Indices of the samples being kept in the original header.
    bool recalculate_info(const ngslib::VCFHeader& hdr, ngslib::VCFRecord& rec);

public:
    /**
     * @brief Constructor that takes command line arguments.
     * @param argc Argument count (relative to the command).
     * @param argv Argument vector (relative to the command).
     * @throws std::runtime_error if arguments are invalid.
     */
    VCFSubsetSamples(int argc, char* argv[]);

    /**
     * @brief Destructor.
     */
    ~VCFSubsetSamples() = default;

    /**
     * @brief Runs the VCF subsetting process.
     * @throws std::runtime_error on file I/O errors or processing errors.
     */
    void run();
};

#endif // __INCLUDE_VCF_SUBSET_SAMPLES_H__