// Implementation file for the VCF Subsetting tool
// Author: Shujia Huang
// Date: 2025-04-27
#include "vcf_subset_samples.h"

// Function to print usage information for the 'subsam' command
void VCFSubsetSamples::print_usage() {
    std::cerr << "Usage: mitoquest subsam [options] -i <input.vcf> -o <output.vcf> [-s <samplelist>] [<sample1> <sample2> ...]\n";
    std::cerr << "\n";
    std::cerr << "Options:\n";
    std::cerr << "  -i, --input FILE    Input VCF/BCF file (required).\n";
    std::cerr << "  -o, --output FILE   Output VCF/BCF file (required).\n";
    std::cerr << "  -s, --sample FILE   List of sample name to keep (one per line).\n";
    std::cerr << "  -O, --output-type   TYPE Output file type [v|z|b|u] (default: Guess format based on output file extension).\n";
    std::cerr << "                          v: VCF, z: compressed VCF (bgzip), b: BCF, u: uncompressed BCF.\n";
    std::cerr << "  --no-update-info    Do not update INFO fields based on extracted samples" << std::endl;
    std::cerr << "  --keep-all-site     Do not remove the POS which only have reference allele in extracted samples" << std::endl;
    std::cerr << "  -h, --help          Print this help message.\n";
    std::cerr << "\n";
    std::cerr << "Version: " << MITOQUEST_VERSION << "\n";
    std::cerr << "\n";
}

// Parses command line arguments
void VCFSubsetSamples::parse_args(int argc, char* argv[]) {
    _output_mode = ""; // Default: determine later

    // Define long options
    static const struct option long_options[] = {
        {"input",       required_argument, 0, 'i'},
        {"output",      required_argument, 0, 'o'},
        {"sample-list", required_argument, 0, 's'},
        {"output-type", required_argument, 0, 'O'},
        {"no-update-info",    no_argument, 0, '1'},
        {"keep-all-site",     no_argument, 0, '2'},
        {"help",              no_argument, 0, 'h'},
        {0, 0, 0, 0} // Terminator
    };

    // Save the complete command line options in VCF header
    _cmdline_string = "##mitoquest_subsam_command=";
    for (size_t i = 0; i < argc; ++i) {
        _cmdline_string += (i > 0) ? " " + std::string(argv[i]) : std::string(argv[i]);
    }

    int option_index = 0;
    int c;
    std::vector<std::string> sv;
    while ((c = getopt_long(argc, argv, "i:o:s:O:h", long_options, &option_index)) != -1) {
        switch (c) {
            case 'i': _input_vcf_path = optarg;  break;
            case 'o': _output_vcf_path = optarg; break;
            case 's': 
                sv = ngslib::get_firstcolumn_from_file(optarg);
                _samples_to_keep.insert(_samples_to_keep.end(), sv.begin(), sv.end());
                break;
            case 'O': _output_mode   = optarg; break;
            case '1': _update_info   = false;  break;
            case '2': _keep_all_site = true;  break;
            case 'h': 
                print_usage();
                exit(EXIT_SUCCESS);
            case '?': 
                // getopt_long already printed an error message.
                print_usage();
                exit(EXIT_FAILURE);
            default: 
                abort(); // Should not happen
        }
    }

    // Collect BAM/CRAM files
    while (optind < argc) {
        _samples_to_keep.push_back(argv[optind++]);
    }

    // Check required arguments
    if (_input_vcf_path.empty()) {
        std::cerr << "Error: Input VCF file (-i) is required.\n";
        print_usage();
        exit(EXIT_FAILURE);
    }
    if (_output_vcf_path.empty()) {
        std::cerr << "Error: Output VCF file (-o) is required.\n";
        print_usage();
        exit(EXIT_FAILURE);
    }
    if (_samples_to_keep.empty()) {
        std::cerr << "Error: At least one sample name must be specified.\n";
        print_usage();
        exit(EXIT_FAILURE);
    }

    // Validate output mode if provided
    if (!_output_mode.empty()) {
        if (_output_mode == "v" || _output_mode == "z" || _output_mode == "b" || _output_mode == "u") {
            _output_mode = "w" + _output_mode;  // Prepend 'w' for htslib mode string
        } else {
            std::cerr << "Error: Invalid output type specified (-O). Use v, z, b, or u.\n";
            print_usage();
            exit(EXIT_FAILURE);
        }

    } else { // Determine output mode if not specified
        // Guess output format based on output file extension using helper function
        std::string lower_fname = _output_vcf_path;
        std::transform(lower_fname.begin(), lower_fname.end(), lower_fname.begin(), ::tolower); // Convert to lowercase

        if (ngslib::suffix_name(lower_fname) == ".bcf") { // Use helper from anonymous namespace
            _output_mode = "wb";
        } else if (ngslib::suffix_name(lower_fname) == ".gz") { // Use helper from anonymous namespace
            _output_mode = "wz";
        } else {
            _output_mode = "w"; // Default to uncompressed VCF for .vcf or unknown
        }
        std::cout << "[INFO] Guessed output mode based on input extension: " << _output_mode << "\n";
    }
}

// Constructor
VCFSubsetSamples::VCFSubsetSamples(int argc, char* argv[]) {
    parse_args(argc, argv);
}

// Recalculate INFO fields
bool VCFSubsetSamples::recalculate_info(const ngslib::VCFHeader& hdr, ngslib::VCFRecord& rec) 
{
    // Requires GT field in FORMAT, and record unpacked for FORMAT and INFO
    if (rec.unpack(BCF_UN_FMT | BCF_UN_INFO) < 0) {
            std::cerr << "Warning: Failed to unpack record for INFO recalculation at "
                      << rec.chrom(hdr) << ":" << (rec.pos() + 1) << ". Skipping INFO update.\n";
            return true;
    }

    int n_alt = rec.n_alt();
    if (n_alt == 0) return true; // No ALT alleles, nothing to count

    // Get genotypes for all samples
    std::vector<std::vector<int>> genotypes;
    int max_ploidy = rec.get_genotypes(hdr, genotypes);
    if (max_ploidy <= 0) {
        // GT field missing or error reading it, cannot recalculate
        std::cerr << "Warning: GT field missing or unreadable at "
                  << rec.chrom(hdr) << ":" << (rec.pos() + 1) 
                  << ". Skipping INFO update.\n";
        return true;
    }

    std::vector<int> ac(n_alt, 0); // Allele count for each ALT allele
    int an = 0;                    // Allele number (total non-missing alleles)
    int ref_ind_count = 0;         // Count of individuals with REF allele
    int hom_ind_count = 0;
    int het_ind_count = 0;
    int available_ind_count = 0;

    for (size_t i = 0; i < genotypes.size(); ++i) {
        // Get the genotype for this sample
        const std::vector<int>& gt = genotypes[i];
        std::vector<int> non_missing_al;

        // Count alleles for this sample
        for (int allele_code : gt) { // Allele code (0=REF, 1=ALT1, ...)
            if (allele_code >= 0) {  // Check if allele is not missing (-1)
                an++;                // Count this allele towards AN
                if (allele_code > 0) { // Is it an ALT allele?
                    if (allele_code - 1 < ac.size()) { // Ensure index is within bounds
                        ac[allele_code - 1]++;         // Increment count for the corresponding ALT allele
                    } else {
                        throw std::runtime_error(
                            "[Error]: Allele code (" + std::to_string(allele_code) + ") "
                            "out of bounds for ALT alleles (" + std::to_string(n_alt) + ") "
                            "at " + rec.chrom(hdr) + ":" + std::to_string(rec.pos() + 1)
                        );
                    }
                }
                non_missing_al.push_back(allele_code); // Add non-missing allele
            }
        }

        // Count homozygous and heterozygous individuals
        if (non_missing_al.size() > 0) {
            // Check if the sample is homozygous or heterozygous
            if (non_missing_al.size() == 1) {  // non-reference
                if (non_missing_al[0] != 0) {
                    hom_ind_count++;
                } else {
                    ref_ind_count++;
                }
            } else if (non_missing_al.size() > 1) {
                het_ind_count++;
            }
            available_ind_count++;
        }
    }

    if ((!_keep_all_site) && (an == 0)) {
        std::cerr << "Warning: Missing call in all samples "<< rec.chrom(hdr) << ":" << (rec.pos() + 1) << "\n";
        return false;
    }

    if ((!_keep_all_site) && (hom_ind_count + het_ind_count == 0)) return false;  // Non variants on this site

    // Update AC, AN, HOM_N, HET_N, Total_N in the record's INFO field
    rec.update_info_int(hdr, "AC", ac.data(), ac.size());
    rec.update_info_int(hdr, "AN", &an, 1);
    rec.update_info_int(hdr, "HOM_N", &hom_ind_count, 1);
    rec.update_info_int(hdr, "HET_N", &het_ind_count, 1);
    rec.update_info_int(hdr, "Total_N", &available_ind_count, 1);

    // Update AF, HOM_PF, HET_PF, SUM_PF in the record's INFO field

    // If AN is 0 (all kept samples had missing genotypes), set AF to missing or 0
    std::vector<float> af(n_alt, 0.0f); // Or std::vector<float> af(n_alt, ngslib::VCFRecord::FLOAT_MISSING);

    // If no available individuals, set PF fields to missing or 0
    float hom_pf = 0.0f; // Or ngslib::VCFRecord::FLOAT_MISSING;
    float het_pf = 0.0f; // Or ngslib::VCFRecord::FLOAT_MISSING;
    float sum_pf = 0.0f; // Or ngslib::VCFRecord::FLOAT_MISSING;
    if (an > 0) {
        for (int i = 0; i < n_alt; ++i) {
            af[i] = static_cast<float>(ac[i]) / an;
            // Handle potential NaN/Inf just in case, though unlikely here
            if (std::isnan(af[i]) || std::isinf(af[i])) {
                af[i] = 0.0f; // Or some other placeholder like VCFRecord::FLOAT_MISSING
            }
        }
        hom_pf = static_cast<double>(hom_ind_count) / available_ind_count;
        het_pf = static_cast<double>(het_ind_count) / available_ind_count;
        sum_pf = static_cast<double>(hom_ind_count + het_ind_count) / available_ind_count;
    }
    rec.update_info_float(hdr, "AF", af.data(), af.size());
    rec.update_info_float(hdr, "HOM_PF", &hom_pf, 1);
    rec.update_info_float(hdr, "HET_PF", &het_pf, 1);
    rec.update_info_float(hdr, "SUM_PF", &sum_pf, 1);

    // Update PT in the record's INFO field
    std::string pt; // plasmic type
    if (ref_ind_count > 0 && hom_ind_count + het_ind_count == 0) {
        pt = "Ref";
    } else if (hom_ind_count > 0 && het_ind_count == 0) {
        pt = "Hom";
    } else if (het_ind_count > 0 && hom_ind_count == 0) {
        pt = "Het";
    } else if (het_ind_count > 0 && hom_ind_count > 0) {
        pt = "Mixed";
    } else {
        pt = "Unknown"; // Fallback case
    }
    rec.update_info_string(hdr, "PT", pt.c_str());

    return true;
}

// Main execution logic
void VCFSubsetSamples::run() {
    try {
        // 1. Open Input VCF
        ngslib::VCFFile reader(_input_vcf_path);
        if (!reader.is_open()) {
            throw std::runtime_error("Failed to open input VCF: " + _input_vcf_path);
        }
        ngslib::VCFHeader original_hdr = reader.header(); // Get a copy we can query

        // Verify requested samples exist in the original header
        std::set<std::string> original_sample_set;
        std::vector<std::string> original_samples_vec = original_hdr.sample_names();
        for(const auto& s : original_samples_vec) {
            original_sample_set.insert(s);
        }

        std::vector<int> sample_indices; // Store 0-based indices of kept samples in original header
        sample_indices.reserve(_samples_to_keep.size());
        for (const auto& sample_name : _samples_to_keep) {
            if (original_sample_set.find(sample_name) == original_sample_set.end()) {
                throw std::runtime_error("Sample '" + sample_name + "' not found in the input VCF header.");
            }
            // Find the original index (needed for recalculate_info)
            int idx = original_hdr.sample_index(sample_name);
            if (idx >= 0) {
                sample_indices.push_back(idx);
            } else {
                // Should not happen based on set check, but defensive coding
                throw std::runtime_error("Internal error: Could not find index for sample '" + sample_name + "'.");
            }
        }

        // 2. Create Subset Header
        ngslib::VCFHeader subset_hdr = original_hdr.subset_samples(_samples_to_keep);
        if (!subset_hdr.is_valid()) {
            throw std::runtime_error("Failed to create subset VCF header.");
        }
        // Optional: Add a header line indicating the subsetting operation
        subset_hdr.add_header_line(_cmdline_string);

        // 3. Open Output VCF
        ngslib::VCFFile outvcf(_output_vcf_path, subset_hdr, _output_mode);
        if (!outvcf.is_open()) {
            throw std::runtime_error("Failed to open output VCF: " + _output_vcf_path);
        }

        // 4. Read, Process, Write Records
        ngslib::VCFRecord rec;
        long record_count = 0;
        while (reader.read(rec) >= 0) { // Read until EOF (-1) or error (< -1)
            record_count++;

            // 在这里添加记录子集化处理
            ngslib::VCFRecord subset_rec = rec.subset_samples(subset_hdr, sample_indices);
            if (!subset_rec.cleanup_alleles(subset_hdr)) { // 先清理不再出现的 ALT 等位基因
                throw std::runtime_error("Error cleaning up alleles in subset record at "
                    + subset_rec.chrom(subset_hdr) + ":" + std::to_string(subset_rec.pos() + 1));
            }  

            // Recalculate INFO fields (AC, AN, AF, ...) based on the subset of samples
            if (_update_info) {
                // Note: This will modify the subset_rec in place
                bool is_valid = recalculate_info(subset_hdr, subset_rec);
                if (!is_valid) {
                    std::cout << "[INFO] No valid genotypes for any kept samples at "
                              << subset_rec.chrom(subset_hdr) << ":" << (subset_rec.pos() + 1)
                              << ". Skipping this record.\n";
                    continue; // Skip this record
                }
            }

            // Note: Subsetting FORMAT fields is handled implicitly by htslib's bcf_subset
            // function when creating the subset header, and bcf_write correctly writes
            // only the FORMAT data for the samples present in the output header.
            // We don't need to manually subset FORMAT fields here.

            // Write the potentially modified record
            if (outvcf.write(subset_rec) < 0) {
                throw std::runtime_error("Error writing VCF record to: " + _output_vcf_path);
            }

            if (record_count % 100000 == 0) {
                std::cout << "[INFO] Processed " << record_count << " records...\r";
            }
        }
        std::cout << "[INFO] Processed " << record_count << " records.\n"; // Final count

        if (reader.io_status() < -1) {
            throw std::runtime_error("Error reading input VCF file.");
        }

        // 5. Cleanup (done by destructors of VcfReader, VcfWriter)
        std::cout << "[INFO] VCF subsetting finished successfully.\n";

    } catch (const std::exception& e) {
        // Rethrow to be caught by main
        throw std::runtime_error("Error during VCF subsetting: " + std::string(e.what()));
    }

    return;
} 

