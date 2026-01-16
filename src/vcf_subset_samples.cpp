// Implementation file for the VCF Subsetting tool
// Author: Shujia Huang
// Date: 2025-04-27
#include "vcf_subset_samples.h"

// Function to print usage information for the 'subsam' command
void VCFSubsetSamples::print_usage() {
    std::cerr << "Usage: mitoquest subsam [options] -i <input.vcf> -o <output.vcf> [-s <samplelist>] [<sample1> <sample2> ...]\n";
    std::cerr << "\nOptions:\n";
    std::cerr << "  -i, --input FILE    Input VCF/BCF file (required).\n";
    std::cerr << "  -o, --output FILE   Output VCF/BCF file (required).\n";
    std::cerr << "  -s, --sample FILE   List of sample name to keep (one per line).\n";
    std::cerr << "  -O, --output-type   TYPE Output file type [v|z|b|u] (default: Guess format based on output file extension).\n";
    std::cerr << "                          v: VCF, z: compressed VCF (bgzip), b: BCF, u: uncompressed BCF.\n";
    std::cerr << "  --no-update-info    Do not update INFO fields based on extracted samples.\n";
    std::cerr << "  --keep-all-site     Do not remove the POS which only have reference allele in extracted samples.\n";
    std::cerr << "  -h, --help          Print this help message.\n\n";
    std::cerr << "Version: " << MITOQUEST_VERSION << "\n";
    std::cerr << std::endl;
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
        // std::cout << "[INFO] Guessed output mode based on input extension: " << _output_mode << "\n";
    }
}

// Constructor
VCFSubsetSamples::VCFSubsetSamples(int argc, char* argv[]) {
    parse_args(argc, argv);
}

// Recalculate INFO fields
// Updated based on logic from MtVariantCaller::_joint_variant_in_pos
bool VCFSubsetSamples::recalculate_info(const ngslib::VCFHeader& hdr, ngslib::VCFRecord& rec) {
    // Requires GT field in FORMAT, and record unpacked for FORMAT and INFO
    if (rec.unpack(BCF_UN_FMT | BCF_UN_INFO) < 0) {
            std::cerr << "Warning: Failed to unpack record for INFO recalculation at "
                      << rec.chrom(hdr) << ":" << (rec.pos() + 1) 
                      << ". Skipping INFO update.\n";
            return true;
    }

    int n_alt = rec.n_alt(); // Number of ALT alleles (without REF), 0 if none
    if (n_alt == 0) return true; // No ALT alleles, nothing to count, keep as is.

    // 1. Get genotypes for all samples
    std::vector<std::vector<int>> genotypes;
    int max_ploidy = rec.get_genotypes(hdr, genotypes);
    if (max_ploidy <= 0) {
        // GT field missing or error reading it, cannot recalculate
        std::cerr << "Warning: GT field missing or unreadable at "
                  << rec.chrom(hdr) << ":" << (rec.pos() + 1) 
                  << ". Skipping INFO update.\n";
        return true;  // Keep original INFO
    }

    int n_samples = genotypes.size();

    // 2. Get Depth (DP) from FORMAT for all samples
    std::vector<int> fmt_dp_vec;
    int n_dp_per_sample = rec.get_format_int(hdr, "DP", fmt_dp_vec);
    if (n_dp_per_sample < 0 || n_dp_per_sample > 1) {
        throw std::runtime_error("[Error]: Error reading DP FORMAT field at "
                                 + rec.chrom(hdr) + ":" + std::to_string(rec.pos() + 1));
    }

    // 3. Get AF from FORMAT for all samples
    std::vector<float> fmt_af_vec;
    int n_af_per_sample = rec.get_format_float(hdr, "AF", fmt_af_vec);
    if (fmt_af_vec.size() != n_samples * n_af_per_sample) {
        throw std::runtime_error(
            "[Error]: Mismatch in AF FORMAT field size at "
            + rec.chrom(hdr) + ":" + std::to_string(rec.pos() + 1)
        );
    }
    if (n_af_per_sample < 0) {
        throw std::runtime_error(
            "[Error]: Error reading AF FORMAT field at "
            + rec.chrom(hdr) + ":" + std::to_string(rec.pos() + 1)
        );

    } else if (n_af_per_sample == 0) {
        // AF field not present, we cannot recalculate allele frequencies
        // Proceed without updating INFO
        std::cerr << "Warning: AF FORMAT field not present at "
                  << rec.chrom(hdr) << ":" << (rec.pos() + 1) 
                  << ". Skipping INFO update.\n";
        return true;  // Keep original INFO

    } else if (n_af_per_sample < n_alt) { // should compare with max_ploidy? No!
        // AF field does not have enough values per sample
        throw std::runtime_error(
            "[Error]: Insufficient AF FORMAT values at " 
            + rec.chrom(hdr) + ":" + std::to_string(rec.pos() + 1) 
            + ". Expected at least " + std::to_string(n_alt) + " per sample, got " 
            + std::to_string(n_af_per_sample) + "."
        );
        return true;  // Keep original INFO
    }

    // --- Statistics Containers ---
    int ref_ind_count = 0;  // Count of individuals with REF allele
    int hom_ind_count = 0;
    int het_ind_count = 0;
    int available_ind_count = 0; // Available individuals number with non-missing GT

    // DP for all non-missing GT samples, equal to available_ind_count if DP present
    std::vector<int> all_available_dp; 
    all_available_dp.reserve(n_samples);

    // VAF collectors: [allele_index][sample_values]
    std::vector<std::vector<double>> alt_all_freqs(n_alt); // any sample with non-missing VAF
    std::vector<std::vector<double>> alt_het_freqs(n_alt); // only het samples

    // --- Loop through Samples ---
    for (size_t i = 0; i < n_samples; ++i) { // `i` is the index of sample
        // Get the genotype for this sample
        const std::vector<int>& gt = genotypes[i];
        std::vector<int> non_missing_al; // Non-missing alleles for this sample
        std::map<int, int> al_to_idx;

        // Process GT
        for (int j(0); j < gt.size(); ++j) { // Allele code (0=REF, 1=ALT1, ...), may be < 1 if error
            al_to_idx[gt[j]] = j;
            if (gt[j] >= 0) {  // Check if allele is not missing (-1)
                non_missing_al.push_back(gt[j]); // Add non-missing allele
            }
        }

        // Logic from MtVariantCaller: only count non-missing genotype
        if (non_missing_al.size() > 0) {
            available_ind_count++;
            // collect DP
            if (n_dp_per_sample > 0) {
                if (std::isnan(fmt_dp_vec[i])) continue; // DP is missing for this sample, skip
                if (fmt_dp_vec[i] != ngslib::VCFRecord::INT_MISSING &&
                    fmt_dp_vec[i] != ngslib::VCFRecord::INT_VECTOR_END) {
                    all_available_dp.push_back(fmt_dp_vec[i]); // Store DP, usually one per sample
                }
            }

            // Determine Genotype Status (Ref/Hom/Het)
            bool is_het_sample = false;
            if (non_missing_al.size() == 1) { // Haploid or one valid allele
                if (non_missing_al[0] != 0) {
                    hom_ind_count++;
                } else {
                    ref_ind_count++;
                }
            } else if (non_missing_al.size() > 1) { // Diploid or more: 0/1, 1/1, 0/2, 0/1/2, etc.
                // Check if all alleles are the same
                bool all_same = std::all_of(
                    non_missing_al.begin(), 
                    non_missing_al.end(), 
                    [&](int a){ return a == non_missing_al[0]; }  // Lambda to check equality
                );
                
                if (all_same) { // Homozygous
                    if (non_missing_al[0] != 0) hom_ind_count++;
                    else ref_ind_count++;
                } else {
                    het_ind_count++;
                    is_het_sample = true;
                }
            }

            // Collect VAFs if AF FORMAT is available (Number = A or R)
            for (int al: non_missing_al) {
                if (al == 0) continue; // Skip REF allele

                // Calculate index in flattened array: sample_idx * n_values + al_to_idx[al]
                size_t af_idx = i * n_af_per_sample + al_to_idx[al];
                float val = fmt_af_vec[af_idx];
                if (std::isnan(val)) continue;
                if (val != ngslib::VCFRecord::FLOAT_MISSING && 
                    val != ngslib::VCFRecord::FLOAT_VECTOR_END) {
                    // store the VAF for this allele
                    alt_all_freqs[al - 1].push_back(val);
                    if (is_het_sample) {
                        alt_het_freqs[al - 1].push_back(val);
                    }
                }
            }
        } // End of non-missing GT processing
    } // End of sample loop
    // --- Check Validity ---
    if ((!_keep_all_site) && all_available_dp.empty()) {
        std::cerr << "[INFO] No valid genotypes for any kept samples at " 
                  << rec.chrom(hdr) << ":" << (rec.pos() + 1) 
                  << ". Skipping this record.\n";
        return false;
    }
    if ((!_keep_all_site) && (hom_ind_count + het_ind_count == 0)) return false;  // Non variants on this site

    // --- Update INFO ---

    // Determine Plasmic Type (PT)
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
    // Update PT in the record's INFO field
    rec.update_info_string(hdr, "PT", pt.c_str());

    // Update AN, REF_N, HET_N, HOM_N in the record's INFO field
    rec.update_info_int(hdr, "AN",    &available_ind_count, 1);
    rec.update_info_int(hdr, "REF_N", &ref_ind_count, 1);
    rec.update_info_int(hdr, "HET_N", &het_ind_count, 1);
    rec.update_info_int(hdr, "HOM_N", &hom_ind_count, 1);

    // Update DP_MEAN and DP_MEDIAN in the record's INFO field
    int dp_mean = (!all_available_dp.empty()) ? static_cast<int>(mean(all_available_dp)) : 0;
    int dp_median = (!all_available_dp.empty()) ? static_cast<int>(median(all_available_dp)) : 0;
    rec.update_info_int(hdr, "DP_MEAN", &dp_mean, 1);
    rec.update_info_int(hdr, "DP_MEDIAN", &dp_median, 1);

    // Update VAF_*
    std::vector<float> v_mean(n_alt), v_median(n_alt), v_mean_het(n_alt), v_median_het(n_alt);
    for (int i = 0; i < n_alt; ++i) {
        v_mean[i] = (!alt_all_freqs[i].empty()) ? static_cast<float>(mean(alt_all_freqs[i])) : 0;
        v_median[i] = (!alt_all_freqs[i].empty()) ? static_cast<float>(median(alt_all_freqs[i])) : 0;
        v_mean_het[i] = (!alt_het_freqs[i].empty()) ? static_cast<float>(mean(alt_het_freqs[i])) : 0;
        v_median_het[i] = (!alt_het_freqs[i].empty()) ? static_cast<float>(median(alt_het_freqs[i])) : 0;
    }
    rec.update_info_float(hdr, "VAF_MEAN", v_mean.data(), n_alt);
    rec.update_info_float(hdr, "VAF_MEDIAN", v_median.data(), n_alt);
    rec.update_info_float(hdr, "VAF_MEAN_HET", v_mean_het.data(), n_alt);
    rec.update_info_float(hdr, "VAF_MEDIAN_HET", v_median_het.data(), n_alt);

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
            if (!subset_rec.cleanup_genotypes(subset_hdr)) { // 清理不再出现的 ALT 等位基因, 并对样本 Genotype 进行更新
                throw std::runtime_error("Error cleaning up genotypes in subset record at "
                    + subset_rec.chrom(subset_hdr) + ":" + std::to_string(subset_rec.pos() + 1));
            }  

            // Recalculate INFO fields (AN, VAF_MEAN, ...) based on the kept samples
            if (_update_info) {
                // Note: This will modify the subset_rec in place
                bool is_valid = recalculate_info(subset_hdr, subset_rec);
                if (!is_valid) {
                    continue; // No genotypes for any kept samples at this POS, Skip it.
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

        std::cout << "[INFO] VCF subsetting finished successfully.\n";
    } catch (const std::exception& e) {
        // Rethrow to be caught by main
        throw std::runtime_error("Error during VCF subsetting: " + std::string(e.what()));
    }

    return;
} 