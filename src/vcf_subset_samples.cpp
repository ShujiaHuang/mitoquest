// Implementation file for the VCF Subsetting tool
// Author: Shujia Huang
// Date: 2025-04-27
#include "vcf_subset_samples.h"

#include <getopt.h>
#include <iostream>
#include <set>
#include <map>
#include <cmath>
#include <stdexcept>

#include "version.h"
#include "algorithm.h"
#include "io/vcf.h"
#include "io/utils.h"

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
            case '2': _keep_all_site = true;   break;
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
    if (rec.unpack(ngslib::VCFRecord::UNPACK_FMT | ngslib::VCFRecord::UNPACK_INFO) < 0) {
            std::cerr << "Warning: Failed to unpack record for INFO recalculation at "
                      << rec.chrom(hdr) << ":" << (rec.pos() + 1) 
                      << ". Skipping INFO update.\n";
            return true;
    }

    int n_alt = rec.n_alt(); // Number of ALT alleles (without REF), 0 if none
    if (n_alt == 0) return _keep_all_site;

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

    // 2. Get Depth (DP) from FORMAT for all samples. DP is optional for
    // generic VCFs; it is only needed for MitoQuest-specific summaries.
    std::vector<int32_t> fmt_dp_vec;
    int n_dp_per_sample = rec.get_format_int(hdr, "DP", fmt_dp_vec);
    const bool has_dp = n_dp_per_sample == 1 && fmt_dp_vec.size() == static_cast<size_t>(n_samples);
    if (n_dp_per_sample > 1) {
        throw std::runtime_error("[Error]: Error reading DP FORMAT field at "
                                 + rec.chrom(hdr) + ":" + std::to_string(rec.pos() + 1));
    }

    // 3. Get AF from FORMAT for all samples
    std::vector<float> fmt_af_vec;
    int n_af_per_sample = rec.get_format_float(hdr, "AF", fmt_af_vec);
    if (n_af_per_sample > 0 && fmt_af_vec.size() != static_cast<size_t>(n_samples * n_af_per_sample)) {
        throw std::runtime_error(
            "[Error]: Mismatch in AF FORMAT field size at "
            + rec.chrom(hdr) + ":" + std::to_string(rec.pos() + 1)
        );
    }
    const bool has_af = n_af_per_sample > 0;
    const ngslib::VCFHeader::FieldNumber af_number_type = hdr.format_number("AF");

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
    std::vector<int32_t> allele_counts(static_cast<size_t>(n_alt) + 1, 0);
    int32_t standard_an = 0;
    int32_t standard_ns = 0;

    // --- Loop through Samples ---
    for (size_t i = 0; i < n_samples; ++i) { // `i` is the index of sample
        // Get the genotype for this sample
        const std::vector<int>& gt = genotypes[i];
        std::vector<int> non_missing_al; // Non-missing alleles for this sample
        std::map<int, size_t> al_to_gt_index;

        // Process GT
        for (size_t j = 0; j < gt.size(); ++j) { // Allele code (0=REF, 1=ALT1, ...), may be < 1 if error
            if (gt[j] >= 0 && gt[j] <= n_alt) {
                non_missing_al.push_back(gt[j]); // Add non-missing allele
                al_to_gt_index[gt[j]] = j;
                ++allele_counts[static_cast<size_t>(gt[j])];
                ++standard_an;
            }
        }

        // Logic from MtVariantCaller: only count non-missing genotype
        if (non_missing_al.size() > 0) {
            available_ind_count++;  // has non-missing
            ++standard_ns;
            // collect DP
            if (has_dp) {
                if (fmt_dp_vec[i] != ngslib::VCFRecord::INT32_MISSING &&
                    fmt_dp_vec[i] != ngslib::VCFRecord::INT32_VECTOR_END) {
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

            // Collect each non-reference allele at most once per sample. AF
            // uses standard Number=A/R indexing when declared that way; only
            // MitoQuest's Number=. layout is aligned to GT positions.
            //
            // MitoQuest caller also writes AF GT-aligned (one value per GT
            // entry, in GT order) even when older headers declare Number=A/R
            // (e.g. tests/data/t.vcf.gz). Detect the actual layout: a
            // per-sample value count that differs from the declared
            // cardinality, or equals the sample's GT entry count, means
            // GT-aligned, and the value must then be looked up by GT
            // position (al_to_gt_index) instead of by allele index.
            const bool af_gt_aligned = [&]() {
                if (af_number_type == ngslib::VCFHeader::FieldNumber::VAR ||
                    af_number_type == ngslib::VCFHeader::FieldNumber::UNKNOWN) return true;
                const int cardinality = (af_number_type == ngslib::VCFHeader::FieldNumber::A) ? n_alt : n_alt + 1;
                if (n_af_per_sample != cardinality) return true;
                for (size_t s = 0; s < n_samples; ++s) {
                    int n_vals = n_af_per_sample;
                    for (int j = 0; j < n_af_per_sample; ++j) {
                        const float v = fmt_af_vec[s * static_cast<size_t>(n_af_per_sample) + static_cast<size_t>(j)];
                        if (ngslib::VCFRecord::is_float_vector_end(v)) { n_vals = j; break; }
                    }
                    if (n_vals > 0 && n_vals == static_cast<int>(genotypes[s].size())) return true;
                }
                return false;
            }();

            std::set<int> seen_alt_alleles;
            for (int al: non_missing_al) {
                if (al == 0 || !seen_alt_alleles.insert(al).second || !has_af) continue;

                size_t af_offset = 0;
                if (af_gt_aligned) {
                    const auto gt_position = al_to_gt_index.find(al);
                    if (gt_position == al_to_gt_index.end()) continue;
                    af_offset = gt_position->second;
                } else if (af_number_type == ngslib::VCFHeader::FieldNumber::A) {
                    af_offset = static_cast<size_t>(al - 1);
                } else {
                    af_offset = static_cast<size_t>(al);  // R layout (FieldNumber::R)
                }
                if (af_offset >= static_cast<size_t>(n_af_per_sample)) continue;
                const size_t af_idx = i * static_cast<size_t>(n_af_per_sample) + af_offset;
                float val = fmt_af_vec[af_idx];
                if (std::isfinite(val) &&
                    !ngslib::VCFRecord::is_float_missing(val) &&
                    !ngslib::VCFRecord::is_float_vector_end(val)) 
                {
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
    // DP is an optional FORMAT field (absent in generic VCFs), so the
    // "no usable sample" gate must key on non-missing GTs, not on DP.
    if ((!_keep_all_site) && available_ind_count == 0) {
        std::cerr << "[INFO] No valid genotypes for any kept samples at " 
                  << rec.chrom(hdr) << ":" << (rec.pos() + 1) 
                  << ". Skipping this record.\n";
        return false;
    }
    if ((!_keep_all_site) && (hom_ind_count + het_ind_count == 0)) return false;  // Non variants on this site

    // --- Update INFO ---
    // Recalculate existing standard VCF summaries from the subset GTs.
    // MitoQuest's AN is explicitly documented as a sample count, so retain
    // that convention only when its companion custom INFO schema is present.
    const bool mitoquest_info_schema = hdr.has_info_tag("PT") && hdr.has_info_tag("REF_N") &&
        hdr.has_info_tag("HET_N") && hdr.has_info_tag("HOM_N") && hdr.has_info_tag("VAF_MEAN");
    if (!mitoquest_info_schema && hdr.has_info_tag("AN")) {
        rec.update_info_int(hdr, "AN", &standard_an, 1);
    }
    if (hdr.has_info_tag("AC")) {
        rec.update_info_int(hdr, "AC", allele_counts.data() + 1, n_alt);
    }
    if (hdr.has_info_tag("AF")) {
        std::vector<float> standard_af(static_cast<size_t>(n_alt), 0.0f);
        if (standard_an > 0) {
            for (int allele = 0; allele < n_alt; ++allele) {
                standard_af[static_cast<size_t>(allele)] = static_cast<float>(
                    static_cast<double>(allele_counts[static_cast<size_t>(allele) + 1]) /
                    static_cast<double>(standard_an));
            }
        }
        rec.update_info_float(hdr, "AF", standard_af.data(), n_alt);
    }
    if (hdr.has_info_tag("NS")) {
        rec.update_info_int(hdr, "NS", &standard_ns, 1);
    }

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
    if (hdr.has_info_tag("PT")) {
        rec.update_info_string(hdr, "PT", pt.c_str());
    }
    if (mitoquest_info_schema) {
        const int32_t custom_an = available_ind_count;
        const int32_t custom_ns = available_ind_count;
        const int32_t custom_ref_n = ref_ind_count;
        const int32_t custom_het_n = het_ind_count;
        const int32_t custom_hom_n = hom_ind_count;

        // MitoQuest AN is a sample count; mirror it into the standard NS
        // tag so downstream tools read a conventional sample-count field.
        if (hdr.has_info_tag("AN")) {
            rec.update_info_int(hdr, "AN", &custom_an, 1);
        }
        if (hdr.has_info_tag("NS")) {
            rec.update_info_int(hdr, "NS", &custom_ns, 1);
        }
        rec.update_info_int(hdr, "REF_N", &custom_ref_n, 1);
        rec.update_info_int(hdr, "HET_N", &custom_het_n, 1);
        rec.update_info_int(hdr, "HOM_N", &custom_hom_n, 1);

        if (has_dp) {
            const int32_t dp_mean = (!all_available_dp.empty())
                ? static_cast<int32_t>(mean(all_available_dp)) : 0;
            const int32_t dp_median = (!all_available_dp.empty())
                ? static_cast<int32_t>(median(all_available_dp)) : 0;
            rec.update_info_int(hdr, "DP_MEAN", &dp_mean, 1);
            rec.update_info_int(hdr, "DP_MEDIAN", &dp_median, 1);
        }

        if (has_af) {
            std::vector<float> v_mean(n_alt), v_mean_het(n_alt);
            for (int i = 0; i < n_alt; ++i) {
                v_mean[i] = (!alt_all_freqs[i].empty())
                    ? static_cast<float>(sum(alt_all_freqs[i]) / available_ind_count) : 0;
                v_mean_het[i] = (!alt_het_freqs[i].empty())
                    ? static_cast<float>(mean(alt_het_freqs[i])) : 0;
            }
            rec.update_info_float(hdr, "VAF_MEAN", v_mean.data(), n_alt);
            rec.update_info_float(hdr, "VAF_MEAN_HET", v_mean_het.data(), n_alt);
        }
    }

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

            // ALT trimming is independent of whether INFO is recalculated.
            // Keep-all-site is the sole opt-in for retaining a record that
            // became reference-only after sample subsetting.
            if (!_keep_all_site && subset_rec.n_alt() == 0) {
                continue;
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