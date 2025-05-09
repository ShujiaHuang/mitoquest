// The C++ codes for VCF/BCF Header
// Author: Shujia Huang
// Date: 2025-04-17
#include "vcf_header.h"

namespace ngslib {

    // main constructor
    VCFHeader::VCFHeader(bcf_hdr_t *hdr_ptr) : _hdr(hdr_ptr, bcf_hdr_destroy) {
        if (!hdr_ptr) {
            // If a null pointer is passed, create an empty header instead
            // to maintain a valid state.
            _hdr.reset(bcf_hdr_init("w"), bcf_hdr_destroy);
            if (!_hdr) {
                 throw std::runtime_error("Failed to initialize empty VCF header.");
            }
        }
    }

    // Default constructor
    VCFHeader::VCFHeader() : _hdr(bcf_hdr_init("w"), bcf_hdr_destroy) {
        if (!_hdr) {
            throw std::runtime_error("Failed to initialize VCF header.");
        }
    }

    // Copy constructor
    VCFHeader::VCFHeader(const VCFHeader& other) : _hdr(other._hdr) {}

    // Move constructor
    VCFHeader::VCFHeader(VCFHeader&& other) noexcept : _hdr(std::move(other._hdr)) {}

    // Copy assignment operator
    VCFHeader& VCFHeader::operator=(const VCFHeader& other) {
        if (this != &other) {
            _hdr = other._hdr; // Share ownership
        }
        return *this;
    }

    // Move assignment operator
    VCFHeader& VCFHeader::operator=(VCFHeader&& other) noexcept {
        if (this != &other) {
            _hdr = std::move(other._hdr);
        }
        return *this;
    }

    // Creates a deep copy of the underlying htslib header.
    VCFHeader VCFHeader::copy_header() const {
        if (!is_valid()) {
            throw std::runtime_error("Cannot copy an invalid VCF header.");
        }
        bcf_hdr_t* dup_hdr = bcf_hdr_dup(_hdr.get());
        if (!dup_hdr) {
            throw std::runtime_error("Failed to duplicate VCF header.");
        }
        // Return a new VCFHeader managing the duplicated raw pointer
        return VCFHeader(dup_hdr);
    }

    // Gets the number of samples in the header.
    int VCFHeader::n_samples() const {
        return is_valid() ? bcf_hdr_nsamples(_hdr.get()) : 0;
    }

    // Gets the list of sample names.
    std::vector<std::string> VCFHeader::sample_names() const {
        std::vector<std::string> names;
        int n = n_samples();
        if (!is_valid() || n == 0) {
            return names;
        }
        // Access samples directly from the header struct
        names.reserve(n);
        for (int i = 0; i < n; ++i) {
            if (_hdr->samples[i]) { // Check for null pointers in the list
                 names.push_back(std::string(_hdr->samples[i]));
            }
        }
        return names;
    }

    // Gets the index of a sample by name.
    int VCFHeader::sample_index(const std::string& sample_name) const {
        if (!is_valid()) return -1;
        // bcf_hdr_id2int returns -1 if not found, which matches our desired return value.
        return bcf_hdr_id2int(_hdr.get(), BCF_DT_SAMPLE, sample_name.c_str());
    }

    // Adds a sample to the header.
    int VCFHeader::add_sample(const std::string& sample_name, bool warn) {
        if (!is_valid()) return -1;
        // bcf_hdr_add_sample returns 0 on success, <0 on failure.
        // It handles warnings internally if warn is true.
        return bcf_hdr_add_sample(_hdr.get(), sample_name.c_str());
        // Note: If warn is true (default in htslib), it prints to stderr on duplicate.
        // If warn is false, it returns -1 on duplicate.
    }

    // Removes a sample from the header by index.
    void VCFHeader::remove_sample(int index) {
        if (!is_valid() || index < 0 || index >= n_samples()) {
             throw std::out_of_range("Sample index out of range for removal.");
        }
        // Get the sample name first
        const char* sample_name = (index < n_samples() && _hdr->samples[index]) ? _hdr->samples[index] : nullptr;
        if (!sample_name) {
             throw std::runtime_error("Could not retrieve sample name for index " + std::to_string(index));
        }
        // bcf_hdr_remove takes type (BCF_HL_*) and the key (sample name for BCF_HL_GEN)
        bcf_hdr_remove(_hdr.get(), BCF_HL_GEN, sample_name);
        // Need to rebuild sample mapping after removal
        if (bcf_hdr_sync(_hdr.get()) < 0) {
             throw std::runtime_error("Failed to sync header after sample removal.");
        }
    }

    // Removes a sample from the header by name.
    void VCFHeader::remove_sample(const std::string& sample_name) {
        int index = sample_index(sample_name);
        if (index < 0) {
            throw std::runtime_error("Sample '" + sample_name + "' not found for removal.");
        }
        remove_sample(index);
    }

    // Subsets the header to include only the specified samples.
    VCFHeader VCFHeader::subset_samples(const std::vector<std::string>& samples_to_keep) const {
        if (!is_valid()) {
            throw std::runtime_error("Cannot subset an invalid VCF header.");
        }

        // Check if all requested samples exist
        // Get current sample names directly
        std::vector<std::string> current_sample_names = sample_names();

        for (const auto& name : samples_to_keep) {
            // Use std::find on the vector of strings
            if (std::find(current_sample_names.begin(), current_sample_names.end(), name) == current_sample_names.end()) {
                 throw std::runtime_error("Sample '" + name + "' requested for subsetting not found in original header.");
            }
        }

        // Create an array of C-style strings for bcf_hdr_subset
        std::vector<std::vector<char>> char_copies;
        std::vector<char*> samples_c;
        char_copies.reserve(samples_to_keep.size());
        samples_c.reserve(samples_to_keep.size());

        for (const auto& s : samples_to_keep) {
            char_copies.push_back(std::vector<char>(s.begin(), s.end()));
            char_copies.back().push_back('\0');  // 添加字符串终止符
            samples_c.push_back(char_copies.back().data());
        }
        std::vector<int> imap(samples_c.size());

        bcf_hdr_t* subset_hdr = bcf_hdr_subset(_hdr.get(), samples_c.size(), samples_c.data(), imap.data());
        if (!subset_hdr) {
            throw std::runtime_error("Failed to create subset VCF header.");
        }

        return VCFHeader(subset_hdr); // Return new VCFHeader managing the subset header
    }

    // Gets the sequence (contig) name for a given ID.
    std::string VCFHeader::seq_name(int id) const {
        if (!is_valid() || id < 0) return "";
        const char* name = bcf_hdr_id2name(_hdr.get(), id);
        return name ? std::string(name) : "";
    }

    // Gets the sequence (contig) ID for a given name.
    int VCFHeader::seq_id(const std::string& name) const {
        if (!is_valid()) return -1;
        // bcf_hdr_name2id returns -1 if not found.
        return bcf_hdr_name2id(_hdr.get(), name.c_str());
    }

    // Gets the length of a sequence (contig) by ID.
    hts_pos_t VCFHeader::seq_length(int id) const {
        if (!is_valid() || id < 0) return 0;

        // Get the header record (hrec) for the contig ID
        const bcf_hrec_t *hrec = bcf_hdr_get_hrec(_hdr.get(), BCF_HL_CTG, "ID", bcf_hdr_id2name(_hdr.get(), id), NULL);
        if (!hrec) return 0; // Contig ID not found or not a contig line

        // Find the 'length' key within the hrec
        for (int i = 0; i < hrec->nkeys; ++i) {
            if (strcmp(hrec->keys[i], "length") == 0) {
                // Found the length key, parse its value
                char *endptr;
                hts_pos_t len = strtoll(hrec->vals[i], &endptr, 10);
                // Check if parsing was successful (endptr should point to the null terminator)
                if (*endptr == '\0') {
                    return len;
                } else {
                    // Parsing failed or value format unexpected
                    return 0;
                }
            }
        }

        // Length key not found
        return 0;
    }

    // Gets the length of a sequence (contig) by name.
    hts_pos_t VCFHeader::seq_length(const std::string& name) const {
        int id = seq_id(name);
        return (id >= 0) ? seq_length(id) : 0;
    }

    // Gets the number of sequences (contigs) defined in the header.
    int VCFHeader::n_seqs() const {
        return is_valid() ? _hdr->n[BCF_DT_CTG] : 0;
    }

    // Gets the names of all sequences (contigs).
    std::vector<std::string> VCFHeader::seq_names() const {
        std::vector<std::string> names;
        if (!is_valid()) return names;
        int n = n_seqs();
        names.reserve(n);
        for (int i = 0; i < n; ++i) {
            const char* name = bcf_hdr_id2name(_hdr.get(), i);
            if (name) {
                names.push_back(std::string(name));
            }
        }
        return names;
    }

    // Adds a metadata line to the header.
    int VCFHeader::add_header_line(const std::string& line, int* id) {
        if (!is_valid()) return -1;
        // bcf_hdr_append adds the line but doesn't parse it immediately.
        // bcf_hdr_parse_line is better for adding structured lines.
        // Let's use bcf_hdr_append for simplicity, assuming the user provides a valid line.
        // For more robust parsing, one might use bcf_hdr_parse_line or specific add functions.
        int ret = bcf_hdr_append(_hdr.get(), line.c_str());
        if (ret != 0) return ret; // Error occurred

        // Need to sync header after appending lines for IDs to be updated
        ret = bcf_hdr_sync(_hdr.get());
        if (ret != 0) return ret; // Error syncing

        // If an ID pointer is provided, try to find the ID of the added line.
        // This is tricky with just bcf_hdr_append. A more complex approach
        // involving parsing the line first might be needed to reliably get the ID.
        // For now, we won't attempt to return the ID reliably via this method.
        if (id) *id = -1; // Indicate ID not reliably determined here

        return 0; // Success
    }

    // Formats the entire header into a string.
    std::string VCFHeader::to_string() const {
        if (!is_valid()) return "";
        kstring_t str = {0, 0, nullptr};
        if (bcf_hdr_format(_hdr.get(), 0, &str) < 0) {
            return "";
        }
        std::string result = (str.s ? std::string(str.s) : "");
        free(str.s);
        return result;
    }

    // Friend function for output stream
    std::ostream& operator<<(std::ostream& os, const VCFHeader& hdr) {
        os << hdr.to_string();
        return os;
    }

} // namespace ngslib