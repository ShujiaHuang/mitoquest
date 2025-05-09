// The C++ codes for VCF/BCF Header
// Author: Shujia Huang
// Date: 2025-04-17
#ifndef __INCLUDE_NGSLIB_VCF_HEADER_H__
#define __INCLUDE_NGSLIB_VCF_HEADER_H__

#include <iostream>
#include <string>
#include <vector>
#include <memory>    // For std::shared_ptr
#include <cstring>   // For strcmp
#include <stdexcept> // For std::runtime_error
#include <algorithm> // For std::find

#include <htslib/kstring.h> // For kstring_t used in bcf_hdr_fmt_text
#include <htslib/vcf.h>
#include <htslib/hts.h> // For hts_version()

#include "utils.h" // For join function

namespace ngslib {

    /**
     * @brief A C++ wrapper class for the htslib VCF header (bcf_hdr_t).
     *
     * This class manages the VCF header information, including metadata lines,
     * contig information, INFO/FORMAT field definitions, and sample names.
     * It provides methods to access, modify, and manage the header data.
     * The underlying bcf_hdr_t is managed by a shared_ptr for safe memory handling.
     */
    class VCFHeader {
    private:
        // Use shared_ptr to manage the lifetime of bcf_hdr_t.
        // This allows multiple VCFHeader objects (e.g., copied ones)
        // to safely share the same underlying htslib header structure.
        // bcf_hdr_destroy is used as the custom deleter.
        std::shared_ptr<bcf_hdr_t> _hdr;

        // Private constructor for internal use (e.g., by VcfReader/Writer)
        // explicit VCFHeader(bcf_hdr_t *hdr_ptr);

    public:
        /**
         * @brief Default constructor. Creates an empty header.
         */
        VCFHeader();
        VCFHeader(bcf_hdr_t *hdr_ptr);  // Constructor for internal use (e.g. by VcfReader)

        /**
         * @brief Copy constructor. Creates a new VCFHeader sharing the underlying bcf_hdr_t.
         * @param other The VCFHeader object to copy from.
         */
        VCFHeader(const VCFHeader& other);

        /**
         * @brief Move constructor.
         * @param other The VCFHeader object to move from.
         */
        VCFHeader(VCFHeader&& other) noexcept;

        /**
         * @brief Copy assignment operator.
         * @param other The VCFHeader object to assign from.
         * @return Reference to this VCFHeader object.
         */
        VCFHeader& operator=(const VCFHeader& other);

        /**
         * @brief Move assignment operator.
         * @param other The VCFHeader object to move assign from.
         * @return Reference to this VCFHeader object.
         */
        VCFHeader& operator=(VCFHeader&& other) noexcept;

        /**
         * @brief Destructor. The underlying bcf_hdr_t is automatically destroyed
         *        when the last shared_ptr pointing to it goes out of scope.
         */
        ~VCFHeader() = default; // Rely on shared_ptr for cleanup

        /**
         * @brief Checks if the header is valid (i.e., not null).
         * @return True if the header is valid, false otherwise.
         */
        bool is_valid() const { return _hdr != nullptr; }
        operator bool() const { return is_valid(); }

        /**
         * @brief Provides direct access to the underlying bcf_hdr_t pointer.
         * Use with caution, as direct manipulation can affect shared instances.
         * @return A const pointer to the bcf_hdr_t struct.
         */
        bcf_hdr_t* hts_header() const { return _hdr.get(); }

        /**
         * @brief Provides direct access to the underlying bcf_hdr_t pointer (non-const).
         * Use with extreme caution, especially if the header is shared.
         * Consider using `copy_header()` if modifications are needed without
         * affecting other shared instances.
         * @return A pointer to the bcf_hdr_t struct.
         */
        bcf_hdr_t* hts_header_writable() { return _hdr.get(); }

        /**
         * @brief Creates a deep copy of the underlying htslib header.
         * Useful when modifications are needed without affecting other VCFHeader
         * instances sharing the original header.
         * @return A new VCFHeader object with a duplicated header.
         */
        VCFHeader copy_header() const;

        /**
         * @brief Gets the number of samples in the header.
         * @return The number of samples, or 0 if the header is invalid.
         */
        int n_samples() const;

        /**
         * @brief Gets the list of sample names.
         * @return A vector of strings containing the sample names. Returns an empty vector if no samples or invalid header.
         */
        std::vector<std::string> sample_names() const;

        /**
         * @brief Gets the index of a sample by name.
         * @param sample_name The name of the sample.
         * @return The 0-based index of the sample, or -1 if not found or header is invalid.
         */
        int sample_index(const std::string& sample_name) const;

        /**
         * @brief Adds a sample to the header.
         * Note: This modifies the underlying header. Use `copy_header()` first
         * if you need to preserve the original shared header.
         * @param sample_name The name of the sample to add.
         * @param warn Warn if sample already exists (passed to bcf_hdr_add_sample).
         * @return 0 on success, -1 on failure (e.g., duplicate name, memory allocation error).
         */
        int add_sample(const std::string& sample_name, bool warn = true);

        /**
         * @brief Removes a sample from the header by index.
         * Note: This modifies the underlying header. Use `copy_header()` first
         * if you need to preserve the original shared header.
         * @param index The 0-based index of the sample to remove.
         */
        void remove_sample(int index);

        /**
         * @brief Removes a sample from the header by name.
         * Note: This modifies the underlying header. Use `copy_header()` first
         * if you need to preserve the original shared header.
         * @param sample_name The name of the sample to remove.
         */
        void remove_sample(const std::string& sample_name);

        /**
         * @brief Subsets the header to include only the specified samples.
         * Creates a new header containing only the specified samples and associated data.
         * @param samples_to_keep A vector of sample names to retain.
         * @return A new VCFHeader object containing only the specified samples.
         * @throws std::runtime_error if a sample in samples_to_keep is not found in the original header.
         */
        VCFHeader subset_samples(const std::vector<std::string>& samples_to_keep) const;

        /**
         * @brief Gets the sequence (contig) name for a given ID.
         * @param id The sequence ID (rid).
         * @return The sequence name, or an empty string if not found or header is invalid.
         */
        std::string seq_name(int id) const;

        /**
         * @brief Gets the sequence (contig) ID for a given name.
         * @param name The sequence name.
         * @return The sequence ID (rid), or -1 if not found or header is invalid.
         */
        int seq_id(const std::string& name) const;

        /**
         * @brief Gets the length of a sequence (contig) by ID.
         * @param id The sequence ID (rid).
         * @return The sequence length, or 0 if not found or header is invalid.
         */
        hts_pos_t seq_length(int id) const;

        /**
         * @brief Gets the length of a sequence (contig) by name.
         * @param name The sequence name.
         * @return The sequence length, or 0 if not found or header is invalid.
         */
        hts_pos_t seq_length(const std::string& name) const;

        /**
         * @brief Gets the number of sequences (contigs) defined in the header.
         * @return The number of sequences, or 0 if header is invalid.
         */
        int n_seqs() const;

        /**
         * @brief Gets the names of all sequences (contigs).
         * @return A vector of strings containing the sequence names.
         */
        std::vector<std::string> seq_names() const;

        /**
         * @brief Adds a metadata line to the header (e.g., INFO, FORMAT, FILTER, contig).
         * Note: This modifies the underlying header. Use `copy_header()` first
         * if you need to preserve the original shared header.
         * @param line The metadata line string (must be correctly formatted).
         * @param id Pointer to store the ID of the newly added line (optional).
         * @return 0 on success, < 0 on error.
         */
        int add_header_line(const std::string& line, int* id = nullptr);

        /**
         * @brief Formats the entire header into a string.
         * @return A string representation of the VCF header. Returns empty string if header is invalid.
         */
        std::string to_string() const;

        /**
         * @brief Gets the htslib version used.
         * @return A string containing the htslib version.
         */
        static std::string htslib_version() {
            return std::string(hts_version());
        }

        /**
         * @brief Get raw pointer for advanced usage (use with caution).
         * @return A pointer to the bcf_hdr_t struct.
         */
        bcf_hdr_t* raw_header() const { return _hdr.get(); }

        // Friend function for output stream
        friend std::ostream& operator<<(std::ostream& os, const VCFHeader& hdr);
        
    }; // class VCFHeader
} // namespace ngslib

#endif // __INCLUDE_NGSLIB_VCF_HEADER_H__