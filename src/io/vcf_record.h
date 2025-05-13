// The C++ codes for VCF/BCF Record
// Author: Shujia Huang
// Date: 2025-04-17
#ifndef __INCLUDE_NGSLIB_VCF_RECORD_H__
#define __INCLUDE_NGSLIB_VCF_RECORD_H__

#include <iostream>
#include <string>
#include <vector>
#include <memory> // For std::shared_ptr
#include <stdexcept>
#include <limits> // For numeric_limits

#include <htslib/vcf.h>
#include "vcf_header.h" // Needs header context
#include "utils.h"      // For potential utility functions

namespace ngslib {

    /**
     * @brief A C++ wrapper class for an htslib VCF/BCF record (bcf1_t).
     *
     * This class manages a single VCF/BCF record, providing methods to access
     * and modify its fields like CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO,
     * and FORMAT/sample data.
     * The underlying bcf1_t is managed by a shared_ptr for safe memory handling,
     * allowing records to be copied and shared efficiently.
     */
    class VCFRecord {
    private:
        // Use shared_ptr to manage the lifetime of bcf1_t.
        // bcf_destroy is used as the custom deleter.
        std::shared_ptr<bcf1_t> _b;

        // Helper to check if the record is valid (shared_ptr is not null)
        bool is_valid_unsafe() const { return _b != nullptr; }

    public:
        // Constants for missing values (matching htslib conventions)
        // Note: float constants cannot be constexpr as bcf_float_* might be functions/macros
        static constexpr int INT_MISSING    = bcf_int32_missing;
        static constexpr int INT_VECTOR_END = bcf_int32_vector_end;
        static const float FLOAT_MISSING;        // Initialized in cpp file
        static const float FLOAT_VECTOR_END;     // Initialized in cpp file
        static const std::string STRING_MISSING; // Initialized in cpp file

        /**
         * @brief Default constructor. Creates an empty, invalid record.
         */
        VCFRecord();
        VCFRecord(bcf1_t *b_ptr); // constructor for internal use (e.g., by VcfReader)

        /**
         * @brief Copy constructor. Shares ownership of the underlying bcf1_t.
         * @param other The VCFRecord object to copy from.
         */
        VCFRecord(const VCFRecord& other);

        /**
         * @brief Move constructor.
         * @param other The VCFRecord object to move from.
         */
        VCFRecord(VCFRecord&& other) noexcept;

        /**
         * @brief Copy assignment operator. Shares ownership.
         * @param other The VCFRecord object to assign from.
         * @return Reference to this VCFRecord object.
         */
        VCFRecord& operator=(const VCFRecord& other);

        /**
         * @brief Move assignment operator.
         * @param other The VCFRecord object to move assign from.
         * @return Reference to this VCFRecord object.
         */
        VCFRecord& operator=(VCFRecord&& other) noexcept;

        /**
         * @brief Destructor. The underlying bcf1_t is automatically destroyed
         *        when the last shared_ptr pointing to it goes out of scope.
         */
        ~VCFRecord() = default; // Rely on shared_ptr for cleanup

        /**
         * @brief Checks if the record is valid (i.e., points to an allocated bcf1_t).
         * @return True if the record is valid, false otherwise.
         */
        bool is_valid() const { return _b != nullptr; }
        operator bool() const { return is_valid(); }

        /**
         * @brief Provides direct access to the underlying bcf1_t pointer.
         * Use with caution, as direct manipulation can affect shared instances.
         * @return A const pointer to the bcf1_t struct. Returns nullptr if invalid.
         */
        bcf1_t* hts_record() const { return _b.get(); }

        /**
         * @brief Provides direct access to the underlying bcf1_t pointer (non-const).
         * Use with extreme caution. Modifying the record directly affects all
         * shared copies. Consider using `copy_record()` first.
         * @return A pointer to the bcf1_t struct. Returns nullptr if invalid.
         */
        bcf1_t* hts_record_writable() { return _b.get(); }

        /**
         * @brief Creates a deep copy of the underlying htslib record.
         * Useful when modifications are needed without affecting other VCFRecord
         * instances sharing the original record.
         * @return A new VCFRecord object with a duplicated record. Returns an invalid record if this one is invalid.
         */
        VCFRecord copy_record() const;

        /**
         * @brief Clears the record, making it ready for reuse (like bcf_clear).
         * Note: This modifies the underlying record. Use `copy_record()` first
         * if you need to preserve the original shared record.
         */
        void clear();

        /**
         * @brief Unpacks specific fields in the BCF record for access.
         * Call this before accessing fields like FILTER, INFO, or FORMAT.
         * @param which Fields to unpack (e.g., BCF_UN_ALL, BCF_UN_FLT, BCF_UN_INFO, BCF_UN_SHR, BCF_UN_FMT).
         *              BCF_UN_SHR unpacks shared fields (FILTER, INFO).
         *              BCF_UN_FMT unpacks FORMAT fields.
         *              BCF_UN_ALL unpacks all.
         * @return 0 on success, negative on error.
         */
        int unpack(int which) const;

        // === Accessors for Core Fields ===

        /**
         * @brief Gets the chromosome ID (rid). Use VCFHeader::seq_name to get the name.
         * @param hdr The VCFHeader associated with this record.
         * @return The chromosome ID, or -1 if invalid.
         */
        int32_t rid(const VCFHeader& hdr) const;

        /**
         * @brief Gets the chromosome name (CHROM).
         * @param hdr The VCFHeader associated with this record.
         * @return The chromosome name string, or empty string if invalid.
         */
        std::string chrom(const VCFHeader& hdr) const;

        /**
         * @brief Gets the 0-based position (POS).
         * @return The position, or -1 if invalid.
         */
        hts_pos_t pos() const;

        /**
         * @brief Gets the record ID (ID field). Multiple IDs are semicolon-separated.
         * @return The ID string, or "." if missing or invalid.
         */
        std::string id() const;

        /**
         * @brief Gets the reference allele (REF).
         * @return The reference allele string, or empty string if invalid.
         */
        std::string ref() const;

        /**
         * @brief Gets the alternate alleles (ALT).
         * @return A vector of strings containing the alternate alleles. Returns empty vector if invalid or no ALTs.
         */
        std::vector<std::string> alt() const;

        /**
         * @brief Gets the number of alternate alleles.
         * @return The number of alternate alleles, or 0 if invalid.
         */
        int n_alt() const;

        /**
         * @brief Gets the quality score (QUAL).
         * @return The quality score. Returns bcf_float_missing if missing or invalid.
         */
        float qual() const;

        /**
         * @brief Gets the number of filters applied. Requires unpack(BCF_UN_FLT).
         * @return The number of filters, or -1 if invalid or not unpacked.
         */
        int n_filter() const;

        /**
         * @brief Gets the filter IDs. Requires unpack(BCF_UN_FLT).
         * Use VCFHeader::id2int(BCF_DT_ID, filter_id) to get the string name.
         * @return A vector of integer filter IDs. Returns empty vector if invalid, no filters, or not unpacked.
         */
        std::vector<int> filter_ids() const;

        /**
         * @brief Gets the filter names as strings. Requires unpack(BCF_UN_FLT).
         * @param hdr The VCFHeader associated with this record.
         * @return A vector of strings containing the filter names (e.g., "PASS", "q10").
         */
        std::vector<std::string> filter_names(const VCFHeader& hdr) const;

        /**
         * @brief Checks if the record has the "PASS" filter or no filters applied. Requires unpack(BCF_UN_FLT).
         * @param hdr The VCFHeader associated with this record.
         * @return True if the record passed filters, false otherwise.
         */
        bool passed_filters(const VCFHeader& hdr) const;

        /**
         * @brief Gets the number of INFO fields present in the record. Requires unpack(BCF_UN_INFO).
         * @return The number of INFO fields, or 0 if invalid or not unpacked.
         */
        int n_info() const;

        /**
         * @brief Gets the number of FORMAT fields present. Requires unpack(BCF_UN_FMT).
         * @return The number of FORMAT fields, or 0 if invalid or not unpacked.
         */
        int n_format() const;

        /**
         * @brief Gets the number of samples for which FORMAT data exists.
         * This should match VCFHeader::n_samples().
         * @return The number of samples, or 0 if invalid.
         */
        int n_samples() const;


        // === INFO Field Accessors ===
        // Note: These require unpack(BCF_UN_INFO) to have been called.
        // They return default/missing values if the tag is not present or the record is invalid.

        /**
         * @brief Gets an integer INFO tag value.
         * @param hdr VCFHeader for tag ID lookup.
         * @param tag The INFO tag name (e.g., "DP").
         * @param values Vector to store the retrieved integer values. Cleared before filling.
         * @return The number of values retrieved (>=0), or a negative value on error (e.g., tag not found, type mismatch).
         *         BCF_ERR_TAG_UNDEF (-1): Tag undefined in header
         *         BCF_ERR_TAG_INVALID (-2): Tag defined for different type
         *         BCF_ERR_N_VALUES (-3): Incorrect number of values
         */
        int get_info_int(const VCFHeader& hdr, const std::string& tag, std::vector<int32_t>& values) const;

        /**
         * @brief Gets a float INFO tag value.
         * @param hdr VCFHeader for tag ID lookup.
         * @param tag The INFO tag name (e.g., "AF").
         * @param values Vector to store the retrieved float values. Cleared before filling.
         * @return The number of values retrieved (>=0), or a negative value on error.
         */
        int get_info_float(const VCFHeader& hdr, const std::string& tag, std::vector<float>& values) const;

        /**
         * @brief Gets a flag INFO tag value.
         * @param hdr VCFHeader for tag ID lookup.
         * @param tag The INFO flag name (e.g., "DB").
         * @return 1 if the flag is present, 0 if not present, negative on error (e.g., tag not defined).
         */
        int get_info_flag(const VCFHeader& hdr, const std::string& tag) const;

        /**
         * @brief Gets a string INFO tag value.
         * @param hdr VCFHeader for tag ID lookup.
         * @param tag The INFO tag name (e.g., "AnnoType").
         * @param values Vector to store the retrieved string values. Cleared before filling.
         * @return The number of values retrieved (>=0), or a negative value on error.
         */
        int get_info_string(const VCFHeader& hdr, const std::string& tag, std::vector<std::string>& values) const;


        // === FORMAT Field Accessors ===
        // Note: These require unpack(BCF_UN_FMT) to have been called.
        // They return default/missing values if the tag is not present, the sample index is invalid, or the record is invalid.

        /**
         * @brief Gets genotype (GT) FORMAT tag values for all samples.
         * Special helper function for the common GT field.
         * @param hdr VCFHeader for tag ID lookup.
         * @param genotypes Vector to store retrieved genotypes. Each inner vector represents a sample's genotype (e.g., {0, 1} for het). 
         *                  Missing alleles are represented by negative values (e.g., bcf_gt_missing). Phasing is stored separately 
         *                  if needed (use bcf_gt_is_phased).
         * @return The ploidy (number of alleles per sample, e.g., 2 for diploid), or a negative value on error.
         */
        int get_genotypes(const VCFHeader& hdr, std::vector<std::vector<int>>& genotypes) const;

        /**
         * @brief Gets integer FORMAT tag values for all samples.
         * @param hdr VCFHeader for tag ID lookup.
         * @param tag The FORMAT tag name (e.g., "GT", "DP").
         * @param values Vector to store the retrieved integer values. It will be resized to n_samples * number_per_sample. Values are stored sample-major (all values for sample 0, then all for sample 1, ...).
         * @return The number of values retrieved per sample (>=0), or a negative value on error.
         */
        int get_format_int(const VCFHeader& hdr, const std::string& tag, std::vector<int32_t>& values) const;

        /**
         * @brief Gets float FORMAT tag values for all samples.
         * @param hdr VCFHeader for tag ID lookup.
         * @param tag The FORMAT tag name (e.g., "PL", "GP").
         * @param values Vector to store the retrieved float values. Resized and stored sample-major.
         * @return The number of values retrieved per sample (>=0), or a negative value on error.
         */
        int get_format_float(const VCFHeader& hdr, const std::string& tag, std::vector<float>& values) const;

        /**
         * @brief Gets string FORMAT tag values for all samples.
         * @param hdr VCFHeader for tag ID lookup.
         * @param tag The FORMAT tag name (e.g., "FT").
         * @param values Vector to store the retrieved string values. Resized and stored sample-major.
         * @return The number of values retrieved per sample (>=0), or a negative value on error.
         */
        int get_format_string(const VCFHeader& hdr, const std::string& tag, std::vector<std::string>& values) const;

        /**
         * @brief Gets the index for a specific tag in the FORMAT field.
         * @param hdr The VCFHeader associated with this record.
         * @param tag The FORMAT tag name (e.g., "GT").
         * @return The ID of the FORMAT field, or -1 if not found.
         */
        int get_format_idx(const VCFHeader& hdr, const std::string& tag) const;

        /**
         * @brief Gets the ploidy for a specific sample. 
         * @param hdr The VCFHeader associated with this record.
         * @return The ploidy, or -1 if invalid.
        */
        int get_sample_ploidy(const VCFHeader& hdr, int sample_idx) const;

        /**
         * @brief Gets the ploidy for all samples.
         * @param hdr The VCFHeader associated with this record.
         * @return A vector of integers representing the ploidy for each sample. Returns empty vector if invalid.
         */
        std::vector<int> get_sample_ploidies(const VCFHeader& hdr) const;

        /**
         * @brief Gets the maximum ploidy across all samples.
         * @param hdr The VCFHeader associated with this record.
         * @return The maximum ploidy, or -1 if invalid.
         */
        int get_max_ploidy(const VCFHeader& hdr) const;


        // === Modifiers ===
        // Note: These modify the underlying record. Use copy_record() first if needed.

        /**
         * @brief Sets the position (0-based).
         * @param pos The new position.
         */
        void set_pos(hts_pos_t pos);

        /**
         * @brief Sets the reference allele.
         * @param ref The reference allele string.
         */
        void set_ref(const std::string& ref);

        /**
         * @brief Sets the reference and alternate alleles. Requires header context.
         * @param hdr The VCFHeader associated with this record.
         * @param alt_alleles Vector of alternate allele strings.
         */
        void set_alt(const VCFHeader& hdr, const std::vector<std::string>& alt_alleles);

        /**
         * @brief Sets the quality score.
         * @param qual The quality score. Use bcf_float_missing for missing QUAL.
         */
        void set_qual(float qual);

        /**
         * @brief Sets the ID field. Semicolon-separated if multiple IDs. Use "." for missing. Requires header context.
         * @param hdr The VCFHeader associated with this record.
         * @param id_str The ID string.
         */
        void set_id(const VCFHeader& hdr, const std::string& id_str);

        /**
         * @brief Adds a filter ID to the record. Requires unpack(BCF_UN_FLT).
         * @param hdr The VCFHeader containing the filter definition.
         * @param filter_name The name of the filter to add (must be defined in header).
         * @return 0 on success, negative on error (e.g., filter not defined).
         */
        int add_filter(const VCFHeader& hdr, const std::string& filter_name);

        /**
         * @brief Sets the filters for the record, overwriting existing ones. Requires unpack(BCF_UN_FLT).
         * @param hdr The VCFHeader containing the filter definitions.
         * @param filter_names Vector of filter names to set. Use an empty vector to clear filters. Use {"PASS"} for PASS.
         * @return 0 on success, negative on error.
         */
        int set_filters(const VCFHeader& hdr, const std::vector<std::string>& filter_names);

        /**
         * @brief Clears all filters from the record. Requires unpack(BCF_UN_FLT). Requires header context.
         * @param hdr The VCFHeader associated with this record.
         * @return 0 on success, negative on error.
         */
        int clear_filters(const VCFHeader& hdr);

        /**
         * @brief Updates an integer INFO tag. Requires unpack(BCF_UN_INFO).
         * @param hdr VCFHeader for tag ID lookup.
         * @param tag The INFO tag name.
         * @param values Pointer to the integer data array.
         * @param n_values Number of values in the array. Use 0 to remove the tag.
         * @return 0 on success, negative on error.
         */
        int update_info_int(const VCFHeader& hdr, const std::string& tag, const int32_t* values, int n_values);

        /**
         * @brief Updates a float INFO tag. Requires unpack(BCF_UN_INFO).
         * @param hdr VCFHeader for tag ID lookup.
         * @param tag The INFO tag name.
         * @param values Pointer to the float data array.
         * @param n_values Number of values in the array. Use 0 to remove the tag.
         * @return 0 on success, negative on error.
         */
        int update_info_float(const VCFHeader& hdr, const std::string& tag, const float* values, int n_values);

        /**
         * @brief Updates a flag INFO tag. Requires unpack(BCF_UN_INFO).
         * @param hdr VCFHeader for tag ID lookup.
         * @param tag The INFO flag name.
         * @param set Flag value (true to set, false to remove).
         * @return 0 on success, negative on error.
         */
        int update_info_flag(const VCFHeader& hdr, const std::string& tag, bool set);

        /**
         * @brief Updates a string INFO tag. Requires unpack(BCF_UN_INFO).
         * @param hdr VCFHeader for tag ID lookup.
         * @param tag The INFO tag name.
         * @param values Pointer to the C-string data. Use NULL or empty string to remove tag.
         * @return 0 on success, negative on error.
         */
        int update_info_string(const VCFHeader& hdr, const std::string& tag, const char* value);

        /**
         * @brief Updates integer FORMAT tag values. Requires unpack(BCF_UN_FMT).
         * @param hdr VCFHeader for tag ID lookup.
         * @param tag The FORMAT tag name (e.g., "DP").
         * @param values Pointer to the integer data array (sample-major order). Size must be n_samples * number_per_sample.
         * @param n_values_per_sample Number of values per sample for this tag.
         * @return 0 on success, negative on error.
         */
        int update_format_int(const VCFHeader& hdr, const std::string& tag, const int32_t* values, int n_values_per_sample);

        /**
         * @brief Updates float FORMAT tag values. Requires unpack(BCF_UN_FMT).
         * @param hdr VCFHeader for tag ID lookup.
         * @param tag The FORMAT tag name (e.g., "PL").
         * @param values Pointer to the float data array (sample-major order).
         * @param n_values_per_sample Number of values per sample for this tag.
         * @return 0 on success, negative on error.
         */
        int update_format_float(const VCFHeader& hdr, const std::string& tag, const float* values, int n_values_per_sample);

        /**
         * @brief Updates string FORMAT tag values. Requires unpack(BCF_UN_FMT).
         * @param hdr VCFHeader for tag ID lookup.
         * @param tag The FORMAT tag name (e.g., "FT").
         * @param values Array of C-strings (one per sample). Size must be n_samples.
         * @return 0 on success, negative on error.
         */
        int update_format_string(const VCFHeader& hdr, const std::string& tag, const char** values);

        /**
         * @brief Cleans up alleles (REF and ALT) to remove invalid or missing values.
         * This is a helper function for update_alleles.
         * @param hdr The VCFHeader associated with this record.
         * @return True if cleanup was successful, false if the record is invalid.
         */
        bool cleanup_alleles(const ngslib::VCFHeader& hdr);

        /**
         * @brief Updates the reference and alternate alleles. Requires header context.
         * @param hdr The VCFHeader associated with this record.
         * @param ref The new reference allele.
         * @param alts Vector of new alternate alleles.
         * @return 0 on success, negative on error.
         */
        int update_alleles(const VCFHeader& hdr, const std::string& ref, const std::vector<std::string>& alts);

        /**
         * @brief Updates genotype (GT) FORMAT tag values. Requires unpack(BCF_UN_FMT).
         * @param hdr VCFHeader for tag ID lookup.
         * @param genotypes Vector of vectors containing genotype data (sample-major, allele-major within sample). Use bcf_gt_unphased() etc. to construct values.
         * @return 0 on success, negative on error.
         */
        int update_genotypes(const VCFHeader& hdr, const std::vector<std::vector<int>>& genotypes);

        /**
         * @brief Subsets the record to include only specified samples
         * 
         * @param hdr The VCF header containing the subset information
         * @param samples_to_keep Vector of sample names to keep
         * @return VCFRecord A new record containing only the specified samples
         * @throws std::runtime_error if record is invalid or subsetting fails
         */
        VCFRecord subset_samples(const VCFHeader& hdr, const std::vector<std::string>& samples_to_keep) const;

        /**
         * @brief Subsets the record to include only specified samples (by index)
         * 
         * @param hdr The VCF header containing the subset information
         * @param sample_indices Vector of sample indices to keep
         * @return VCFRecord A new record containing only the specified samples
         * @throws std::runtime_error if record is invalid or subsetting fails
         */
        VCFRecord subset_samples(const VCFHeader& hdr, std::vector<int>& sample_indices) const;

        // Friend function for output stream (basic representation)
        friend std::ostream& operator<<(std::ostream& os, const VCFRecord& rec);

    }; // class VCFRecord

} // namespace ngslib

#endif // __INCLUDE_NGSLIB_VCF_RECORD_H__