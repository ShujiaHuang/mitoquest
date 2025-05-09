// The C++ codes for VCF/BCF Record
// Author: Roo (Inspired by Shujia Huang's BAM wrappers)
// Date: 2025-04-17

#include "vcf_record.h"
#include <htslib/vcf.h>
#include <htslib/kstring.h> // For string manipulation if needed
#include <vector>
#include <string>
#include <cstring> // For strcmp, strlen
#include <stdexcept>
#include <limits> // Required for numeric_limits

namespace ngslib {

    // Initialize static const members
    const std::string VCFRecord::STRING_MISSING = ".";
    // Initialize float constants here as they might not be true compile-time constants
    const float VCFRecord::FLOAT_MISSING = bcf_float_missing;
    const float VCFRecord::FLOAT_VECTOR_END = bcf_float_vector_end;

    // --- Constructors & Destructor ---

    // Private constructor
    VCFRecord::VCFRecord(bcf1_t *b_ptr) : _b(b_ptr, bcf_destroy) {
        // Takes ownership of b_ptr via shared_ptr with custom deleter
        if (!b_ptr) {
             // If a null pointer is passed, create an empty, invalid record.
             // This differs from VCFHeader where we created a minimal valid one.
             // For records, usually you get them from a reader, so null implies invalidity.
             _b.reset(); // Ensure shared_ptr is null
        }
    }

    // Default constructor
    VCFRecord::VCFRecord() : _b(nullptr) {} // Creates an invalid record

    // Copy constructor
    VCFRecord::VCFRecord(const VCFRecord& other) : _b(other._b) {} // Shares ownership

    // Move constructor
    VCFRecord::VCFRecord(VCFRecord&& other) noexcept : _b(std::move(other._b)) {}

    // Copy assignment operator
    VCFRecord& VCFRecord::operator=(const VCFRecord& other) {
        if (this != &other) {
            _b = other._b; // Share ownership
        }
        return *this;
    }

    // Move assignment operator
    VCFRecord& VCFRecord::operator=(VCFRecord&& other) noexcept {
        if (this != &other) {
            _b = std::move(other._b);
        }
        return *this;
    }

    // --- Core Methods ---

    VCFRecord VCFRecord::copy_record() const {
        if (!is_valid_unsafe()) return VCFRecord(); // Return invalid if source is invalid
        bcf1_t* dup_b = bcf_dup(_b.get());
        if (!dup_b) {
            throw std::runtime_error("Failed to duplicate VCF record.");
        }
        // Create a new VCFRecord taking ownership of the duplicated raw pointer
        return VCFRecord(dup_b);
    }

    void VCFRecord::clear() {
        if (is_valid_unsafe()) {
            bcf_clear(_b.get());
        }
    }

    int VCFRecord::unpack(int which) {
        if (!is_valid_unsafe()) return -1; // Or some other error code
        return bcf_unpack(_b.get(), which);
    }

    // --- Accessors for Core Fields ---
    VCFRecord VCFRecord::subset_samples(const VCFHeader& hdr, const std::vector<std::string>& samples_to_keep) const {
        if (!is_valid()) {
            throw std::runtime_error("Cannot subset an invalid VCF record");
        }
    
        // Create index mapping
        std::vector<int> sample_indices;
        sample_indices.reserve(samples_to_keep.size());
        
        for (const auto& name : samples_to_keep) {
            int idx = hdr.sample_index(name);
            if (idx < 0) {
                throw std::runtime_error("Sample '" + name + "' not found in VCF record");
            }
            sample_indices.push_back(idx);
        }
    
        return subset_samples(hdr, sample_indices);
    }
    
    VCFRecord VCFRecord::subset_samples(const VCFHeader& hdr, std::vector<int>& sample_indices) const {
        if (!is_valid()) {
            throw std::runtime_error("Cannot subset an invalid VCF record");
        }
    
        // Create a copy of the record for subsetting
        VCFRecord subset_rec = copy_record();
        
        // Call bcf_subset with indices
        // Note: bcf_subset modifies the record in place, so we need to pass the copy
        // Signature: int bcf_subset(const bcf_hdr_t *h, bcf1_t *v, int n, int *imap);
        if (bcf_subset(hdr.hts_header(), subset_rec._b.get(), sample_indices.size(), sample_indices.data()) != 0) {
            throw std::runtime_error("Failed to subset record at " + chrom(hdr) + ":" + std::to_string(pos()+1));
        }
    
        return subset_rec;
    }

    bool VCFRecord::update_alleles(const VCFHeader& hdr, const std::string& ref, const std::vector<std::string>& alts) {
        if (!is_valid_unsafe() || !hdr.is_valid()) return false;
        
        // 确保记录已解包
        if (!(_b->unpacked & BCF_UN_STR)) {
            if (bcf_unpack(_b.get(), BCF_UN_STR) < 0) {
                return false;
            }
        }
    
        // 计算所有等位基因的数量（REF + ALTs）
        const int n_alleles = 1 + alts.size();
        
        // 创建一个指针数组来存储所有等位基因字符串
        std::vector<const char*> alleles(n_alleles);
        
        // 设置 REF 等位基因
        alleles[0] = ref.c_str();
        
        // 设置 ALT 等位基因
        for (size_t i = 0; i < alts.size(); ++i) {
            alleles[i + 1] = alts[i].c_str();
        }
    
        // 一次性更新所有等位基因，使用传入的 header
        int ret = bcf_update_alleles(
            hdr.hts_header(),  // 使用传入的 header
            _b.get(),
            alleles.data(),
            n_alleles
        );
    
        return (ret >= 0);
    }

    int32_t VCFRecord::rid(const VCFHeader& hdr) const {
        if (!is_valid_unsafe() || !hdr.is_valid()) return -1;
        return _b->rid;
    }

    std::string VCFRecord::chrom(const VCFHeader& hdr) const {
        if (!is_valid_unsafe() || !hdr.is_valid()) return "";
        return hdr.seq_name(_b->rid);
    }

    hts_pos_t VCFRecord::pos() const {
        // POS is 0-based internally in bcf1_t
        return is_valid_unsafe() ? _b->pos : -1;
    }

    std::string VCFRecord::id() const {
        if (!is_valid_unsafe()) return STRING_MISSING;
        // Ensure shared fields are unpacked before accessing id
        // Although bcf_get_id doesn't strictly require it, it's good practice
        // bcf_unpack(_b.get(), BCF_UN_SHR); // Caller should ensure this if needed
        return (_b->d.id ? std::string(_b->d.id) : STRING_MISSING);
    }

    std::string VCFRecord::ref() const {
        if (!is_valid_unsafe() || _b->n_allele == 0) return "";
        // bcf_unpack(_b.get(), BCF_UN_STR); // Alleles are usually available without unpack
        return (_b->d.allele[0] ? std::string(_b->d.allele[0]) : "");
    }

    std::vector<std::string> VCFRecord::alt() const {
        std::vector<std::string> alleles;
        if (!is_valid_unsafe() || _b->n_allele <= 1) return alleles;
        // bcf_unpack(_b.get(), BCF_UN_STR); // Alleles usually available
        alleles.reserve(_b->n_allele - 1);
        for (int i = 1; i < _b->n_allele; ++i) {
             alleles.push_back(_b->d.allele[i] ? std::string(_b->d.allele[i]) : "");
        }
        return alleles;
    }

     int VCFRecord::n_alt() const {
         return is_valid_unsafe() ? (_b->n_allele > 0 ? _b->n_allele - 1 : 0) : 0;
     }

    float VCFRecord::qual() const {
        if (!is_valid_unsafe()) return FLOAT_MISSING; // Or std::numeric_limits<float>::quiet_NaN();
        return bcf_float_is_missing(_b->qual) ? FLOAT_MISSING : _b->qual;
    }

    int VCFRecord::n_filter() const {
        // Requires BCF_UN_FLT
        if (!is_valid_unsafe()) return -1;
        return _b->d.n_flt; // Returns number of filters, or -1 if not unpacked
    }

    std::vector<int> VCFRecord::filter_ids() const {
        std::vector<int> ids;
        int n = n_filter(); // Checks validity and unpack status implicitly
        if (n <= 0) return ids; // No filters or not unpacked/invalid
        ids.assign(_b->d.flt, _b->d.flt + n);
        return ids;
    }

    std::vector<std::string> VCFRecord::filter_names(const VCFHeader& hdr) const {
        std::vector<std::string> names;
        if (!hdr.is_valid()) return names;
        std::vector<int> ids = filter_ids(); // Checks validity/unpack
        names.reserve(ids.size());
        for (int id : ids) {
            const char* name = bcf_hdr_int2id(hdr.hts_header(), BCF_DT_ID, id);
            names.push_back(name ? std::string(name) : "?"); // Use "?" for unknown IDs
        }
        return names;
    }

     bool VCFRecord::passed_filters(const VCFHeader& hdr) const {
        int n = n_filter();
        if (n < 0) return false; // Not unpacked or invalid
        if (n == 0) return true;  // No filters applied means PASS

        // Check if the only filter is PASS
        if (n == 1) {
            const char* name = bcf_hdr_int2id(hdr.hts_header(), BCF_DT_ID, _b->d.flt[0]);
            return (name && strcmp(name, "PASS") == 0);
        }
        // Multiple filters, or filters other than PASS
        return false;
    }

    int VCFRecord::n_info() const {
        // Requires BCF_UN_INFO
        if (!is_valid_unsafe()) return 0;
        // Access internal structure - check if info is populated
        return _b->n_info;
    }

    int VCFRecord::n_format() const {
        // Requires BCF_UN_FMT
        if (!is_valid_unsafe()) return 0;
        return _b->n_fmt;
    }

    int VCFRecord::n_samples() const {
        // Doesn't require unpack, gets from core struct
        return is_valid_unsafe() ? _b->n_sample : 0;
    }

    // --- INFO Field Accessors ---

    int VCFRecord::get_info_int(const VCFHeader& hdr, const std::string& tag, std::vector<int32_t>& values) const {
        values.clear();
        if (!is_valid_unsafe() || !hdr.is_valid()) return BCF_ERR_TAG_INVALID; // Indicate error

        int32_t* buffer = nullptr;
        int n_values = 0; // htslib expects pointer to int for count
        int ret = bcf_get_info_int32(hdr.hts_header(), _b.get(), tag.c_str(), &buffer, &n_values);

        if (ret < 0) { // Error or tag not found
            if (buffer) free(buffer);
            return ret; // Return htslib error code
        }

        if (n_values > 0 && buffer) {
            values.assign(buffer, buffer + n_values);
        }
        if (buffer) free(buffer);
        return n_values; // Return number of values read
    }

    int VCFRecord::get_info_float(const VCFHeader& hdr, const std::string& tag, std::vector<float>& values) const {
        values.clear();
        if (!is_valid_unsafe() || !hdr.is_valid()) return BCF_ERR_TAG_INVALID;

        float* buffer = nullptr;
        int n_values = 0;
        int ret = bcf_get_info_float(hdr.hts_header(), _b.get(), tag.c_str(), &buffer, &n_values);

        if (ret < 0) {
            if (buffer) free(buffer);
            return ret;
        }

        if (n_values > 0 && buffer) {
            values.assign(buffer, buffer + n_values);
        }
        if (buffer) free(buffer);
        return n_values;
    }

    int VCFRecord::get_info_flag(const VCFHeader& hdr, const std::string& tag) const {
        if (!is_valid_unsafe() || !hdr.is_valid()) return BCF_ERR_TAG_INVALID; // Indicate error

        // Need to unpack INFO first
        if (!(_b->unpacked & BCF_UN_INFO)) {
             // Attempt to unpack if not already done (const_cast needed for unpack)
             // Note: Modifying unpack status in a const method is tricky.
             // Ideally, unpack should be called by the user beforehand.
             // We'll proceed assuming it might be unpacked, but this isn't ideal const-correctness.
             // A better design might involve a non-const getter or requiring prior unpack.
             if (bcf_unpack(const_cast<bcf1_t*>(_b.get()), BCF_UN_INFO) != 0) {
                 // Consider logging a warning here instead of returning an error code,
                 // as the failure might be due to the const context.
                 // For now, return an error code consistent with htslib.
                 return BCF_ERR_TAG_INVALID; // Failed to unpack
             }
        }

        int tag_id = bcf_hdr_id2int(hdr.hts_header(), BCF_DT_ID, tag.c_str());
        if (tag_id < 0) return 0; // Tag not defined in header means flag is not present

        // Iterate through INFO fields to find the flag
        for (int i = 0; i < _b->n_info; ++i) {
            bcf_info_t *info = &_b->d.info[i];
            if (info->key == tag_id) {
                // Check if it's actually a flag type (length 0 or 1 with type BCF_BT_NULL)
                 // Also check the type defined in the header
                 int type = bcf_hdr_id2type(hdr.hts_header(), BCF_HL_INFO, tag_id);
                 return (type == BCF_HT_FLAG) ? 1 : 0; // Return 1 if flag found and type matches, 0 otherwise
            }
        }
        return 0; // Flag not found in the record's INFO fields
    }

    int VCFRecord::get_info_string(const VCFHeader& hdr, const std::string& tag, std::vector<std::string>& values) const {
        values.clear();
        if (!is_valid_unsafe() || !hdr.is_valid()) return BCF_ERR_TAG_INVALID;

        char* buffer = nullptr;
        int buffer_size = 0; // htslib will realloc
        int n_values = 0; // For string, this often indicates if the tag exists and has length > 0

        // Use bcf_get_info_string. It returns length of string, or <0 on error.
        // It handles allocation.
        int ret = bcf_get_info_string(hdr.hts_header(), _b.get(), tag.c_str(), &buffer, &buffer_size);

        if (ret < 0) { // Error or tag not found
             if (buffer) free(buffer);
             return ret; // Return htslib error code (-1: undef, -2: type err, -3: missing)
        }

        if (ret > 0 && buffer) { // ret is the length of the string data
            // VCF strings can be comma-separated lists. We need to split them.
            // For simplicity here, we treat it as a single string value.
            // A more robust implementation would split by comma if Number != 1.
            values.push_back(std::string(buffer));
            n_values = 1; // Indicate one string value retrieved
        }
        // If ret == 0, it means empty string value, technically present.

        if (buffer) free(buffer);
        return n_values; // Return number of strings pushed (0 or 1 in this simple case)
    }


    // --- FORMAT Field Accessors ---
    int VCFRecord::get_genotypes(const VCFHeader& hdr, std::vector<std::vector<int>>& genotypes) const {
        genotypes.clear();
        if (!is_valid_unsafe() || !hdr.is_valid() || n_samples() == 0) 
            return BCF_ERR_TAG_INVALID;
    
        // 确保 FORMAT 字段已解包
        if (!((_b->unpacked & BCF_UN_FMT))) {
            if (bcf_unpack(_b.get(), BCF_UN_FMT) < 0) {
                return BCF_ERR_TAG_INVALID;
            }
        }
    
        // 查找 GT 字段
        int gt_fmt_id = -1;
        for (int i = 0; i < _b->n_fmt; ++i) {
            if (_b->d.fmt[i].id == bcf_hdr_id2int(hdr.hts_header(), BCF_DT_ID, "GT")) {
                gt_fmt_id = i;
                break;
            }
        }
        
        if (gt_fmt_id < 0) return BCF_ERR_TAG_INVALID;  // GT 字段不存在
    
        bcf_fmt_t *fmt = &_b->d.fmt[gt_fmt_id];
        int n_samp = n_samples();
        genotypes.resize(n_samp);
    
        // 处理每个样本
        for (int i = 0; i < n_samp; ++i) {
            // 获取当前样本的基因型数据起始位置, 这比使用 bcf_get_genotypes 更灵活，更适合处理可变倍性的情况
            uint8_t *curr_sample = fmt->p + i * fmt->size;

            // 根据最大倍性分配空间并填充数据
            genotypes[i].reserve(fmt->size);
            for (int j = 0; j < fmt->size; ++j) {
                // 处理基因型数据
                int allele = static_cast<int>(curr_sample[j]);
                if (allele == bcf_int8_vector_end + 256) break; // 处理向量结束,此处不可以用 continue       
                
                if (allele == bcf_int8_missing) { // 处理缺失值: '.'
                    genotypes[i].push_back(-1);   // 使用 -1 表示缺失
                } else {
                    genotypes[i].push_back(bcf_gt_allele(allele));  // 存等位基因值
                }
            }
        }

        return fmt->size; // max ploidy
    }

    int VCFRecord::get_format_int(const VCFHeader& hdr, const std::string& tag, std::vector<int32_t>& values) const {
        values.clear();
        if (!is_valid_unsafe() || !hdr.is_valid() || n_samples() == 0) return BCF_ERR_TAG_INVALID;

        int32_t* buffer = nullptr;
        int n_values_per_sample = 0; // htslib expects pointer to int
        int ret = bcf_get_format_int32(hdr.hts_header(), _b.get(), tag.c_str(), &buffer, &n_values_per_sample);

        if (ret < 0) {
            if (buffer) free(buffer);
            return ret; // Return htslib error code
        }

        int total_values = n_samples() * n_values_per_sample;
        if (total_values > 0 && buffer) {
            values.assign(buffer, buffer + total_values);
        }
        if (buffer) free(buffer);
        return n_values_per_sample; // Return number of values *per sample*
    }

    int VCFRecord::get_format_float(const VCFHeader& hdr, const std::string& tag, std::vector<float>& values) const {
        values.clear();
        if (!is_valid_unsafe() || !hdr.is_valid() || n_samples() == 0) return BCF_ERR_TAG_INVALID;

        float* buffer = nullptr;
        int n_values_per_sample = 0;
        int ret = bcf_get_format_float(hdr.hts_header(), _b.get(), tag.c_str(), &buffer, &n_values_per_sample);

        if (ret < 0) {
            if (buffer) free(buffer);
            return ret;
        }

        int total_values = n_samples() * n_values_per_sample;
        if (total_values > 0 && buffer) {
            values.assign(buffer, buffer + total_values);
        }
        if (buffer) free(buffer);
        return n_values_per_sample;
    }

    int VCFRecord::get_format_string(const VCFHeader& hdr, const std::string& tag, std::vector<std::string>& values) const {
        values.clear();
        if (!is_valid_unsafe() || !hdr.is_valid() || n_samples() == 0) return BCF_ERR_TAG_INVALID;

        char** buffer = nullptr; // Array of strings
        int n_values_per_sample = 0; // Number of strings per sample (usually 1)
        int ret = bcf_get_format_string(hdr.hts_header(), _b.get(), tag.c_str(), &buffer, &n_values_per_sample);

        if (ret < 0) {
            // Note: bcf_get_format_string allocates buffer[0] contiguously for all strings.
            // Only need to free buffer[0] and buffer itself.
            if (buffer) {
                if (buffer[0]) free(buffer[0]);
                free(buffer);
            }
            return ret;
        }

        int n_samp = n_samples();
        values.resize(n_samp * n_values_per_sample); // Resize the output vector

        if (buffer && buffer[0]) {
             // Data is in buffer[0], laid out contiguously, separated by terminators?
             // Or is buffer an array of char* pointers? The docs are a bit ambiguous.
             // Let's assume buffer is char** pointing to individual strings per sample.
             // This seems more consistent with bcf_update_format_string.
             for (int i = 0; i < n_samp; ++i) {
                 // Assuming n_values_per_sample is 1 for strings
                 if (buffer[i]) {
                     values[i] = std::string(buffer[i]);
                 } else {
                     values[i] = STRING_MISSING; // Or empty string
                 }
             }
        }

        if (buffer) {
            // Free according to htslib example for bcf_get_format_string
             if (buffer[0]) free(buffer[0]);
             free(buffer);
        }
        return n_values_per_sample; // Usually 1 for strings
    }

    // --- Modifiers ---

    void VCFRecord::set_pos(hts_pos_t pos) {
        if (is_valid_unsafe()) {
            _b->pos = pos;
        }
    }

    void VCFRecord::set_ref(const std::string& ref) {
        // Setting REF requires setting ALT simultaneously using bcf_update_alleles
        // This is a simplified placeholder; use set_alt which handles both.
        // If only REF needs changing without changing ALT count, it's complex.
        throw std::logic_error("set_ref requires setting ALT alleles simultaneously. Use set_alt.");
    }

    // Add hdr parameter
    void VCFRecord::set_alt(const VCFHeader& hdr, const std::vector<std::string>& alt_alleles) {
        if (!is_valid_unsafe() || !hdr.is_valid()) return; // Check header validity too

        // Need the REF allele first
        std::string current_ref = ref();
        if (current_ref.empty()) {
             throw std::runtime_error("Cannot set ALT alleles: current REF allele is missing or invalid.");
        }

        // Create the combined allele list (REF + ALTs)
        std::vector<const char*> alleles_c;
        alleles_c.push_back(current_ref.c_str());
        for (const auto& alt : alt_alleles) {
            alleles_c.push_back(alt.c_str());
        }

        int n_alleles = alleles_c.size();
        // Correct signature: bcf_update_alleles(const bcf_hdr_t *hdr, bcf1_t *line, const char **alleles, int nals)
        int ret = bcf_update_alleles(hdr.hts_header(), _b.get(), alleles_c.data(), n_alleles);
        if (ret != 0) {
            throw std::runtime_error("Failed to update alleles in VCF record. Error code: " + std::to_string(ret));
        }
    }

    void VCFRecord::set_qual(float qual) {
        if (is_valid_unsafe()) {
            _b->qual = (qual == FLOAT_MISSING) ? bcf_float_missing : qual;
        }
    }

    // Add hdr parameter and revert to bcf_update_id (deprecated)
    void VCFRecord::set_id(const VCFHeader& hdr, const std::string& id_str) {
        if (!is_valid_unsafe() || !hdr.is_valid()) return;
        // Signature: bcf_update_id(const bcf_hdr_t *hdr, bcf1_t *line, const char *id)
        const char* id_to_set = (id_str == STRING_MISSING || id_str.empty()) ? nullptr : id_str.c_str();
        int ret = bcf_update_id(hdr.hts_header(), _b.get(), id_to_set);
         if (ret != 0) {
            throw std::runtime_error("Failed to update ID in VCF record. Error code: " + std::to_string(ret));
        }
    }

    // Revert to bcf_add_filter (deprecated)
    int VCFRecord::add_filter(const VCFHeader& hdr, const std::string& filter_name) {
        if (!is_valid_unsafe() || !hdr.is_valid()) return -1;
        int filter_id = bcf_hdr_id2int(hdr.hts_header(), BCF_DT_ID, filter_name.c_str());
        if (filter_id < 0) return -1; // Filter not defined in header
        // Signature: bcf_add_filter(const bcf_hdr_t *hdr, bcf1_t *line, int filter_id)
        // Requires BCF_UN_FLT
        return bcf_add_filter(hdr.hts_header(), _b.get(), filter_id);
    }

    int VCFRecord::set_filters(const VCFHeader& hdr, const std::vector<std::string>& filter_names) {
        if (!is_valid_unsafe() || !hdr.is_valid()) return -1;
        // Requires BCF_UN_FLT
        std::vector<int32_t> filter_ids;
        filter_ids.reserve(filter_names.size());
        for(const auto& name : filter_names) {
            int id = bcf_hdr_id2int(hdr.hts_header(), BCF_DT_ID, name.c_str());
            if (id < 0) return -1; // One of the filters not defined
            filter_ids.push_back(id);
        }
        // Correct signature: bcf_update_filter(const bcf_hdr_t *hdr, bcf1_t *line, const int *filter_ids, int n_filters)
        // Requires BCF_UN_FLT
        // Use int32_t* for htslib - need to cast or copy if filter_ids is std::vector<int>
        // Safest is to copy to a temporary int32_t vector if needed, but htslib often uses int internally here.
        // Let's assume int* is acceptable based on common usage, but be wary.
        return bcf_update_filter(hdr.hts_header(), _b.get(), filter_ids.empty() ? nullptr : filter_ids.data(), filter_ids.size());
    }

    // Add hdr parameter
    int VCFRecord::clear_filters(const VCFHeader& hdr) {
        if (!is_valid_unsafe() || !hdr.is_valid()) return -1;
        // Correct signature: bcf_update_filter(const bcf_hdr_t *hdr, bcf1_t *line, const int *filter_ids, int n_filters)
        // Requires BCF_UN_FLT
        return bcf_update_filter(hdr.hts_header(), _b.get(), nullptr, 0);
    }

    int VCFRecord::update_info_int(const VCFHeader& hdr, const std::string& tag, const int32_t* values, int n_values) {
        if (!is_valid_unsafe() || !hdr.is_valid()) return -1;
        // Correct signature: bcf_update_info_int32(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const int32_t *values, int n)
        // Requires BCF_UN_INFO
        return bcf_update_info_int32(hdr.hts_header(), _b.get(), tag.c_str(), values, n_values);
    }

    int VCFRecord::update_info_float(const VCFHeader& hdr, const std::string& tag, const float* values, int n_values) {
        if (!is_valid_unsafe() || !hdr.is_valid()) return -1;
        // Correct signature: bcf_update_info_float(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const float *values, int n)
        // Requires BCF_UN_INFO
        return bcf_update_info_float(hdr.hts_header(), _b.get(), tag.c_str(), values, n_values);
    }

    int VCFRecord::update_info_flag(const VCFHeader& hdr, const std::string& tag, bool set) {
        if (!is_valid_unsafe() || !hdr.is_valid()) return -1;
        // Correct signature: bcf_update_info_flag(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const void *values, int n)
        // Requires BCF_UN_INFO
        // Pass NULL data pointer and set n=1 to set flag, n=0 to remove flag
        return bcf_update_info_flag(hdr.hts_header(), _b.get(), tag.c_str(), nullptr, set ? 1 : 0);
    }

    int VCFRecord::update_info_string(const VCFHeader& hdr, const std::string& tag, const char* value) {
        if (!is_valid_unsafe() || !hdr.is_valid()) return -1;
        // Correct signature: bcf_update_info_string(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const char *values)
        // Requires BCF_UN_INFO
        // Pass NULL value to remove the tag
        return bcf_update_info_string(hdr.hts_header(), _b.get(), tag.c_str(), (value && strlen(value) > 0) ? value : nullptr);
    }

    int VCFRecord::update_format_int(const VCFHeader& hdr, const std::string& tag, const int32_t* values, int n_values_per_sample) {
        if (!is_valid_unsafe() || !hdr.is_valid()) return -1;
        // Correct signature: bcf_update_format_int32(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const int32_t *values, int n)
        // Requires BCF_UN_FMT
        return bcf_update_format_int32(hdr.hts_header(), _b.get(), tag.c_str(), values, n_samples() * n_values_per_sample);
    }

    int VCFRecord::update_format_float(const VCFHeader& hdr, const std::string& tag, const float* values, int n_values_per_sample) {
        if (!is_valid_unsafe() || !hdr.is_valid()) return -1;
        // Correct signature: bcf_update_format_float(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const float *values, int n)
        // Requires BCF_UN_FMT
        return bcf_update_format_float(hdr.hts_header(), _b.get(), tag.c_str(), values, n_samples() * n_values_per_sample);
    }

    int VCFRecord::update_format_string(const VCFHeader& hdr, const std::string& tag, const char** values) {
        if (!is_valid_unsafe() || !hdr.is_valid()) return -1;
        // Correct signature: bcf_update_format_string(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const char **values, int n_sample)
        // Requires BCF_UN_FMT
        return bcf_update_format_string(hdr.hts_header(), _b.get(), tag.c_str(), values, n_samples());
    }

    int VCFRecord::update_genotypes(const VCFHeader& hdr, const int32_t* genotypes, int ploidy) {
        if (!is_valid_unsafe() || !hdr.is_valid()) return -1;
        // Correct signature: bcf_update_genotypes(const bcf_hdr_t *hdr, bcf1_t *line, const int32_t *gt_arr, int n_sample_ploidy)
        // Requires BCF_UN_FMT
        // GT tag ID is not needed for this specific function
        return bcf_update_genotypes(hdr.hts_header(), _b.get(), genotypes, n_samples() * ploidy);
    }

    // --- Output Stream ---
    std::ostream& operator<<(std::ostream& os, const VCFRecord& rec) {
        if (!rec.is_valid()) {
            os << "[Invalid VCFRecord]";
            return os;
        }
        // This provides a very basic representation.
        // For full VCF line, need header context and bcf_format
        os << "VCFRecord(rid=" << rec._b->rid
           << ", pos=" << rec.pos()
           << ", n_allele=" << rec._b->n_allele
           << ", qual=" << rec.qual()
           << ", n_sample=" << rec.n_samples()
           << ")";
        return os;
    }

} // namespace ngslib