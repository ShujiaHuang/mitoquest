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

    int VCFRecord::unpack(int which) const {
        if (!is_valid_unsafe()) return -1; // Or some other error code
        return bcf_unpack(_b.get(), which);
    }

    // --- Accessors for Core Fields ---
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
        unpack(BCF_UN_SHR); // Caller should ensure this if needed
        return (_b->d.id ? std::string(_b->d.id) : STRING_MISSING);
    }

    std::string VCFRecord::ref() const {
        if (!is_valid_unsafe() || _b->n_allele == 0) return "";
        unpack(BCF_UN_STR); // Alleles are usually available without unpack
        return (_b->d.allele[0] ? std::string(_b->d.allele[0]) : "");
    }

    std::vector<std::string> VCFRecord::alt() const {
        std::vector<std::string> alleles;
        if (!is_valid_unsafe() || _b->n_allele <= 1) return alleles;
        unpack(BCF_UN_STR); // Alleles usually available
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
            //  if (bcf_unpack(const_cast<bcf1_t*>(_b.get()), BCF_UN_INFO) != 0) {
             if (unpack(BCF_UN_INFO) != 0) {
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
            return -1;
        
        // 无需在此确保 FORMAT 字段是否已解包, bcf_get_genotypes 会处理
        // if (!((_b->unpacked & BCF_UN_FMT))) {
        //     if (unpack(BCF_UN_FMT) < 0) {
        //         return -1;
        //     }
        // }
        int* gt_arr = nullptr;
        int n_gt_arr = 0;
        int ret = bcf_get_genotypes(hdr.hts_header(), _b.get(), &gt_arr, &n_gt_arr);

        // 使用 RAII 确保 gt_arr 被释放
        std::unique_ptr<int, void(*)(void*)> gt_guard(gt_arr, free);
        if (ret < 0) {
            return -1; // Error or tag not found
        }

        int n_samp = n_samples();
        int ploidy = ret / n_samp;

        genotypes.resize(n_samp);
        // 对每个样本单独处理
        for (int i = 0; i < n_samp; ++i) {
            int smp_ploidy = get_sample_ploidy(hdr, i);

            // 根据实际倍性调整该样本的向量大小
            genotypes[i].reserve(smp_ploidy); 

             // 填充基因型数据
            for (int j = 0; j < smp_ploidy; ++j) {
                int32_t allele_val = gt_arr[i * ploidy + j];

                // if true, the sample has smaller ploidy
                if (allele_val == bcf_int32_vector_end) {
                    break; // 处理向量结束
                }

                genotypes[i].push_back(bcf_gt_allele(allele_val));
            }
        }
    
        return ploidy; // This is the max ploidy
    }

    int VCFRecord::get_format_int(const VCFHeader& hdr, const std::string& tag, std::vector<int32_t>& values) const {
        values.clear();
        if (!is_valid_unsafe() || !hdr.is_valid() || n_samples() == 0) return BCF_ERR_TAG_INVALID;

        int32_t* buffer = nullptr;
        int buffer_capacity = 0;  // 缓冲区容量（元素个数）
        int total_values = bcf_get_format_int32(hdr.hts_header(), _b.get(), tag.c_str(), &buffer, &buffer_capacity);

        // 错误处理
        if (total_values < 0) {
            if (buffer) free(buffer);
            
            /**
             * @brief 
             * return -1;  // 头文件中没有这个 FORMAT 字段
             * return -2;  // 类型不匹配
             * return -3;  // 该记录中不存在这个标签，或标签被标记为删除
             * return -4;  // 内存分配失败
             */
            return total_values; // Return htslib error code
        }

        // 成功：复制数据到 vector
        if (total_values > 0 && buffer) {
            values.assign(buffer, buffer + total_values);
        }
        if (buffer) free(buffer);

        int n_values_per_sample = (n_samples() > 0) ? (total_values / n_samples()) : 0;
        return n_values_per_sample; // Return number of values *per-sample*
    }

    int VCFRecord::get_format_float(const VCFHeader& hdr, const std::string& tag, std::vector<float>& values) const {
        values.clear();
        if (!is_valid_unsafe() || !hdr.is_valid() || n_samples() == 0) {
            return BCF_ERR_TAG_INVALID;
        }

        float* buffer = nullptr;
        int buffer_capacity = 0;  // 缓冲区容量（元素个数）
        int total_values = bcf_get_format_float(hdr.hts_header(), _b.get(), tag.c_str(), &buffer, &buffer_capacity);

        // 错误处理
        if (total_values < 0) {
            if (buffer) free(buffer);
            return total_values;  // 返回错误码
        }

        // 成功：复制数据到 vector
        if (total_values > 0 && buffer) {
            values.assign(buffer, buffer + total_values);
        }
        if (buffer) free(buffer);
        
        // 返回每个样本的值数量
        int n_values_per_sample = (n_samples() > 0) ? (total_values / n_samples()) : 0;
        return n_values_per_sample;
    }

    /**
     * @brief 
     * 
     * @param hdr 
     * @param tag 
     * @param values 
     * @return int 
     * 
     * Test code:
     * 
        std::vector<std::string> filters;
        int ret = record.get_format_string(header, "FT", filters);

        if (ret > 0) {
            std::cout << "Chars per sample: " << ret << std::endl;
            for (size_t i = 0; i < filters.size(); i++) {
                std::cout << "Sample " << i << ": [" << filters[i] << "]" << std::endl;
            }
        }
            
    // VCF 数据：
    // Sample0: "PASS" (4 chars + '\0')
    // Sample1: "FAIL" (4 chars + '\0')  
    // Sample2: "."    (1 char  + '\0')

    // bcf_get_format_string 返回：
    total_bytes = 15  // 3 个样本 × 5 个字符
    n_chars_per_sample = 15 / 3 = 5

    // buffer[0] 指向的数据：
    // 索引: 0    1   2   3   4   5   6   7   8   9    10  11  12  13  14
    //      'P' 'A' 'S' 'S' '\0' 'F' 'A' 'I' 'L' '\0' '.' '\0' ?   ?   ?
    //      └───── 样本0 ─────┘  └───── 样本1 ─────┘  └─ 样本2 ─┘
     *
     */
    int VCFRecord::get_format_string(const VCFHeader& hdr, const std::string& tag, std::vector<std::string>& values) const {
        values.clear();
        if (!is_valid_unsafe() || !hdr.is_valid() || n_samples() == 0) {
            return BCF_ERR_TAG_INVALID;
        }

        char** buffer = nullptr;
        int buffer_capacity = 0;
        int total_bytes = bcf_get_format_string(hdr.hts_header(), _b.get(), tag.c_str(), &buffer, &buffer_capacity);

        // 错误处理
        if (total_bytes < 0) {
            if (buffer) {
                if (*buffer) free(*buffer);  // 释放 buffer[0]
                free(buffer);
            }
            return total_bytes;
        }

        // 解析字符串数据
        int n_samp = n_samples();
        int n_chars_per_sample = 0;
        
        if (total_bytes > 0 && buffer && *buffer) {
            n_chars_per_sample = total_bytes / n_samp;
            values.reserve(n_samp);  // 预分配空间
            
            char* ptr = *buffer;  // 指向实际数据的起始位置
            for (int i = 0; i < n_samp; i++) {
                // 找到实际字符串长度（遇到 '\0' 或达到最大长度）
                int len = 0;
                while (len < n_chars_per_sample && ptr[len] != '\0') {
                    len++;
                }
                
                // 处理缺失值
                if (len == 1 && ptr[0] == '.') {
                    values.emplace_back(STRING_MISSING); // Or empty string: ""
                } else {
                    values.emplace_back(ptr, len);
                }
                
                ptr += n_chars_per_sample;  // 移动到下一个样本
            }
        }
        
        // 释放内存
        if (buffer) {
            if (*buffer) free(*buffer);  // 先释放实际数据
            free(buffer);                // 再释放指针
        }
        
        return n_chars_per_sample;
    }


    int VCFRecord::get_format_idx(const VCFHeader& hdr, const std::string& tag) const {
        if (!is_valid_unsafe() || !hdr.is_valid()) {
            return -1;
        }
    
        // 确保 FORMAT 字段已解包
        if (!(_b->unpacked & BCF_UN_FMT)) {
            if (unpack(BCF_UN_FMT) < 0) {
                return -1;
            }
        }
    
        // 获取 tag 字段在 FORMAT 中的 ID. 这里 tag 是一个字符串，可能是 "GT" 或其他格式字段
        int fmt_id = bcf_hdr_id2int(hdr.hts_header(), BCF_DT_ID, tag.data());
        if (fmt_id < 0) {
            return -1;  // 字段未定义
        }
    
        // 在记录的 FORMAT 字段中查找 tag
        for (int i = 0; i < _b->n_fmt; ++i) {
            if (_b->d.fmt[i].id == fmt_id) {
                return i;  // 返回 tag 在 FORMAT 中的索引位置
            }
        }
    
        return -1;  // tag 字段不存在于当前记录中
    }

    int VCFRecord::get_sample_ploidy(const VCFHeader& hdr, int sample_idx) const {
        // 安全性检查
        if (!is_valid_unsafe() || !hdr.is_valid() || 
            sample_idx < 0 || sample_idx >= n_samples()) {
            return -1;
        }
    
        // 获取 GT 字段索引
        int gt_idx = get_format_idx(hdr, "GT");
        if (gt_idx < 0) return -1;
    
        // 获取 fmt 结构
        bcf_fmt_t *fmt_gt = &_b->d.fmt[gt_idx];
        if (!fmt_gt) return -1;
    
        // 计算实际倍性（找到第一个 vector_end 或遍历完所有位置）
        int i, j, actual_ploidy = 0;
        #define BRANCH(type_t, convert, vector_end) { \
            /* 获取当前样本的基因型数据 */ \
            uint8_t *ptr = fmt_gt->p + sample_idx*fmt_gt->size; \
            for (j=0; j < fmt_gt->n; j++) { \
                if (convert(&ptr[j * sizeof(type_t)]) == vector_end) break; \
            } \
            actual_ploidy = j; \
        }
        switch (fmt_gt->type) {
            case BCF_BT_INT8:  BRANCH(int8_t,  le_to_i8,  bcf_int8_vector_end); break;
            case BCF_BT_INT16: BRANCH(int16_t, le_to_i16, bcf_int16_vector_end); break;
            case BCF_BT_INT32: BRANCH(int32_t, le_to_i32, bcf_int32_vector_end); break;
            default:
                throw std::runtime_error("Unexpected case: " + std::to_string(fmt_gt->type) + "\n"); 
        }

        return actual_ploidy;
    }

    std::vector<int> VCFRecord::get_sample_ploidies(const VCFHeader& hdr) const {
        std::vector<int> ploidies;
        if (!is_valid_unsafe() || !hdr.is_valid()) {
            return ploidies;
        }
    
        int n_samp = n_samples();
        ploidies.reserve(n_samp);
    
        for (int i = 0; i < n_samp; ++i) {
            ploidies.push_back(get_sample_ploidy(hdr, i));
        }
    
        return ploidies;
    }

    int VCFRecord::get_max_ploidy(const VCFHeader& hdr) const {
        if (!is_valid_unsafe() || !hdr.is_valid()) {
            return -1;
        }
    
        int max_ploidy = 0;
        int n_samp = n_samples();
    
        for (int i = 0; i < n_samp; ++i) {
            int ploidy = get_sample_ploidy(hdr, i);
            max_ploidy = std::max(max_ploidy, ploidy);
        }
    
        return max_ploidy;
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

    bool VCFRecord::cleanup_genotypes(const ngslib::VCFHeader& hdr) {
        // 1. 基本检查
        if (!is_valid_unsafe() || !hdr.is_valid()) return false;
    
        // 2. 确保记录已解包
        if (!(_b->unpacked & BCF_UN_STR)) {
            if (unpack(BCF_UN_STR) < 0) {
                return false;
            }
        }
    
        // 3. 获取当前的等位基因信息
        std::string current_ref = ref();
        std::vector<std::string> current_alts = alt();
        if (current_ref.empty()) return false;
    
        // 4. 获取基因型数据
        std::vector<std::vector<int>> genotypes;
        int max_ploidy = get_genotypes(hdr, genotypes);
        if (max_ploidy <= 0) return false;
    
        // 5. 标记使用的等位基因
        std::vector<bool> allele_used(1 + current_alts.size(), false);
        allele_used[0] = true;  // REF 总是保留
        for (const auto& sample_gt : genotypes) {
            for (int gt : sample_gt) {
                if (gt >= 0 && gt < allele_used.size()) {
                    allele_used[gt] = true;
                }
            }
        }
    
        // 6. 检查是否需要清理
        bool all_alleles_used = std::all_of(allele_used.begin(), allele_used.end(), [](bool v){ return v; });
        bool all_single_base = (current_ref.length() == 1) && std::all_of(
            current_alts.begin(), current_alts.end(), [](const std::string& alt)
            { return alt.length() == 1; }
        );
                              
        if (all_alleles_used && all_single_base) {
            return true;  // 所有等位基因都在使用且都是单碱基，无需清理
        }
    
        // 7. 创建新的等位基因列表
        std::vector<std::string> new_alts;
        std::vector<int> allele_map(allele_used.size(), -1);
        allele_map[0] = 0;  // REF 映射到自身
        int new_index = 1;
        
        for (size_t i = 0; i < current_alts.size(); ++i) {
            if (allele_used[i + 1]) {
                new_alts.push_back(current_alts[i]);
                allele_map[i + 1] = new_index++;
            }
        }
    
        // 8. 更新基因型数据
        std::vector<std::vector<int>> new_genotypes;
        new_genotypes.reserve(genotypes.size());
        for (const auto& sample_gt : genotypes) {
            std::vector<int> new_gt;
            new_gt.reserve(sample_gt.size());
            for (int gt : sample_gt) {  

                if (allele_used[gt] && gt >= 0 &&  // gt may be missing
                    gt < allele_used.size()) 
                {
                    new_gt.push_back(allele_map[gt]);  // 映射到新的等位基因索引
                } else {
                    new_gt.push_back(-1);  // 标记缺失值
                }
            }
            new_genotypes.push_back(std::move(new_gt));
        }
    
        // 9. 更新记录
        // 先更新等位基因列表
        if (!new_alts.empty() && update_alleles(hdr, current_ref, new_alts) < 0) {
            return false;
        }
    
        // 然后更新基因型数据
        if (update_genotypes(hdr, new_genotypes) < 0) {
            return false;
        }
    
        return true;
    }

    int VCFRecord::update_alleles(const VCFHeader& hdr, const std::string& ref, 
                                  const std::vector<std::string>& alts) {
        // 1. 安全性检查
        if (!is_valid_unsafe() || !hdr.is_valid() || alts.empty()) return -1;
        
        // 2. 确保记录已解包
        if (!(_b->unpacked & BCF_UN_STR)) {
            if (unpack(BCF_UN_STR) < 0) {
                return -1;
            }
        }
    
        // 3. 简化 REF 和 ALT 序列
        std::string simplified_ref = ref;
        std::vector<std::string> simplified_alts = alts;
        
        // 3.1 找到所有序列共同的前缀长度
        size_t prefix_len = 0;
        while (prefix_len < ref.length()) {
            char c = ref[prefix_len];
            bool has_prefix = true;
            for (const auto& alt : alts) {
                // 只在相同长度的序列中找共同前缀
                if ((alt.length() != ref.length()) || (prefix_len >= alt.length()) || (alt[prefix_len] != c)) {
                    has_prefix = false;
                    break;
                }
            }

            if (has_prefix) { // 如果有共同前缀，增加长度
                prefix_len++;
            } else { // 否则，退出循环
                break;
            }
        }
        
        // 3.2 找到所有序列共同的后缀长度
        size_t suffix_len = 0;
        while (prefix_len + suffix_len < ref.length()) {
            size_t pos = ref.length() - 1 - suffix_len;
            char c = ref[pos];
            bool has_suffix = true;
            for (const auto& alt : alts) {
                if ((prefix_len + suffix_len >= alt.length()) || (alt[alt.length() - 1 - suffix_len] != c)) {
                    has_suffix = false;
                    break;
                }
            }

            if (has_suffix) {
                suffix_len++;
            } else {
                break;
            }
        }
    
        // 3.3 如果有共同前缀或后缀，删除它们（前缀要保留至少一个碱基）
        if (prefix_len > 0 || suffix_len > 0) {
            prefix_len = prefix_len > 1 ? prefix_len - 1 : 0;  // 前缀要保留至少一个碱基
            size_t new_len = ref.length() - prefix_len - suffix_len;
            simplified_ref = ref.substr(prefix_len, new_len);
            
            for (size_t i = 0; i < alts.size(); ++i) {
                new_len = alts[i].length() - prefix_len - suffix_len;
                simplified_alts[i] = alts[i].substr(prefix_len, new_len);
            }
        }
    
        // 4. 准备等位基因数组
        const int n_alleles = 1 + simplified_alts.size();
        std::vector<const char*> alleles(n_alleles);
        
        // 5. 设置 REF 和 ALT 等位基因
        alleles[0] = simplified_ref.c_str();
        for (size_t i = 0; i < simplified_alts.size(); ++i) {
            alleles[i + 1] = simplified_alts[i].c_str();
        }
    
        // 6. 调用 htslib 函数更新等位基因
        int ret = bcf_update_alleles(
            hdr.hts_header(),
            _b.get(),
            alleles.data(),
            n_alleles
        );

        // 7. 如果前缀被删除，需要更新位置
        // 此时的 prefix_len 已在 3.3 代码中被更为原始 prefix_len - 1，所以和 0 对比即可
        if (prefix_len > 0) {
            _b->pos += prefix_len;  // 更新变异位置
        }
    
        return (ret >= 0) ? 0 : -1;
    }

    int VCFRecord::update_genotypes(const VCFHeader& hdr, const std::vector<std::vector<int>>& genotypes) {
        // 1. 安全性检查
        if (!is_valid_unsafe() || !hdr.is_valid()) return -1;
        
        // 2. 确保 FORMAT 字段已解包
        if (!((_b->unpacked & BCF_UN_FMT))) {
            if (unpack(BCF_UN_FMT) < 0) {
                return -1;
            }
        }
        
        // 3. 找到 GT 字段
        int gt_idx = get_format_idx(hdr, "GT");
        if (gt_idx < 0) return -1;  // GT 未在头文件中定义
        
        bcf_fmt_t *fmt_gt = &_b->d.fmt[gt_idx];
        int original_ploidy = fmt_gt->n; // 记录原始大小
        
        // 4. 计算需要的总空间
        int n_samp = n_samples();
        if (genotypes.size() != n_samp) return -1; // n_samples 与 genotypes.size() 不匹配
        
        // 5. 找到最大倍性
        int max_ploidy = 0;
        for (const auto& sample_gt : genotypes) {
            max_ploidy = std::max(max_ploidy, static_cast<int>(sample_gt.size()));
        }
        if (max_ploidy <= 0) return -1;
        
        // 按照最大倍性为样本分配空间
        std::vector<int32_t> gt_arr(n_samp * max_ploidy);
        
        // 6. 为每个样本设置基因型
        for (int i = 0; i < n_samp; ++i) {
            const auto& sample_gt = genotypes[i];
            int32_t* curr_sample = gt_arr.data() + i * max_ploidy;
            
            // 设置该样本的基因型
            for (size_t j = 0; j < max_ploidy; ++j) {
                if (j < sample_gt.size()) {
                    // 实际的基因型值
                    if (sample_gt[j] < 0) {
                        curr_sample[j] = bcf_gt_missing;  // 标记缺失的基因型，输出为 ‘.’
                    } else {
                        curr_sample[j] = bcf_gt_is_phased(sample_gt[j]) ? bcf_gt_phased(sample_gt[j]) : bcf_gt_unphased(sample_gt[j]);
                    }
                } else {
                    // 对于倍性较低的样本，用向量结束标记填充
                    curr_sample[j] = bcf_int32_vector_end;
                }
            }
        }

        // 7. 更新记录中的基因型数据
        int ret = bcf_update_genotypes(hdr.hts_header(), _b.get(), gt_arr.data(), n_samp * max_ploidy);

        // 检查最大倍性的大小是否发生了改变
        if (ret >= 0 && _b->n_fmt > 0) {
            bcf_fmt_t *fmt = &_b->d.fmt[gt_idx];
            if (get_max_ploidy(hdr) != original_ploidy) {
                std::cerr << "[INFO]: Ploidy changed from " << original_ploidy << " to " 
                          << get_max_ploidy(hdr) << " at " << chrom(hdr) << ":" << pos() + 1 
                          << "\n";
            }
        }
        return (ret >= 0) ? max_ploidy : ret;
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
           << ", pos="         << rec.pos() + 1
           << ", n_allele="    << rec._b->n_allele
           << ", qual="        << rec.qual()
           << ", n_sample="    << rec.n_samples()
           << ")";
        return os;
    }

} // namespace ngslib