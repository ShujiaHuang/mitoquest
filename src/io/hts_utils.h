/**
 * @brief 
 * @file hts_utils.h
 * @author Shujia Huang
 * 
 * @date 2025-03-06
 * 
 */

#ifndef _HTS_UTILS_H_
#define _HTS_UTILS_H_

#include <string>

#include <htslib/hts.h>

namespace hts {
    enum class Format {
        unknown = htsExactFormat::unknown_format,
        binary_format = htsExactFormat::binary_format,
        text_format   = htsExactFormat::text_format,
        sam  = htsExactFormat::sam,
        bam  = htsExactFormat::bam,
        bai  = htsExactFormat::bai,
        cram = htsExactFormat::cram,
        crai = htsExactFormat::crai,
        vcf  = htsExactFormat::vcf,
        bcf  = htsExactFormat::bcf,
        csi  = htsExactFormat::csi,
        gzi  = htsExactFormat::gzi,
        tbi  = htsExactFormat::tbi,
        bed  = htsExactFormat::bed
    };

    // Helper function to get format from filename
    inline Format get_format(const std::string& filename) {
        const char* fname = filename.c_str();
        htsFile *fp = hts_open(fname, "r");
        if (!fp) {
            return Format::unknown;
        }

        const htsFormat *fmt = hts_get_format(fp);
        if (!fmt) {
            return Format::unknown;
        }
        Format format = static_cast<Format>(fmt->format);
        hts_close(fp);
        
        return format;
    }

    // Helper function to check if file is CRAM
    inline bool is_cram(const std::string& filename) {
        return get_format(filename) == Format::cram;
    }
}

#endif

