/**
 * @file mt_caller_utils.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2025-02-24
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#ifndef _MT_CALLER_UTILS_H_
#define _MT_CALLER_UTILS_H_

#include <string>
#include <vector>
#include <htslib/bgzf.h>
#include <htslib/tbx.h>

#include "io/fasta.h"
#include "io/utils.h"

 std::string vcf_header_define(const std::string &ref_file_path, const std::vector<std::string> &samples);
 void merge_file_by_line(const std::vector<std::string> &infiles, const std::string &outfile,
                         std::string header, bool is_remove_tempfile);

 #endif
