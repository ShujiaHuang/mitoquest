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

#include "external/robin_hood.h"  // robin_hood::unordered_map

struct GenomeRegion {
    std::string ref_id;
    uint32_t start;
    uint32_t end;

    GenomeRegion() : ref_id(""), start(0), end(0) {};
    GenomeRegion(const std::string& rid, uint32_t s, uint32_t e) : ref_id(rid), start(s), end(e) {
        if (start > end) {
            throw std::invalid_argument("[ERROR] start postion is larger than end position in "
                                        "GenomeRegion: " + ref_id + ":" + 
                                        std::to_string(start) + "-" + std::to_string(end));
        }
    };
};

struct AlignBase {
    std::string base; // read base (could be a base or indels sequence)
    int base_qual;    // base quality, get mean quality of seq if it is indels
    int rp;           // base position in read, the first base is 1, record the leftmost position for indels 

    int mapq;         // mapping quality
    char map_strand;  // mapping reference strand, should be one of '-', '+' or '*'
};

struct AlignInfo {
    std::string ref_id;
    uint32_t    ref_pos;
    std::string ref_base; // reference base
    std::vector<AlignBase> read_bases;  // mapped read bases

    AlignInfo() : ref_id(""), ref_pos(0), ref_base("") {};
    AlignInfo(const std::string& rid, uint32_t pos, const std::string& rb) 
        : ref_id(rid), ref_pos(pos), ref_base(rb) {};

};

typedef robin_hood::unordered_map<std::string, AlignInfo> Ref2BaseMap;
typedef robin_hood::unordered_map<uint32_t, Ref2BaseMap> PosMap;  // give a short name for this type
typedef std::vector<PosMap> PosMapVector;  // record the alignment information for each position in a region

 std::string vcf_header_define(const std::string &ref_file_path, const std::vector<std::string> &samples);
 void merge_file_by_line(const std::vector<std::string> &infiles, const std::string &outfile,
                         std::string header, bool is_remove_tempfile);

 #endif
