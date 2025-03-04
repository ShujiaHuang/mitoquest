/**
 * @file mt_caller_utils.h
 * @brief  Utilities for mtDNA variant caller 
 * 
 * @author Shujia Huang (hshujia@qq.com)
 * @date 2025-02-24
 * 
 */
#ifndef _MT_CALLER_UTILS_H_
#define _MT_CALLER_UTILS_H_

#include <cmath>   // use log function
#include <string>
#include <vector>
#include <htslib/bgzf.h>
#include <htslib/tbx.h>

#include "io/fasta.h"
#include "io/utils.h"
#include "external/robin_hood.h"  // robin_hood::unordered_map, robin_hood::unordered_set

#include "algorithm.h"

// define data types for variant calling
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
    std::string ref_base;
    std::string read_base; // read base (could be a base or indels sequence)
    char base_qual;   // read base quality, get mean quality of seq if it is indels
    int rpr;          // read position rank, the first base is 1, record the leftmost position for indels 

    int mapq;         // mapping quality
    char map_strand;  // mapping reference strand, should be one of '-', '+' or '*'
};

struct AlignInfo {
    std::string ref_id;
    uint32_t    ref_pos;
    std::vector<AlignBase> align_bases;

    AlignInfo() : ref_id(""), ref_pos(0) {};
    AlignInfo(const std::string& rid, uint32_t pos) : ref_id(rid), ref_pos(pos) {};
};

typedef robin_hood::unordered_map<uint32_t, AlignInfo> PosMap;  // key: ref_pos, value: AlignInfo
typedef std::vector<PosMap> PosMapVector; // record the alignment information for each position in a region

typedef struct {
    int ref_fwd, ref_rev;
    int alt_fwd, alt_rev;
    double fs;
    double sor;
} StrandBiasInfo;

struct VariantInfo {
    std::string ref_id;
    uint32_t ref_pos;
    int total_depth;

    std::vector<std::string> ref_bases;  // REF, could be a base or indels sequence  
    std::vector<std::string> alt_bases;  // ALT, could be a base or indels sequence
    std::vector<std::string> var_types;  // REF, SNV, INS, DEL, or MNV
    std::vector<int> depths;             // depth for each type of base
    std::vector<double> freqs;           // frequency for each type of base
    std::vector<StrandBiasInfo> strand_bias;

    VariantInfo() : ref_id(""), ref_pos(0), total_depth(0) {};
    VariantInfo(const std::string& rid, uint32_t pos) : ref_id(rid), ref_pos(pos), total_depth(0) {};
};

// get the total depth for a reference position
int get_total_depth(const AlignInfo &align_infor);

/**
 * @brief Mann-Whitney-Wilcoxon Rank Sum Test for REF and ALT array.
 * 
 * @param ref_base A reference base
 * @param alt_bases_string 
 * @param bases 
 * @param values 
 * @return double   Phred-scale value of ranksum-test pvalue
 * 
 * Note: There's some difference between scipy.stats.ranksums with R's wilcox.test:
 *       https://stackoverflow.com/questions/12797658/pythons-scipy-stats-ranksums-vs-rs-wilcox-test
 * 
 */
double ref_vs_alt_ranksumtest(const char ref_base,
                              const std::string alt_bases_string,
                              const std::vector<char> &bases,
                              const std::vector<int> &values);

double ref_vs_alt_ranksumtest(const char ref_base, 
                              const std::string alt_bases_string,
                              const std::vector<char> &bases,
                              const std::vector<char> &values);

StrandBiasInfo strand_bias(const char ref_base, 
                           const std::string alt_bases_string,
                           const std::vector<char> &bases,
                           const std::vector<char> &strands);

std::string vcf_header_define(const std::string &ref_file_path, const std::vector<std::string> &samples);
void merge_file_by_line(const std::vector<std::string> &infiles, const std::string &outfile,
                        std::string header, bool is_remove_tempfile);

 #endif
