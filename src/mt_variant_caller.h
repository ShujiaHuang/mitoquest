/**
 * @file mt_variant_caller.h
 * @brief human mtDNA variant caller. 
 * @author Shujia Huang
 * @date 2025-02-12
 * 
 */
#ifndef _MT_VARIANT_CALLER_H_
#define _MT_VARIANT_CALLER_H_

#include <string>
#include <vector>
#include <ctime>  // clock, time_t

#include "external/thread_pool.h"

#include "io/fasta.h"
#include "io/bam.h"
#include "io/utils.h"
#include "mt_caller_utils.h"

static const bool IS_DELETE_CACHE = true;

class MtVariantCaller {
public:
    struct Config {
        std::string reference_file;
        std::vector<std::string> bam_files;  // input multiple BAM files
        std::string calling_regions;         // calling regions
        std::string output_file;             // output VCF file

        int min_mapq  = 0;   // a mapping quality score less than this value will be filtered
        int min_baseq = 20;  // a base quality score less than this value will be filtered
        float heteroplasmy_threshold = 0.2;
        int thread_count = 1;
        int chunk_size   = 1000;              // Process this many bases per thread
        bool filename_has_samplename = false; // use filename as sample name
    };

    ngslib::Fasta reference;  // reference fasta object

    explicit MtVariantCaller(const Config& config);
    ~MtVariantCaller();

    // Main processing function
    void print_calling_interval();
    bool run();

private:
    // Prevent copying (C++11 style)
    MtVariantCaller(const MtVariantCaller&) = delete;
    MtVariantCaller& operator=(const MtVariantCaller&) = delete;

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
        std::string base;
        int base_qual;    // base quality
        int rp;           // base position in read
        int mapq;         // mapping quality
        char map_strand;  // mapping reference strand, should be one of '-', '+' or '.'

        // 注意这个赋值顺序要和结构体定义的顺序一致
        AlignBase() : base(""), base_qual(0), rp(0), mapq(0), map_strand('.') {};
        AlignBase(const std::string& b, int bq, int rp, int mq, char ms) 
            : base(b), base_qual(bq), rp(rp), mapq(mq), map_strand(ms) {}
    };

    struct Variant {
        std::string ref_id;
        uint32_t ref_pos;
        std::string ref_base;   // reference base
        std::string alt_base;   // alternative base

        float frequency;        // frequency of non-reference base
        int quality;
        int depth;
    };

    // Member variables
    Config _config;
    std::vector<std::string> _samples_id;          // sample ID of all alignment files (BAM/CRAM/SAM)
    std::vector<GenomeRegion> _calling_intervals;  // vector of calling regions

    // Helper methods
    void _get_calling_interval();  // load the calling region from input
    void _get_sample_id_from_bam();
    GenomeRegion _make_genome_region(std::string gregion);
    bool _caller_process();  // main process function

    // bool _pileup(const GenomeRegion& region);
    // void _pileup_worker(const GenomeRegion& region, const ngslib::Bam& bam, const std::string& sample_id);

};

#endif // _MT_VARIANT_CALLER_H_
