/**
 * @file mt_variant_caller.h
 * @brief human mtDNA variant caller. 
 * @author Shujia Huang
 * @date 2025-02-12
 * 
 */
#ifndef _MT_VARIANT_CALLER_H_
#define _MT_VARIANT_CALLER_H_

#include <getopt.h>
#include <string>
#include <vector>
#include <set>
#include <ctime>  // clock, time_t

#include "version.h"
#include "io/fasta.h"
#include "io/bam.h"
#include "io/iobgzf.h"
#include "io/utils.h"
#include "external/thread_pool.h"

#include "basetype.h"
#include "mt_utils.h"

static const bool IS_DELETE_CACHE = true;

class MtVariantCaller {
public:
    struct Config {
        std::string reference_file;         // input reference fasta file
        std::vector<std::string> bam_files; // input BAM files
        std::string calling_regions;        // input calling regions
        std::string output_file;            // output VCF file

        int min_mapq;                       // a mapping quality score less than this value will be filtered
        int min_baseq;                      // a base quality score less than this value will be filtered
        float heteroplasmy_threshold;
        int thread_count;
        int chunk_size;               // Process this many bases per thread
        bool pairs_map_only;          // only use the paired reads which mapped to the same chromosome
        bool proper_pairs_only;       // only use properly paired reads
        bool filename_has_samplename; // use filename as sample name
    };
    ngslib::Fasta reference;          // reference fasta object

    explicit MtVariantCaller(int argc, char* argv[]);
    ~MtVariantCaller() = default; /*析构函数的实现，不需要特殊操作，放空*/

    // Main processing function
    void usage(const Config &config);
    void run() { _caller_process(); }

private:
    // Prevent copying (C++11 style)
    MtVariantCaller(const MtVariantCaller&) = delete;
    MtVariantCaller& operator=(const MtVariantCaller&) = delete;

    std::string _cmdline_string;
    Config _config;                                // command line options
    std::vector<std::string> _samples_id;          // sample ID of all alignment files (BAM/CRAM/SAM)
    std::vector<ngslib::GenomeRegion> _calling_intervals;  // vector of calling regions

    // Helper methods
    void _get_calling_interval();                  // load the calling region from input
    void _print_calling_interval();
    void _get_sample_id_from_bam();
    ngslib::GenomeRegion _make_genome_region(std::string gregion);

    void _caller_process();  // main process function
    bool _fetch_base_in_region(const ngslib::GenomeRegion genome_region, std::vector<PosVariantMap> &samples_pileup_v);

    // integrate the variant information of all samples in the region
    bool _variant_discovery(const std::vector<PosVariantMap> &samples_pileup_v, 
                            const ngslib::GenomeRegion genome_region,
                            const std::string out_vcf_fn);
};

PosVariantMap call_pileup_in_sample(const std::string sample_bam_fn, 
                                    const std::string &fa_seq,
                                    const ngslib::GenomeRegion gr,
                                    const MtVariantCaller::Config &config);

void seek_position(const std::string &fa_seq,   // must be the whole chromosome sequence
                   const std::vector<ngslib::BamRecord> &sample_map_reads,
                   const ngslib::GenomeRegion gr,
                   const int min_baseq,
                   const double min_af,
                   PosMap &sample_posinfo_map);

VariantInfo basetype_caller_unit(const AlignInfo &pos_align_info, const double min_af);

/**
 * @brief Get the Pileup object
 * 
 * @param bt BaseType
 * @param smp_bi BaseType::BatchInfo 
 * @return VariantInfo 
 * 
 */
VariantInfo get_pos_pileup(const BaseType &bt, const BaseType::BatchInfo *smp_bi);
VCFRecord call_variant_in_pos(std::vector<VariantInfo> variant_infos, const double hf_cutoff);

#endif // _MT_VARIANT_CALLER_H_

