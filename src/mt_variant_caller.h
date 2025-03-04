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

#include "basetype.h"
#include "mt_caller_utils.h"

#include "io/fasta.h"
#include "io/bam.h"
#include "external/thread_pool.h"

static const bool IS_DELETE_CACHE = true;



class MtVariantCaller {
public:
    struct Config {
        std::string reference_file;          // input reference fasta file
        std::vector<std::string> bam_files;  // input BAM files
        std::string calling_regions;         // input calling regions
        std::string output_file;             // output VCF file

        int min_mapq  = 0;   // a mapping quality score less than this value will be filtered
        int min_baseq = 20;  // a base quality score less than this value will be filtered
        float heteroplasmy_threshold = 0.2;
        int thread_count = 1;
        int chunk_size   = 1000;              // Process this many bases per thread
        bool pairs_map_only    = false;       // only use the paired reads which mapped to the same chromosome
        bool proper_pairs_only = false;       // only use properly paired reads
        bool filename_has_samplename = false; // use filename as sample name
    };
    ngslib::Fasta reference;  // reference fasta object

    explicit MtVariantCaller(const Config& config);
    ~MtVariantCaller();

    // Main processing function
    bool run();

private:
    // Prevent copying (C++11 style)
    MtVariantCaller(const MtVariantCaller&) = delete;
    MtVariantCaller& operator=(const MtVariantCaller&) = delete;

    // Member variables
    Config _config;                                // command line options
    std::vector<std::string> _samples_id;          // sample ID of all alignment files (BAM/CRAM/SAM)
    std::vector<GenomeRegion> _calling_intervals;  // vector of calling regions

    // Helper methods
    void _get_calling_interval();  // load the calling region from input
    void _print_calling_interval();
    void _get_sample_id_from_bam();
    GenomeRegion _make_genome_region(std::string gregion);

    bool _caller_process();  // main process function
    bool _fetch_base_in_region(const GenomeRegion genome_region,
                               PosMapVector &batchsamples_posinfomap_vector);

};

PosMap call_variant_in_sample(const std::string sample_bam_fn, 
                              const std::string &fa_seq,
                              const GenomeRegion gr,
                              const MtVariantCaller::Config &config);

void seek_position(const std::string &fa_seq,   // must be the whole chromosome sequence
                   const std::vector<ngslib::BamRecord> &sample_map_reads,
                   const GenomeRegion gr,
                   PosMap &sample_posinfo_map);

VariantInfo variant_caller_unit(const AlignInfo &pos_align_info, double min_af);

/**
 * @brief Get the variant object
 * 
 * @param bt BaseType
 * @param smp_bi BatchInfo 
 * @return VariantInfo 
 */
VariantInfo get_variant(const BaseType &bt, const BatchInfo *smp_bi);

#endif // _MT_VARIANT_CALLER_H_

