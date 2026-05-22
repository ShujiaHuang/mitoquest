/**
 * @file mt_copynum.h
 * @brief mtDNA copy-number estimation from aligned BAM/CRAM files.
 *
 * Estimates per-chromosome relative copy number by comparing each
 * chromosome's length-normalized fragment ratio to the average of
 * the autosomes. The mitochondrial chromosome is weighted by 2 so
 * that the reported value matches the conventional "copies per
 * diploid cell" interpretation.
 *
 * @author Shujia Huang (hshujia@qq.com)
 * @date 2025-01-02
 */
#ifndef _MT_COPYNUM_H_
#define _MT_COPYNUM_H_

#include <getopt.h>
#include <string>
#include <vector>

#include "external/robin_hood.h"  // robin_hood::unordered_set

#include "version.h"
#include "io/bam.h"
#include "io/fasta.h"
#include "io/utils.h"      // ngslib::GenomeRegion

#include "algorithm.h"     // mean, standard_error
#include "mt_utils.h"      // SeqType, is_autosomal, is_mitochondrial

class MtCopyNumber {
public:
    struct Config {
        std::string reference_file; // reference genome (required for CRAM, used for GC content)
        std::string input_bam;      // input BAM/CRAM file
        std::string output_file;    // output TSV file (empty => stdout)

        int min_mapq;       // a mapping quality score less than this value will be filtered
        int thread_count;   // number of worker threads
        SeqType seq_type;   // sequencing type: AUTO|PE|SE
    };

    // Mean + 95% confidence interval of a vector of values.
    struct Statistics {
        double mean;
        double ci_lower;
        double ci_upper;

        Statistics() : mean(0.0), ci_lower(0.0), ci_upper(0.0) {}
    };

    // Per-chromosome bookkeeping.
    struct ChromosomeData {
        std::string name;
        uint32_t    length;
        int64_t     count;             // counted fragments
        double      gc_content;        // GC fraction over the whole chromosome
        double      normalized_ratio;  // (count / total_count) / (length / total_length)
        Statistics  cn_stats;          // copy number relative to autosomes (mean + CI95)
        // Track read names for paired-end de-duplication of fragments.
        robin_hood::unordered_set<std::string> processed_reads;

        ChromosomeData(const std::string &n, uint32_t l)
            : name(n), length(l), count(0), gc_content(0.0), normalized_ratio(0.0) {}
    };

    explicit MtCopyNumber(int argc, char* argv[]);
    // Direct-config constructor (mainly for tests / library callers).
    explicit MtCopyNumber(Config config);
    ~MtCopyNumber() = default;

    // Main processing function
    void run();

    // -----------------------------------------------------------------
    // Pure-math helpers, exposed as static so they can be unit-tested
    // and re-used without instantiating the class.
    // -----------------------------------------------------------------

    // GC fraction over a nucleotide string. Counts G/C / (non-N bases);
    // returns 0.0 when the sequence has no informative bases.
    static double compute_gc_content(const std::string &seq);

    // Mean + 95% CI (normal approximation, z = 1.96) over a value vector.
    // Returns all-zero Statistics for an empty input. With a single value
    // the standard error is undefined, so CI lower/upper collapse to mean.
    static Statistics compute_statistics(const std::vector<double> &values);

    // In-place: fill .normalized_ratio and .cn_stats for every contig.
    // The autosomes form the diploid baseline; mtDNA is weighted by 2 so
    // the value matches "copies per diploid cell".
    // Throws std::runtime_error if the input is empty or all-zero.
    static void compute_normalized_ratios(std::vector<ChromosomeData> &chromosomes);

private:
    // Prevent copying (C++11 style)
    MtCopyNumber(const MtCopyNumber&) = delete;
    MtCopyNumber& operator=(const MtCopyNumber&) = delete;

    std::string _cmdline_string;
    Config _config;

    // Helper methods
    static void usage();
    void _parse_args(int argc, char* argv[]);

    // Auto-detect SE/PE by inspecting the first ~1000 mapped reads.
    SeqType _detect_seq_type() const;

    // List all chromosomes from the BAM/CRAM header.
    std::vector<ChromosomeData> _load_chromosomes() const;

    // Per-chromosome fragment counter; safe to run from a worker thread
    // because each call opens its own ngslib::Bam handle.
    void _count_chrom_fragments(ChromosomeData &chrom_data, SeqType seq_type) const;

    // Compute GC content for each chromosome via ngslib::Fasta.
    void _calculate_gc_content(std::vector<ChromosomeData> &chromosomes) const;

    // Write TSV results to _config.output_file (or stdout if empty).
    void _write_results(const std::vector<ChromosomeData> &chromosomes,
                        SeqType seq_type) const;
};

#endif  // _MT_COPYNUM_H_
