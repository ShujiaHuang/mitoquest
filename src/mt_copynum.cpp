/**
 * @file mt_copynum.cpp
 * @brief Implementation of the 'copynum' subcommand of mitoquest.
 *
 * @author Shujia Huang (hshujia@qq.com)
 * @date 2025-01-02
 */
#include "mt_copynum.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <thread>

#include "external/thread_pool.h"  // ThreadPool

// Print usage information for the 'copynum' command.
void MtCopyNumber::usage() {
    std::cerr << "Usage: mitoquest copynum [options] <input.bam/cram>\n\n"
                 "Description:\n"
                 "  Estimate per-chromosome relative copy number from a sorted/indexed BAM\n"
                 "  or CRAM file. The autosomal chromosomes are used as the diploid baseline\n"
                 "  (copy number = 2); the mitochondrial chromosome is reported in the same\n"
                 "  scale, i.e. the expected number of mtDNA molecules per diploid cell.\n"
                 "\nOptions:\n"
                 "  -r, --reference FILE   Reference genome FASTA file (required; needed for\n"
                 "                         CRAM decoding and GC content calculation).\n"
                 "  -o, --output    FILE   Output TSV file (default: stdout).\n"
                 "  -q, --mapq      INT    Minimum mapping quality score [0].\n"
                 "  -t, --threads   INT    Number of worker threads [hardware_concurrency].\n"
                 "  -s, --seqtype   STR    Sequencing type: auto|pe|se [auto].\n"
                 "  -h, --help             Print this help message.\n\n"
              << "Version: " << MITOQUEST_VERSION << "\n"
              << std::endl;
}

// Parse command line arguments specific to the 'copynum' command.
void MtCopyNumber::_parse_args(int argc, char* argv[]) {
    // Defaults
    _config.reference_file.clear();
    _config.input_bam.clear();
    _config.output_file.clear();
    _config.min_mapq     = 0;
    _config.thread_count = std::thread::hardware_concurrency();
    if (_config.thread_count <= 0) _config.thread_count = 2;
    _config.seq_type     = SeqType::AUTO;

    // Save the complete command line for the output header.
    _cmdline_string = "#mitoquest_copynum_command=";
    for (int i = 0; i < argc; ++i) {
        _cmdline_string += (i > 0) ? " " + std::string(argv[i]) : std::string(argv[i]);
    }

    static const struct option long_options[] = {
        {"reference", required_argument, 0, 'r'},
        {"output",    required_argument, 0, 'o'},
        {"mapq",      required_argument, 0, 'q'},
        {"threads",   required_argument, 0, 't'},
        {"seqtype",   required_argument, 0, 's'},
        {"help",      no_argument,       0, 'h'},
        {0, 0, 0, 0}  // Terminator
    };

    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "r:o:q:t:s:h", long_options, &option_index)) != -1) {
        switch (c) {
            case 'r': _config.reference_file = optarg;         break;
            case 'o': _config.output_file    = optarg;         break;
            case 'q': _config.min_mapq       = std::stoi(optarg); break;
            case 't': _config.thread_count   = std::stoi(optarg); break;
            case 's': {
                std::string v = optarg;
                std::transform(v.begin(), v.end(), v.begin(), ::tolower);
                if      (v ==   "pe") _config.seq_type = SeqType::PE;
                else if (v ==   "se") _config.seq_type = SeqType::SE;
                else if (v == "auto") _config.seq_type = SeqType::AUTO;
                else {
                    std::cerr << "Error: Invalid sequencing type '" << optarg
                              << "'. Use 'auto', 'pe', or 'se'.\n";
                    usage();
                    exit(EXIT_FAILURE);
                }
                break;
            }
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
            case '?':
                // getopt_long already printed an error message.
                usage();
                exit(EXIT_FAILURE);
            default:
                abort();  // should not happen
        }
    }

    // The single positional argument is the input BAM/CRAM file.
    if (optind >= argc) {
        std::cerr << "Error: Input BAM/CRAM file is required.\n";
        usage();
        exit(EXIT_FAILURE);
    }
    _config.input_bam = argv[optind++];
    if (optind < argc) {
        std::cerr << "Error: Only one input BAM/CRAM file is allowed.\n";
        usage();
        exit(EXIT_FAILURE);
    }

    // Sanity checks
    if (_config.reference_file.empty()) {
        std::cerr << "Error: Reference genome (-r/--reference) is required.\n";
        usage();
        exit(EXIT_FAILURE);
    }
    if (!ngslib::is_readable(_config.input_bam)) {
        throw std::runtime_error("[copynum] Input file not readable: " + _config.input_bam);
    }
    if (!ngslib::is_readable(_config.reference_file)) {
        throw std::runtime_error("[copynum] Reference file not readable: " + _config.reference_file);
    }
    if (_config.thread_count <= 0) _config.thread_count = 1;
    if (_config.min_mapq < 0)      _config.min_mapq     = 0;
}

// Constructor
MtCopyNumber::MtCopyNumber(int argc, char* argv[]) {
    _parse_args(argc, argv);
}

// Direct-config constructor (does NOT validate file paths so callers can
// build a Config in tests without touching the filesystem).
MtCopyNumber::MtCopyNumber(Config config) : _config(std::move(config)) {
    if (_config.thread_count <= 0) _config.thread_count = 1;
    if (_config.min_mapq    <  0) _config.min_mapq     = 0;
    _cmdline_string = "#mitoquest_copynum_command=<programmatic>";
}

// Auto-detect sequencing type from the first ~1000 mapped reads.
SeqType MtCopyNumber::_detect_seq_type() const {
    static const int READ_LIMIT = 1000;

    ngslib::Bam bf(_config.input_bam, "r", _config.reference_file);
    ngslib::BamRecord rec;
    int n_mapped = 0;
    // No fetch() => read sequentially from the start of the file.
    while (bf.read(rec) >= 0 && n_mapped < READ_LIMIT) {
        if (!rec.is_mapped()) continue;
        if (rec.is_paired()) return SeqType::PE;
        ++n_mapped;
    }
    return SeqType::SE;  // 1000 mapped reads, none paired
}

// List all chromosomes from the BAM/CRAM header.
std::vector<MtCopyNumber::ChromosomeData> MtCopyNumber::_load_chromosomes() const {
    std::vector<ChromosomeData> chromosomes;

    ngslib::Bam bf(_config.input_bam, "r", _config.reference_file);
    ngslib::BamHeader &hdr = bf.header();
    sam_hdr_t *h = hdr.h();
    if (!h) {
        throw std::runtime_error("[copynum] Failed to read header from " + _config.input_bam);
    }
    chromosomes.reserve(h->n_targets);
    for (int i = 0; i < h->n_targets; ++i) {
        chromosomes.emplace_back(hdr.seq_name(i), hdr.seq_length(i));
    }
    return chromosomes;
}

// Count fragments for a single chromosome. Each call opens its own ngslib::Bam
// handle so the function is safe to invoke from worker threads in parallel.
void MtCopyNumber::_count_chrom_fragments(ChromosomeData &chrom_data,
                                          SeqType seq_type) const {
    ngslib::Bam bf(_config.input_bam, "r", _config.reference_file);
    ngslib::BamHeader &hdr = bf.header();

    // Skip silently if the chromosome is absent in the BAM index (defensive).
    if (hdr.name2id(chrom_data.name) < 0) {
        return;
    }
    if (!bf.fetch(chrom_data.name, 0, chrom_data.length)) {
        return;
    }

    ngslib::BamRecord rec;
    while (bf.read(rec) >= 0) {
        // Skip unmapped, secondary, and supplementary alignments.
        if (!rec.is_mapped() || rec.is_secondary() || rec.is_supplementary()) {
            continue;
        }
        if (rec.mapq() < _config.min_mapq) {
            continue;
        }

        if (seq_type == SeqType::SE) {
            // Single-end: count every primary alignment.
            ++chrom_data.count;
        } else {
            // Paired-end: count each fragment exactly once via read1.
            // For unpaired reads in a PE BAM (e.g. orphan mates), fall back
            // to qname-based de-duplication so they are not double counted.
            if (rec.is_read1()) {
                ++chrom_data.count;
            } else if (!rec.is_paired()) {
                std::string qn = rec.qname();
                if (chrom_data.processed_reads.insert(qn).second) {
                    ++chrom_data.count;
                }
            }
        }
    }

    // Free the read-name set early; we no longer need it.
    chrom_data.processed_reads.clear();
}

// GC fraction over a nucleotide string. Counts G/C among non-N bases.
double MtCopyNumber::compute_gc_content(const std::string &seq) {
    int64_t gc    = 0;
    int64_t valid = 0;
    for (char c : seq) {
        char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        if (u == 'N') continue;
        ++valid;
        if (u == 'G' || u == 'C') ++gc;
    }
    return (valid > 0) ? static_cast<double>(gc) / valid : 0.0;
}

// Mean + 95% CI (normal approximation, z = 1.96).
MtCopyNumber::Statistics MtCopyNumber::compute_statistics(const std::vector<double> &values) {
    Statistics stats;
    if (values.empty()) return stats;

    stats.mean = mean(values);
    if (values.size() < 2) {
        stats.ci_lower = stats.mean;
        stats.ci_upper = stats.mean;
        return stats;
    }

    double se      = standard_error(values);
    double ci_half = 1.96 * se;  // 95% CI under the normal approximation
    stats.ci_lower = stats.mean - ci_half;
    stats.ci_upper = stats.mean + ci_half;
    return stats;
}

// Length-normalized fragment ratio + per-chromosome copy number relative to
// the average autosomal normalized ratio.
void MtCopyNumber::compute_normalized_ratios(std::vector<ChromosomeData> &chromosomes) {
    if (chromosomes.empty()) {
        throw std::runtime_error("[copynum] No chromosomes provided.");
    }

    int64_t total_fragments     = 0;
    int64_t total_genome_length = 0;
    for (const auto &chrom : chromosomes) {
        total_fragments     += chrom.count;
        total_genome_length += chrom.length;
    }

    if (total_fragments == 0 || total_genome_length == 0) {
        throw std::runtime_error("[copynum] No fragments counted; check input BAM/CRAM and "
                                 "filtering options (e.g. --mapq).");
    }

    for (auto &chrom : chromosomes) {
        double fragment_ratio  = static_cast<double>(chrom.count)  / total_fragments;
        double length_ratio    = static_cast<double>(chrom.length) / total_genome_length;
        chrom.normalized_ratio = (length_ratio > 0.0) ? fragment_ratio / length_ratio : 0.0;
    }

    // Each chromosome's copy number is reported relative to the average
    // autosomal normalized ratio. The mitochondrial chromosome is weighted
    // by 2 so the value matches "copies per diploid cell".
    for (auto &chrom : chromosomes) {
        double copy_w = is_mitochondrial(chrom.name) ? 2.0 : 1.0;
        std::vector<double> ratios;
        ratios.reserve(chromosomes.size());
        for (const auto &chrtmp : chromosomes) {
            if (is_autosomal(chrtmp.name) && chrtmp.normalized_ratio > 0.0) {
                ratios.push_back(copy_w * chrom.normalized_ratio / chrtmp.normalized_ratio);
            }
        }
        chrom.cn_stats = compute_statistics(ratios);
    }
}

// Compute GC content for each chromosome via ngslib::Fasta.
void MtCopyNumber::_calculate_gc_content(std::vector<ChromosomeData> &chromosomes) const {
    ngslib::Fasta reference(_config.reference_file);

    for (auto &chrom : chromosomes) {
        if (!reference.has_seq(chrom.name)) {
            chrom.gc_content = 0.0;
            continue;
        }
        // Use the simple chrom-name form to fetch the full sequence.
        std::string seq = reference.fetch(chrom.name);
        chrom.gc_content = compute_gc_content(seq);
    }
}

// Write TSV results to _config.output_file (or stdout if empty).
void MtCopyNumber::_write_results(const std::vector<ChromosomeData> &chromosomes,
                                  SeqType seq_type) const {
    std::ofstream ofs;
    std::ostream *out = &std::cout;
    if (!_config.output_file.empty()) {
        ofs.open(_config.output_file);
        if (!ofs.is_open()) {
            throw std::runtime_error("[copynum] Could not open output file: " + _config.output_file);
        }
        out = &ofs;
    }

    *out << _cmdline_string << "\n";
    *out << "#Sequencing type: "
         << (seq_type == SeqType::PE ? "paired-end" : "single-end") << "\n";
    *out << "#Chromosome\tFragments\tChrom_Length\tGC_Content\t"
            "Fragment_Normalized_Ratio\tCopyNum\tCopyNum-CI95-Lower\tCopyNum-CI95-Upper\n";

    for (const auto &chrom : chromosomes) {
        *out << chrom.name   << "\t"
             << chrom.count  << "\t"
             << chrom.length << "\t"
             << std::fixed   << std::setprecision(4) << chrom.gc_content        << "\t"
             << std::fixed   << std::setprecision(4) << chrom.normalized_ratio  << "\t"
             << std::fixed   << std::setprecision(4) << chrom.cn_stats.mean     << "\t"
             << std::fixed   << std::setprecision(4) << chrom.cn_stats.ci_lower << "\t"
             << std::fixed   << std::setprecision(4) << chrom.cn_stats.ci_upper << "\n";
    }
}

// Main entry point
void MtCopyNumber::run() {
    // 1. Resolve sequencing type.
    SeqType seq_type = _config.seq_type;
    if (seq_type == SeqType::AUTO) {
        seq_type = _detect_seq_type();
        std::cout << "[INFO] Auto-detected sequencing type: "
                  << (seq_type == SeqType::PE ? "paired-end" : "single-end")
                  << std::endl;
    }

    // 2. Enumerate chromosomes from the BAM header.
    std::vector<ChromosomeData> chromosomes = _load_chromosomes();
    if (chromosomes.empty()) {
        throw std::runtime_error("[copynum] No chromosomes found in BAM/CRAM header.");
    }

    // 3. Count fragments per chromosome in parallel (one task per chromosome).
    {
        ThreadPool pool(_config.thread_count);
        std::vector<std::future<void>> results;
        results.reserve(chromosomes.size());
        for (auto &chrom : chromosomes) {
            results.emplace_back(pool.submit(std::bind(&MtCopyNumber::_count_chrom_fragments,
                                                       this,
                                                       std::ref(chrom),
                                                       seq_type)));
        }
        // Wait for all tasks (and propagate any exception).
        for (auto &f : results) f.get();
    }

    // 4. Compute length-normalized ratios and copy-number statistics.
    compute_normalized_ratios(chromosomes);

    // 5. Compute GC content per chromosome.
    _calculate_gc_content(chromosomes);

    // 6. Emit TSV results.
    _write_results(chromosomes, seq_type);
}
