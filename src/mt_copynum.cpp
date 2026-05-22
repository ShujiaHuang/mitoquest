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
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_map>

#include "external/thread_pool.h"  // ThreadPool
#include "io/utils.h"              // ngslib::is_readable, ngslib::split

namespace {

// Split a string on `sep`, dropping empty tokens.
std::vector<std::string> split_nonempty(const std::string &s, char sep) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == sep) {
            if (!cur.empty()) out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

// Trim surrounding ASCII whitespace in place.
void trim_inplace(std::string &s) {
    auto not_space = [](unsigned char c){ return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
}

// Parse a single "chr:start-end" / "chr:start" / "chr" token into a
// 1-based inclusive GenomeRegion. `chrom_length` (if provided) supplies
// the default end when only `chr` or `chr:start` is given.
ngslib::GenomeRegion parse_one_region_token(
    const std::string &token_in,
    const std::function<uint32_t(const std::string &)> &chrom_length)
{
    std::string token = token_in;
    trim_inplace(token);
    if (token.empty()) {
        throw std::invalid_argument("[copynum] Empty region token.");
    }

    std::string chrom;
    uint32_t    start = 1;
    uint32_t    end   = std::numeric_limits<uint32_t>::max();

    auto colon = token.find(':');
    if (colon == std::string::npos) {
        chrom = token;
    } else {
        chrom = token.substr(0, colon);
        std::string coord = token.substr(colon + 1);
        trim_inplace(coord);
        if (coord.empty()) {
            // "chr:" => whole chromosome
        } else {
            auto dash = coord.find('-');
            std::string s_str = (dash == std::string::npos) ? coord : coord.substr(0, dash);
            std::string e_str = (dash == std::string::npos) ? std::string() : coord.substr(dash + 1);
            trim_inplace(s_str);
            trim_inplace(e_str);
            try {
                if (!s_str.empty()) start = static_cast<uint32_t>(std::stoul(s_str));
            } catch (const std::exception &) {
                throw std::invalid_argument("[copynum] Bad region start in '" + token_in + "'");
            }
            if (!e_str.empty()) {
                try {
                    end = static_cast<uint32_t>(std::stoul(e_str));
                } catch (const std::exception &) {
                    throw std::invalid_argument("[copynum] Bad region end in '" + token_in + "'");
                }
            }
        }
    }

    if (chrom.empty()) {
        throw std::invalid_argument("[copynum] Missing chromosome in '" + token_in + "'");
    }

    if (end == std::numeric_limits<uint32_t>::max() && chrom_length) {
        uint32_t L = chrom_length(chrom);
        if (L > 0) end = L;
    }
    if (start == 0) start = 1;   // be tolerant: treat 0-based 0 as 1-based 1
    if (end == 0)   end   = start;
    if (start > end) {
        throw std::invalid_argument("[copynum] start > end in region '" + token_in + "'");
    }
    return ngslib::GenomeRegion(chrom, start, end);
}

}  // anonymous namespace

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
                 "  -L, --regions   STR    Restrict counting / GC / length normalization to\n"
                 "                         the listed regions, given either as a comma-\n"
                 "                         separated list (e.g. 'chrM:1-300,chrM:16000-16569')\n"
                 "                         or as a path to a file (one region per line, either\n"
                 "                         'chr:start-end' samtools form or a BED-style\n"
                 "                         'chr<TAB>start<TAB>end' triple; '#' starts a\n"
                 "                         comment). Chromosomes not covered by any region\n"
                 "                         are still measured over their full length.\n"
                 "                         Typical use: exclude NUMT-affected zones from the\n"
                 "                         mtDNA copy-number estimate.\n"
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
    _config.regions_arg.clear();
    _config.regions.clear();

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
        {"regions",   required_argument, 0, 'L'},
        {"help",      no_argument,       0, 'h'},
        {0, 0, 0, 0}  // Terminator
    };

    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "r:o:q:t:s:L:h", long_options, &option_index)) != -1) {
        switch (c) {
            case 'r': _config.reference_file = optarg;         break;
            case 'o': _config.output_file    = optarg;         break;
            case 'q': _config.min_mapq       = std::stoi(optarg); break;
            case 't': _config.thread_count   = std::stoi(optarg); break;
            case 'L': _config.regions_arg    = optarg;         break;
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

    // Parse regions argument (file path or comma-separated list).
    // The chromosome-length resolver runs later in run() once the BAM
    // header is known; here we accept open-ended tokens ("chrM" or
    // "chrM:1") and rely on _attach_regions_to_chroms to clamp them.
    if (!_config.regions_arg.empty()) {
        try {
            _config.regions = parse_regions_arg(_config.regions_arg, nullptr);
        } catch (const std::exception &e) {
            std::cerr << "Error: Failed to parse -L/--regions argument: "
                      << e.what() << "\n";
            usage();
            exit(EXIT_FAILURE);
        }
        if (_config.regions.empty()) {
            std::cerr << "Error: -L/--regions was given but no valid region was parsed.\n";
            usage();
            exit(EXIT_FAILURE);
        }
    }
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

    // If -L/--regions was supplied, attach the user-supplied intervals to
    // each chromosome and recompute its `effective_length`. Chromosomes
    // not covered by any region keep their default (full-length) behaviour.
    _attach_regions_to_chroms(chromosomes);

    return chromosomes;
}

// Distribute Config::regions into the matching ChromosomeData entries.
// Each region is clamped to the chromosome length and overlapping / adjacent
// intervals on the same chromosome are merged (1-based inclusive). After
// this call, `effective_length` is the sum of the merged interval lengths
// (or the full chromosome length when no region targets that chromosome).
void MtCopyNumber::_attach_regions_to_chroms(
    std::vector<ChromosomeData> &chromosomes) const
{
    if (_config.regions.empty()) return;

    // Build name -> index lookup once.
    std::unordered_map<std::string, size_t> idx;
    idx.reserve(chromosomes.size() * 2);
    for (size_t i = 0; i < chromosomes.size(); ++i) {
        idx.emplace(chromosomes[i].name, i);
    }

    bool any_match = false;
    for (const auto &r : _config.regions) {
        auto it = idx.find(r.chrom);
        if (it == idx.end()) {
            std::cerr << "[WARN] -L region references chromosome '" << r.chrom
                      << "' which is not in the BAM header; skipping." << std::endl;
            continue;
        }
        auto &cd = chromosomes[it->second];
        uint32_t s = r.start;
        uint32_t e = std::min<uint32_t>(r.end, cd.length);
        if (s < 1)          s = 1;
        if (s > cd.length)  continue;          // wholly past end
        if (s > e)          continue;          // empty after clamp
        cd.regions.emplace_back(cd.name, s, e);
        any_match = true;
    }

    if (!any_match) {
        throw std::runtime_error("[copynum] None of the regions supplied via -L/--regions "
                                 "matched any chromosome in the BAM/CRAM header.");
    }

    // Merge overlapping / adjacent intervals on each chromosome and
    // recompute effective_length.
    for (auto &cd : chromosomes) {
        if (cd.regions.empty()) continue;  // unrestricted => keep length
        std::sort(cd.regions.begin(), cd.regions.end(),
                  [](const ngslib::GenomeRegion &a, const ngslib::GenomeRegion &b){
                      return a.start < b.start;
                  });
        std::vector<ngslib::GenomeRegion> merged;
        merged.reserve(cd.regions.size());
        for (const auto &r : cd.regions) {
            if (!merged.empty() && r.start <= merged.back().end + 1) {
                merged.back().end = std::max(merged.back().end, r.end);
            } else {
                merged.push_back(r);
            }
        }
        cd.regions.swap(merged);

        uint64_t eff = 0;
        for (const auto &r : cd.regions) {
            eff += static_cast<uint64_t>(r.end) - r.start + 1;
        }
        cd.effective_length = (eff > std::numeric_limits<uint32_t>::max())
                                ? std::numeric_limits<uint32_t>::max()
                                : static_cast<uint32_t>(eff);
    }
}

// Count fragments for a single chromosome. Each call opens its own ngslib::Bam
// handle so the function is safe to invoke from worker threads in parallel.
//
// When `chrom_data.regions` is non-empty, only fragments whose 5' end (read1
// for PE, or any primary alignment for SE / orphan mates) falls inside one
// of those intervals are counted. This is what makes -L/--regions usable as
// a NUMT-exclusion filter for the mtDNA contig.
void MtCopyNumber::_count_chrom_fragments(ChromosomeData &chrom_data,
                                          SeqType seq_type) const {
    ngslib::Bam bf(_config.input_bam, "r", _config.reference_file);
    ngslib::BamHeader &hdr = bf.header();

    // Skip silently if the chromosome is absent in the BAM index (defensive).
    if (hdr.name2id(chrom_data.name) < 0) {
        return;
    }

    // Walk either the single full-chromosome interval or the user-supplied
    // sub-intervals. Both branches go through the same counting loop below.
    std::vector<ngslib::GenomeRegion> intervals;
    if (chrom_data.regions.empty()) {
        intervals.emplace_back(chrom_data.name, 1u, chrom_data.length);
    } else {
        intervals = chrom_data.regions;
    }

    auto count_one = [&](const ngslib::GenomeRegion &iv) {
        // ngslib::Bam::fetch(seq, beg, end) takes a 0-based half-open
        // interval ([beg, end)). The interval here is 1-based inclusive
        // so convert by subtracting 1 from the start.
        if (!bf.fetch(chrom_data.name,
                      static_cast<hts_pos_t>(iv.start) - 1,
                      static_cast<hts_pos_t>(iv.end))) {
            return;
        }

        ngslib::BamRecord rec;
        while (bf.read(rec) >= 0) {
            if (!rec.is_mapped() || rec.is_secondary() || rec.is_supplementary()) {
                continue;
            }
            if (rec.mapq() < _config.min_mapq) continue;

            // Anchor each fragment by its 5' end (map_ref_start_pos is
            // 0-based) so that a fragment is counted at most once across
            // the union of intervals, even though the BAM iterator may
            // return reads whose alignment extends across an interval edge.
            hts_pos_t anchor0 = rec.map_ref_start_pos();
            if (anchor0 < 0) continue;
            uint32_t anchor1 = static_cast<uint32_t>(anchor0) + 1;
            if (anchor1 < iv.start || anchor1 > iv.end) continue;

            if (seq_type == SeqType::SE) {
                ++chrom_data.count;
            } else {
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
    };

    for (const auto &iv : intervals) {
        count_one(iv);
    }

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
        // Use effective_length so that region-restricted contigs (e.g. mtDNA
        // with NUMT zones masked out via -L) contribute only the bases that
        // were actually measured.
        total_genome_length += chrom.effective_length;
    }

    if (total_fragments == 0 || total_genome_length == 0) {
        throw std::runtime_error("[copynum] No fragments counted; check input BAM/CRAM and "
                                 "filtering options (e.g. --mapq, --regions).");
    }

    for (auto &chrom : chromosomes) {
        double fragment_ratio  = static_cast<double>(chrom.count)            / total_fragments;
        double length_ratio    = static_cast<double>(chrom.effective_length) / total_genome_length;
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

// Compute GC content for each chromosome via ngslib::Fasta. When -L/--regions
// restricts a chromosome to a set of intervals, the GC fraction is computed
// (and length-weighted) over those intervals only; otherwise the whole
// chromosome is used.
void MtCopyNumber::_calculate_gc_content(std::vector<ChromosomeData> &chromosomes) const {
    ngslib::Fasta reference(_config.reference_file);

    for (auto &chrom : chromosomes) {
        if (!reference.has_seq(chrom.name)) {
            chrom.gc_content = 0.0;
            continue;
        }
        if (chrom.regions.empty()) {
            std::string seq = reference.fetch(chrom.name);
            chrom.gc_content = compute_gc_content(seq);
            continue;
        }
        // Length-weighted GC across all measured sub-intervals.
        int64_t gc = 0, valid = 0;
        for (const auto &iv : chrom.regions) {
            // Fasta::fetch(chr, start, end) takes 1-based start, 1-based end inclusive
            // per faidx semantics in this codebase.
            std::string seq = reference.fetch(chrom.name, iv.start, iv.end);
            for (char c : seq) {
                char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                if (u == 'N') continue;
                ++valid;
                if (u == 'G' || u == 'C') ++gc;
            }
        }
        chrom.gc_content = (valid > 0) ? static_cast<double>(gc) / valid : 0.0;
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
    if (!_config.regions.empty()) {
        *out << "#Regions argument: " << _config.regions_arg << "\n";
    }
    *out << "#Chromosome\tFragments\tChrom_Length\tGC_Content\t"
            "Fragment_Normalized_Ratio\tCopyNum\tCopyNum-CI95-Lower\tCopyNum-CI95-Upper\t"
            "Effective_Length\tRegions_Used\n";

    for (const auto &chrom : chromosomes) {
        // Compact representation of the per-chromosome intervals (or "." if
        // the whole chromosome was used).
        std::string regions_str = ".";
        if (!chrom.regions.empty()) {
            std::ostringstream oss;
            for (size_t i = 0; i < chrom.regions.size(); ++i) {
                if (i) oss << ",";
                oss << chrom.regions[i].start << "-" << chrom.regions[i].end;
            }
            regions_str = oss.str();
        }
        *out << chrom.name   << "\t"
             << chrom.count  << "\t"
             << chrom.length << "\t"
             << std::fixed   << std::setprecision(4) << chrom.gc_content        << "\t"
             << std::fixed   << std::setprecision(4) << chrom.normalized_ratio  << "\t"
             << std::fixed   << std::setprecision(4) << chrom.cn_stats.mean     << "\t"
             << std::fixed   << std::setprecision(4) << chrom.cn_stats.ci_lower << "\t"
             << std::fixed   << std::setprecision(4) << chrom.cn_stats.ci_upper << "\t"
             << chrom.effective_length << "\t"
             << regions_str << "\n";
    }
}

// Parse a `-L/--regions` value into a vector of 1-based inclusive
// GenomeRegion intervals. Accepts either a file path or a comma-separated
// list of `chr:start-end` tokens; the file form also accepts BED-style
// `chr<TAB>start<TAB>end` lines (0-based half-open, converted to 1-based
// inclusive).
std::vector<ngslib::GenomeRegion> MtCopyNumber::parse_regions_arg(
    const std::string &arg,
    const std::function<uint32_t(const std::string &)> &chrom_length)
{
    std::vector<ngslib::GenomeRegion> out;
    if (arg.empty()) return out;

    if (ngslib::is_readable(arg)) {
        std::ifstream ifs(arg);
        if (!ifs.is_open()) {
            throw std::runtime_error("[copynum] Could not open regions file: " + arg);
        }
        std::string line;
        while (std::getline(ifs, line)) {
            // strip comments and surrounding whitespace
            auto hash = line.find('#');
            if (hash != std::string::npos) line.erase(hash);
            trim_inplace(line);
            if (line.empty()) continue;

            // Auto-detect BED form: at least one tab AND at least three
            // tab-separated fields whose 2nd and 3rd parse as integers.
            if (line.find('\t') != std::string::npos) {
                std::vector<std::string> cols;
                std::string cur;
                for (char c : line) {
                    if (c == '\t') { cols.push_back(cur); cur.clear(); }
                    else           { cur.push_back(c); }
                }
                if (!cur.empty()) cols.push_back(cur);
                if (cols.size() >= 3) {
                    try {
                        uint32_t s0 = static_cast<uint32_t>(std::stoul(cols[1])); // 0-based start
                        uint32_t e0 = static_cast<uint32_t>(std::stoul(cols[2])); // 0-based half-open end
                        if (e0 <= s0) {
                            throw std::invalid_argument("[copynum] empty BED interval in '" + line + "'");
                        }
                        std::string chrom = cols[0];
                        trim_inplace(chrom);
                        out.emplace_back(chrom, s0 + 1u, e0);  // BED -> 1-based inclusive
                        continue;
                    } catch (const std::invalid_argument &) {
                        // Not actually BED; fall through to samtools-style parser below.
                    } catch (const std::out_of_range &) {
                        throw std::runtime_error("[copynum] BED interval out of range: " + line);
                    }
                }
            }
            // samtools-style "chr:start-end" per-line form.
            out.push_back(parse_one_region_token(line, chrom_length));
        }
    } else {
        for (auto tok : split_nonempty(arg, ',')) {
            trim_inplace(tok);
            if (tok.empty()) continue;  // ignore whitespace-only tokens
            out.push_back(parse_one_region_token(tok, chrom_length));
        }
    }
    return out;
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
