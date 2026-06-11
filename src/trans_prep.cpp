/**
 * @file trans_prep.cpp
 * @brief Implementation of the `mitoquest trans-prep` subcommand.
 *
 * @author Shujia Huang (hshujia@qq.com)
 * @date 2026-05-28
 */
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#include "trans_prep.h"
#include "io/utils.h"           // ngslib::is_readable

namespace {

// Trim surrounding ASCII whitespace in place.
void trim_inplace(std::string &s) {
    auto not_space = [](unsigned char c){ return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
}

// Returns true when the FAM token denotes "missing parent / unknown".
bool is_missing_id(const std::string& s) {
    return s.empty() || s == "0" || s == "." || s == "-9";
}

// Format a double with up to `precision` digits, trimming trailing zeros.
std::string fmt_double(double v, int precision = 6) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << v;
    std::string out = oss.str();
    // strip trailing zeros and a possibly orphaned trailing dot
    auto dot = out.find('.');
    if (dot != std::string::npos) {
        size_t last = out.find_last_not_of('0');
        if (last == dot) last--;            // keep at least "0"
        out.erase(last + 1);
    }
    return out;
}

}  // anonymous namespace

// ---------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------

TransmissionPrep::TransmissionPrep(int argc, char* argv[]) {
    _parse_args(argc, argv);
}

TransmissionPrep::TransmissionPrep(Config config) : _config(std::move(config)) {
    if (_config.vcf_path.empty() || _config.fam_path.empty()) {
        throw std::invalid_argument("[trans-prep] vcf_path and fam_path are required.");
    }
}

// ---------------------------------------------------------------------
// Static helpers
// ---------------------------------------------------------------------

const char* TransmissionPrep::tsv_header() {
    return "CHROM\tPOS\tREF\tALT\tFAM_ID\tGRANDMOTHER_ID\tMOTHER_ID\tCHILD_ID\t"
           "GRANDMOTHER_DP\tGRANDMOTHER_AD_REF\tGRANDMOTHER_AD_ALT\tGRANDMOTHER_VAF\t"
           "HAS_G\t"
           "MOTHER_DP\tMOTHER_AD_REF\tMOTHER_AD_ALT\tMOTHER_VAF\t"
           "CHILD_DP\tCHILD_AD_REF\tCHILD_AD_ALT\tCHILD_VAF\tQC";
}

std::string TransmissionPrep::format_row(const PairRecord& r) {
    std::ostringstream oss;
    oss << r.chrom    << "\t" << r.pos       << "\t" << r.ref       << "\t" << r.alt       << "\t"
        << r.fam_id   << "\t"
        << (r.grandmother_id.empty() ? "NA" : r.grandmother_id) << "\t"
        << r.mother_id << "\t" << r.child_id << "\t"
        << r.g_dp     << "\t" << r.g_ad_ref  << "\t" << r.g_ad_alt  << "\t" << fmt_double(r.g_vaf) << "\t"
        << r.has_g    << "\t"
        << r.m_dp     << "\t" << r.m_ad_ref  << "\t" << r.m_ad_alt  << "\t" << fmt_double(r.m_vaf) << "\t"
        << r.c_dp     << "\t" << r.c_ad_ref  << "\t" << r.c_ad_alt  << "\t" << fmt_double(r.c_vaf) << "\t"
        << r.qc;
    return oss.str();
}

std::vector<TransmissionPrep::Trio>
TransmissionPrep::parse_fam(const std::string& fam_path,
                            const ngslib::VCFHeader& hdr,
                            MatchingStats& stats) {
    std::vector<Trio> valid_trios;
    std::ifstream fam_file(fam_path);
    if (!fam_file.is_open()) {
        throw std::runtime_error("[trans-prep] Failed to open FAM file: " + fam_path);
    }
    stats.total_vcf_samples = hdr.n_samples();

    std::string line;
    while (std::getline(fam_file, line)) {
        // Strip CR for files with Windows line endings.
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        Trio t;
        if (!(iss >> t.fam_id >> t.child_id >> t.father_id >> t.mother_id)) {
            std::cerr << "[trans-prep][Warning] Malformed FAM line, skipping: "
                      << line << "\n";
            continue;
        }
        stats.total_fam_lines++;

        if (is_missing_id(t.mother_id)) {
            stats.ignored_no_mother_in_fam++;
            continue;
        }

        t.child_idx  = hdr.sample_index(t.child_id);
        t.mother_idx = hdr.sample_index(t.mother_id);

        const bool child_in_vcf  = (t.child_idx  >= 0);
        const bool mother_in_vcf = (t.mother_idx >= 0);

        if (child_in_vcf && mother_in_vcf) {
            valid_trios.push_back(t);
            stats.valid_mother_child_pairs++;
        } else {
            if (!child_in_vcf)  stats.ignored_missing_child++;
            if (!mother_in_vcf) stats.ignored_missing_mother++;
        }
    }
    return valid_trios;
}

std::vector<TransmissionPrep::Trio>
TransmissionPrep::parse_gm_fam(const std::string& fam_path,
                               const ngslib::VCFHeader& hdr,
                               MatchingStats& stats) {
    // Same structure as parse_fam but updates the gm_* counters.
    std::vector<Trio> valid;
    std::ifstream fam_file(fam_path);
    if (!fam_file.is_open()) {
        throw std::runtime_error("[trans-prep] Failed to open GM FAM file: " + fam_path);
    }

    std::string line;
    while (std::getline(fam_file, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        Trio t;
        if (!(iss >> t.fam_id >> t.child_id >> t.father_id >> t.mother_id)) {
            std::cerr << "[trans-prep][Warning] Malformed GM FAM line, skipping: "
                      << line << "\n";
            continue;
        }
        stats.gm_total_fam_lines++;

        if (is_missing_id(t.mother_id)) {
            stats.gm_ignored_no_mother++;
            continue;
        }

        t.child_idx  = hdr.sample_index(t.child_id);   // this is the MC-mother
        t.mother_idx = hdr.sample_index(t.mother_id);  // this is the grandmother

        const bool child_in_vcf  = (t.child_idx  >= 0);
        const bool mother_in_vcf = (t.mother_idx >= 0);

        if (child_in_vcf && mother_in_vcf) {
            valid.push_back(t);
        } else {
            if (!child_in_vcf)  stats.gm_ignored_missing_child++;
            if (!mother_in_vcf) stats.gm_ignored_missing_mother++;
        }
    }
    return valid;
}

void TransmissionPrep::resolve_gm_for_trios(std::vector<Trio>& trios,
                                            const std::vector<Trio>& gm_trios,
                                            MatchingStats& stats) {
    if (gm_trios.empty()) return;

    // Build a lookup: mc_mother_id -> (grandmother_id, grandmother_idx).
    // When several GM FAM lines refer to the same MC-mother, keep the first
    // valid one (the GM FAM is expected to be 1:1 per mother).
    std::unordered_map<std::string, std::pair<std::string, int>> gm_by_mother;
    gm_by_mother.reserve(gm_trios.size());
    for (const auto& gm : gm_trios) {
        auto it = gm_by_mother.find(gm.child_id);
        if (it == gm_by_mother.end()) {
            gm_by_mother.emplace(gm.child_id,
                                 std::make_pair(gm.mother_id, gm.mother_idx));
        }
    }

    for (auto& t : trios) {
        auto it = gm_by_mother.find(t.mother_id);
        if (it != gm_by_mother.end()) {
            t.grandmother_id  = it->second.first;
            t.grandmother_idx = it->second.second;
            ++stats.gm_matched_trios;
        }
    }
}

void TransmissionPrep::_validate_vcf_format(const ngslib::VCFHeader& hdr) {
    bcf_hdr_t* h = hdr.hts_header();
    if (!h) throw std::runtime_error("[trans-prep] Invalid VCF header.");

    int gt_id = bcf_hdr_id2int(h, BCF_DT_ID, "GT");
    int dp_id = bcf_hdr_id2int(h, BCF_DT_ID, "DP");
    int ad_id = bcf_hdr_id2int(h, BCF_DT_ID, "AD");
    if (gt_id < 0 || !bcf_hdr_idinfo_exists(h, BCF_HL_FMT, gt_id)) {
        throw std::runtime_error("[trans-prep] Input VCF lacks FORMAT/GT.");
    }
    if (dp_id < 0 || !bcf_hdr_idinfo_exists(h, BCF_HL_FMT, dp_id)) {
        throw std::runtime_error("[trans-prep] Input VCF lacks FORMAT/DP. "
                                 "Use a VCF produced by `mitoquest caller`.");
    }
    if (ad_id < 0 || !bcf_hdr_idinfo_exists(h, BCF_HL_FMT, ad_id)) {
        throw std::runtime_error("[trans-prep] Input VCF lacks FORMAT/AD. "
                                 "Use a VCF produced by `mitoquest caller`.");
    }
    // We do *not* enforce a particular Number= cardinality on AD.  The
    // `mitoquest caller` declares AD as `Number=A` but emits a per-sample
    // GT-aligned layout: AD[i] is the read depth of the allele at GT
    // position i.  Other tooling may declare `Number=R` (REF + all ALTs)
    // or `Number=.` (variable).  We always interpret AD via FORMAT/GT, so
    // any declaration is acceptable here.
}

// ---------------------------------------------------------------------
// CLI handling
// ---------------------------------------------------------------------

void TransmissionPrep::usage() {
    std::cerr << "Usage: mitoquest trans-prep [options] -v <multi.vcf.gz> -f <pedigree.fam>\n\n"
                 "Description:\n"
                 "  Extract mother-child mtDNA allele transmission pairs from a multi-\n"
                 "  sample VCF produced by `mitoquest caller`, joined to a PLINK FAM\n"
                 "  file describing the trios. Only mother-child relationships are\n"
                 "  retained; fathers are ignored (mtDNA is maternally inherited).\n"
                 "\n"
                 "  The TSV emitted here is the input of `mitoquest ne-estimate`.\n"
                 "\n"
                 "  When the optional --gm-fam file is provided, each MC pair's mother\n"
                 "  is looked up as the CHILD in that FAM; if found, the corresponding\n"
                 "  MOTHER is the pedigree grandmother, and the output row gains the\n"
                 "  GRANDMOTHER_* columns and HAS_G = 1 (a 3-generation G-M-C trio).\n"
                 "  Pairs without a matched grandmother are emitted with HAS_G = 0\n"
                 "  and NA/0 fills, so the downstream `ne-estimate` module can fit\n"
                 "  Ne from a mixed cohort of 2-gen and 3-gen pedigrees.\n"
                 "\nRequired options:\n"
                 "  -v, --vcf       FILE   Input multi-sample VCF/BCF.\n"
                 "  -f, --fam       FILE   PLINK FAM file (whitespace-delimited):\n"
                 "                         Family_ID Child_ID Father_ID Mother_ID Sex Phenotype\n"
                 "                         Missing parents may be encoded as 0 / . / -9.\n"
                 "\nOptional options:\n"
                 "  -g, --gm-fam    FILE   Second PLINK FAM file describing\n"
                 "                         grandmother -> mother transmissions.  Its\n"
                 "                         CHILD_ID must equal an MC FAM's MOTHER_ID\n"
                 "                         for the corresponding pair to gain trio\n"
                 "                         columns.  When absent, every row has\n"
                 "                         HAS_G = 0 and the GRANDMOTHER_* columns\n"
                 "                         are filled with NA / 0.\n"
                 "  -o, --output    FILE   Output TSV file (default: stdout).\n"
                 "  -d, --min-depth INT    Minimum read depth (DP) required for BOTH the\n"
                 "                         mother and the child at a site.  Sites failing\n"
                 "                         this threshold are emitted with QC=LOW_DEPTH and\n"
                 "                         not exported as PASS rows [500].\n"
                 "      --require-pass     Keep only sites whose VCF FILTER is PASS [on].\n"
                 "      --no-require-pass  Disable the FILTER=PASS gate.\n"
                 "      --snv-only         Restrict to single-nucleotide variants [on].\n"
                 "      --no-snv-only      Disable the SNV-only restriction.\n"
                 "  -h, --help             Print this help message.\n\n"
                 "AD interpretation:\n"
                 "  FORMAT/AD is decoded *per sample* using FORMAT/GT.  AD[i] is\n"
                 "  the read depth of the allele at GT position i; alleles\n"
                 "  absent from a sample's GT contribute 0 reads for that\n"
                 "  sample.  This matches the layout emitted by `mitoquest\n"
                 "  caller` and also handles standard Number=R AD correctly.\n"
                 "  Pairs where mother *or* child has GT='.' at a site are\n"
                 "  dropped from that site (counted on STDERR).\n\n"
                 "[WARNING] Garbage In, Garbage Out:\n"
                 "  Sample IDs are matched between FAM and VCF case-sensitively.  Trios\n"
                 "  with mismatched IDs are silently dropped.  Always inspect the\n"
                 "  matching report on stderr.\n\n"
              << "Version: " << MITOQUEST_VERSION << "\n"
              << std::endl;
}

void TransmissionPrep::_parse_args(int argc, char* argv[]) {
    // Defaults
    _config.vcf_path.clear();
    _config.fam_path.clear();
    _config.gm_fam_path.clear();
    _config.output_file.clear();
    _config.min_depth    = 500;
    _config.require_pass = true;
    _config.snv_only     = true;

    _cmdline_string = "#mitoquest_trans_prep_command=";
    for (int i = 0; i < argc; ++i) {
        _cmdline_string += (i > 0) ? " " + std::string(argv[i]) : std::string(argv[i]);
    }

    static const struct option long_options[] = {
        {"vcf",              required_argument, 0, 'v'},
        {"fam",              required_argument, 0, 'f'},
        {"gm-fam",           required_argument, 0, 'g'},
        {"output",           required_argument, 0, 'o'},
        {"min-depth",        required_argument, 0, 'd'},
        {"require-pass",     no_argument,       0,  1 },
        {"no-require-pass",  no_argument,       0,  2 },
        {"snv-only",         no_argument,       0,  3 },
        {"no-snv-only",      no_argument,       0,  4 },
        {"help",             no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    int c;
    optind = 1;
    while ((c = getopt_long(argc, argv, "v:f:g:o:d:h", long_options, &option_index)) != -1) {
        switch (c) {
            case 'v': _config.vcf_path     = optarg;            break;
            case 'f': _config.fam_path     = optarg;            break;
            case 'g': _config.gm_fam_path  = optarg;            break;
            case 'o': _config.output_file  = optarg;            break;
            case 'd': _config.min_depth    = std::stoi(optarg); break;
            case 1:   _config.require_pass = true;              break;
            case 2:   _config.require_pass = false;             break;
            case 3:   _config.snv_only     = true;              break;
            case 4:   _config.snv_only     = false;             break;
            case 'h': usage(); std::exit(EXIT_SUCCESS);
            case '?':
            default:  usage(); std::exit(EXIT_FAILURE);
        }
    }

    if (_config.vcf_path.empty() || _config.fam_path.empty()) {
        std::cerr << "[trans-prep] Missing required option: -v/--vcf and/or -f/--fam.\n";
        usage();
        std::exit(EXIT_FAILURE);
    }
    if (!ngslib::is_readable(_config.vcf_path.c_str())) {
        throw std::runtime_error("[trans-prep] VCF file not readable: " + _config.vcf_path);
    }
    if (!ngslib::is_readable(_config.fam_path.c_str())) {
        throw std::runtime_error("[trans-prep] FAM file not readable: " + _config.fam_path);
    }
    if (!_config.gm_fam_path.empty() &&
        !ngslib::is_readable(_config.gm_fam_path.c_str())) {
        throw std::runtime_error("[trans-prep] GM FAM file not readable: " + _config.gm_fam_path);
    }
    if (_config.min_depth < 0) {
        throw std::runtime_error("[trans-prep] --min-depth must be >= 0.");
    }
}

void TransmissionPrep::_print_matching_report(std::ostream& os) const {
    os << "\n========== mitoquest trans-prep matching report ==========\n"
       << "Total samples in VCF header        : " << _stats.total_vcf_samples         << "\n"
       << "Total valid lines in MC FAM file   : " << _stats.total_fam_lines           << "\n"
       << "----------------------------------------------------------\n"
       << "Successfully matched mother-child  : " << _stats.valid_mother_child_pairs  << "\n"
       << "Ignored (child missing in VCF)     : " << _stats.ignored_missing_child     << "\n"
       << "Ignored (mother missing in VCF)    : " << _stats.ignored_missing_mother    << "\n"
       << "Ignored (no mother defined in FAM) : " << _stats.ignored_no_mother_in_fam  << "\n";
    if (!_config.gm_fam_path.empty()) {
        os << "----------------------------------------------------------\n"
           << "GM FAM file                        : " << _config.gm_fam_path          << "\n"
           << "Total valid lines in GM FAM file   : " << _stats.gm_total_fam_lines      << "\n"
           << "MC pairs matched to a grandmother  : " << _stats.gm_matched_trios        << "\n"
           << "Ignored (no grandmother in GM FAM) : " << _stats.gm_ignored_no_mother    << "\n"
           << "Ignored (MC-mother missing in VCF) : " << _stats.gm_ignored_missing_child<< "\n"
           << "Ignored (grandmother missing in VCF): " << _stats.gm_ignored_missing_mother << "\n";
    }
    os << "==========================================================\n\n";
}

// ---------------------------------------------------------------------
// run()
// ---------------------------------------------------------------------

long long TransmissionPrep::run() {
    ngslib::VCFFile reader(_config.vcf_path);
    if (!reader.is_open()) {
        throw std::runtime_error("[trans-prep] Failed to open VCF: " + _config.vcf_path);
    }
    ngslib::VCFHeader& hdr = reader.header();
    _validate_vcf_format(hdr);

    std::vector<Trio> trios = parse_fam(_config.fam_path, hdr, _stats);

    // Optional 3-gen extension: parse the GM FAM and resolve each MC trio's
    // mother against it.  Tries whose mother appears as a CHILD in the GM FAM
    // gain a matched grandmother (grandmother_idx >= 0) and are emitted with
    // HAS_G = 1 plus the GRANDMOTHER_* columns populated from the VCF.
    std::vector<Trio> gm_trios;
    if (!_config.gm_fam_path.empty()) {
        gm_trios = parse_gm_fam(_config.gm_fam_path, hdr, _stats);
        resolve_gm_for_trios(trios, gm_trios, _stats);
    }

    _print_matching_report(std::cerr);

    if (trios.empty()) {
        throw std::runtime_error("[trans-prep] No valid mother-child pairs found.");
    }

    // Open output (file or stdout).
    std::ofstream out_stream;
    std::ostream* out = &std::cout;
    if (!_config.output_file.empty()) {
        out_stream.open(_config.output_file);
        if (!out_stream.is_open()) {
            throw std::runtime_error("[trans-prep] Failed to open output: "
                                     + _config.output_file);
        }
        out = &out_stream;
    }

    // Header and provenance.
    if (!_cmdline_string.empty()) (*out) << _cmdline_string << "\n";
    (*out) << tsv_header() << "\n";

    long long total_records   = 0;
    long long total_emitted   = 0;
    long long total_low_depth = 0;

    // Per-record buffers (reused across iterations).
    std::vector<int32_t> dp_all;     // length = n_samples (one DP per sample)
    std::vector<int32_t> ad_all;     // length = n_samples * ad_per_sample,
                                     //   padded with bcf_int32_vector_end
                                     //   when individual samples emit fewer
                                     //   AD values (the caller's GT-aligned
                                     //   AD layout is per-sample variable).

    // GT buffer + ploidy width for the current record.  htslib stores GT
    // as int32 codes; missing alleles are signalled by
    // `bcf_gt_is_missing(code) == 1`.  See htslib/vcf.h.
    int32_t* gt_arr   = nullptr;
    int      gt_arr_n = 0;
    int      gt_ploidy = 0;

    ngslib::VCFRecord rec;
    while (reader.read(rec) >= 0) {
        rec.unpack(BCF_UN_ALL);
        if (_config.require_pass && !rec.passed_filters(hdr)) continue;
        if (rec.n_alt() == 0) continue;

        bcf1_t* b = rec.hts_record();
        if (!b) continue;

        // SNV-only filter.  We only require that the site contains *at least*
        // one SNV ALT — mixed multi-allelic sites such as A>G,GT are kept and
        // the per-ALT loop below skips the indel ALT(s) individually.
        if (_config.snv_only) {
            if (!(bcf_get_variant_types(b) & VCF_SNP)) continue;
            if (rec.ref().size() != 1) continue;
        }
        total_records++;

        // Hoisted DP / AD / GT extraction (exactly once per record).
        // Note: ngslib::VCFRecord::get_format_int returns the per-sample
        // value count (not the total), and writes the flat per-sample-major
        // buffer into the supplied vector.
        int dp_per_sample = rec.get_format_int(hdr, "DP", dp_all);
        if (dp_per_sample != 1) continue;       // DP missing or multi-valued

        int ad_per_sample = rec.get_format_int(hdr, "AD", ad_all);
        const int n_samples = hdr.n_samples();
        if (ad_per_sample <= 0 || n_samples <= 0) continue;

        const int n_alleles = b->n_allele;

        // FORMAT/GT is mandatory under the GT-aligned AD interpretation.
        int n_gt = bcf_get_genotypes(hdr.hts_header(), b, &gt_arr, &gt_arr_n);
        if (n_gt <= 0) continue;                 // GT field absent at this record
        gt_ploidy = n_gt / n_samples;
        if (gt_ploidy <= 0) continue;

        // GT-aligned AD lookup.  Returns:
        //   bcf_int32_missing  if GT is '.', or AD/GT lengths disagree, or
        //                      a sample-level malformedness is detected
        //                      (caller drops the pair at this site).
        //   0                  if `target_allele_idx` is not present in this
        //                      sample's GT (i.e. the caller did not record
        //                      any reads supporting that allele in this
        //                      sample's call).
        //   AD value           the read depth of `target_allele_idx` for
        //                      this sample, taken from the GT-matched AD
        //                      slot.
        auto sample_allele_depth = [&](int s_idx, int target_allele_idx) -> int32_t {
            const int32_t* g = gt_arr        + s_idx * gt_ploidy;
            const int32_t* d = ad_all.data() + s_idx * ad_per_sample;

            // Walk this sample's GT positions.
            int gt_len = 0;
            for (int p = 0; p < gt_ploidy; ++p) {
                if (g[p] == bcf_int32_vector_end) break;
                if (bcf_gt_is_missing(g[p]))      return bcf_int32_missing;
                ++gt_len;
            }
            if (gt_len == 0) return bcf_int32_missing;

            // Walk this sample's AD values.
            int ad_len = 0;
            for (int p = 0; p < ad_per_sample; ++p) {
                if (d[p] == bcf_int32_vector_end) break;
                ++ad_len;
            }
            // Strict: AD length must equal GT length, otherwise we cannot
            // unambiguously map AD positions to alleles.
            if (ad_len != gt_len) return bcf_int32_missing;

            // Find the GT position whose allele index matches the target.
            for (int p = 0; p < gt_len; ++p) {
                if (bcf_gt_allele(g[p]) == target_allele_idx) {
                    if (d[p] == bcf_int32_missing) return bcf_int32_missing;
                    return d[p];
                }
            }
            return 0;     // target allele not in GT ⇒ 0 supporting reads
        };

        const std::string chrom = rec.chrom(hdr);
        const int32_t     pos1  = static_cast<int32_t>(rec.pos() + 1);
        const std::string ref   = rec.ref();
        const std::vector<std::string> alts = rec.alt();

        for (int alt_idx = 1; alt_idx < n_alleles; ++alt_idx) {
            const std::string& alt = alts[alt_idx - 1];
            // SNV alt-length check (in case multi-allelic site mixes SNV/indel).
            if (_config.snv_only && alt.size() != 1) continue;

            for (const Trio& t : trios) {
                const int m_idx = t.mother_idx;
                const int c_idx = t.child_idx;
                const int g_idx = t.grandmother_idx;

                const int32_t m_dp = dp_all[m_idx];
                const int32_t c_dp = dp_all[c_idx];
                if (m_dp == bcf_int32_missing || c_dp == bcf_int32_missing) continue;
                if (m_dp <= 0 || c_dp <= 0) continue;

                // GT-aligned AD lookup.  Any GT='.' or AD/GT mismatch in
                // either sample drops the pair at this site (counted).
                const int32_t m_ref = sample_allele_depth(m_idx, 0);
                const int32_t m_alt = sample_allele_depth(m_idx, alt_idx);
                const int32_t c_ref = sample_allele_depth(c_idx, 0);
                const int32_t c_alt = sample_allele_depth(c_idx, alt_idx);
                if (m_ref == bcf_int32_missing || m_alt == bcf_int32_missing ||
                    c_ref == bcf_int32_missing || c_alt == bcf_int32_missing) {
                    ++_stats.pair_site_dropped_gt_missing;
                    continue;
                }

                PairRecord row;
                row.chrom     = chrom;
                row.pos       = pos1;
                row.ref       = ref;
                row.alt       = alt;
                row.fam_id    = t.fam_id;
                row.mother_id = t.mother_id;
                row.child_id  = t.child_id;
                row.m_dp      = m_dp;
                row.c_dp      = c_dp;
                row.m_ad_ref  = m_ref;
                row.m_ad_alt  = m_alt;
                row.c_ad_ref  = c_ref;
                row.c_ad_alt  = c_alt;
                row.m_vaf     = (m_dp > 0) ? static_cast<double>(m_alt) / m_dp : 0.0;
                row.c_vaf     = (c_dp > 0) ? static_cast<double>(c_alt) / c_dp : 0.0;

                // Trio columns.  HAS_G = 1 requires both the pedigree link AND
                // a usable DP / AD / GT at the grandmother sample; otherwise
                // we fall back to HAS_G = 0 with NA / 0 fills so the row is
                // still valid as a 2-gen MC observation.
                if (g_idx >= 0) {
                    const int32_t g_dp  = dp_all[g_idx];
                    const int32_t g_ref = sample_allele_depth(g_idx, 0);
                    const int32_t g_alt = sample_allele_depth(g_idx, alt_idx);
                    if (g_dp  != bcf_int32_missing && g_dp > 0 &&
                        g_ref != bcf_int32_missing &&
                        g_alt != bcf_int32_missing) {
                        row.grandmother_id = t.grandmother_id;
                        row.g_dp      = g_dp;
                        row.g_ad_ref  = g_ref;
                        row.g_ad_alt  = g_alt;
                        row.g_vaf     = static_cast<double>(g_alt) / g_dp;
                        row.has_g     = 1;
                    }
                }

                const bool low_depth = (m_dp < _config.min_depth) ||
                                       (c_dp < _config.min_depth);
                row.qc = low_depth ? "LOW_DEPTH" : "PASS";
                if (low_depth) ++total_low_depth;

                (*out) << format_row(row) << "\n";
                if (!low_depth) ++total_emitted;
            }
        }
    }

    if (out_stream.is_open()) out_stream.close();
    if (gt_arr) free(gt_arr);

    std::cerr << "[trans-prep] Processed "  << total_records  << " variant records.\n"
              << "[trans-prep] Wrote "      << (total_emitted + total_low_depth)
              << " rows ("                  << total_emitted  << " PASS, "
              << total_low_depth            << " LOW_DEPTH).\n"
              << "[trans-prep] (pair, site) dropped due to GT='.' or AD/GT mismatch: "
              << _stats.pair_site_dropped_gt_missing << "\n";

    return total_emitted;
}
