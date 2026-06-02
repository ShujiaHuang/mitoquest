/**
 * @file trans_prep.h
 * @brief Extract mother-child mtDNA allele transmission pairs from a multi-
 *        sample VCF (produced by `mitoquest caller`) plus a PLINK-format FAM
 *        file describing the trios.
 *
 * The output TSV is consumed by `mitoquest ne-estimate` to fit the
 * mitochondrial bottleneck size (Ne) via Beta-Binomial maximum likelihood.
 *
 * Only mother-child relationships are retained; fathers do not transmit
 * mtDNA in mammals and are deliberately ignored.
 *
 * @author Shujia Huang (hshujia@qq.com)
 * @date 2026-05-28
 */
#ifndef _MT_TRANS_PREP_H_
#define _MT_TRANS_PREP_H_

#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>

#include "version.h"
#include "io/vcf.h"
#include "io/vcf_header.h"
#include "io/vcf_record.h"

class TransmissionPrep {
public:
    struct Config {
        std::string vcf_path;        // input multi-sample VCF/BCF
        std::string fam_path;        // input PLINK FAM file (mother-child pairs)
        std::string gm_fam_path;     // optional PLINK FAM file for grandmother-mother
                                     // pairs; when set, the trio-detection logic
                                     // looks up the grandmother of each MC pair's
                                     // mother here and emits the trio columns.
                                     // Empty => no trio extension (legacy output).
        std::string output_file;     // output TSV file (empty => stdout)
        int  min_depth;              // minimum DP for both mother and child
        bool require_pass;           // skip records whose FILTER is not PASS
        bool snv_only;               // restrict to bi-allelic / multi-allelic SNVs
    };

    // PLINK FAM record tied to VCF sample indices.
    struct Trio {
        std::string fam_id;
        std::string child_id;
        std::string father_id;       // kept for diagnostics only
        std::string mother_id;
        int child_idx;               // -1 if not in VCF header
        int mother_idx;              // -1 if not in VCF header

        // Populated by resolve_gm_for_trios() when the mother of this trio
        // appears as the CHILD in a separate GM FAM line; otherwise empty /
        // -1.  When grandmother_idx >= 0, the VCF row carries the G-M-C
        // trio data and HAS_G = 1.
        std::string grandmother_id;
        int         grandmother_idx  = -1;

        Trio() : child_idx(-1), mother_idx(-1) {}
    };

    // Statistics for the FAM <-> VCF matching report (printed to stderr).
    struct MatchingStats {
        int total_fam_lines           = 0;
        int total_vcf_samples         = 0;
        int valid_mother_child_pairs  = 0;
        int ignored_no_mother_in_fam  = 0;
        int ignored_missing_child     = 0;
        int ignored_missing_mother    = 0;
        // (pair, site) observations dropped because the mother and/or the
        // child has GT='.' at the record (or because GT and AD lengths
        // disagree, which prevents GT-aligned AD lookup).
        long long pair_site_dropped_gt_missing = 0;

        // GM FAM matching counters (zero when --gm-fam is not provided).
        int gm_total_fam_lines        = 0;
        int gm_matched_trios          = 0;  // MC pair whose mother has a GM in VCF
        int gm_ignored_no_mother      = 0;
        int gm_ignored_missing_child  = 0;  // "child" in GM FAM = the MC mother
        int gm_ignored_missing_mother = 0;  // "mother" in GM FAM = the grandmother
    };

    // One emitted output row.  The grandmother columns are filled only when
    // `grandmother_id` is non-empty; otherwise they are written as "NA" /
    // 0 / 0.0 / "NA" so that the TSV always has the same width (HAS_G = 0
    // tells downstream consumers to ignore the G columns).
    struct PairRecord {
        std::string chrom;
        int32_t     pos;             // 1-based
        std::string ref;
        std::string alt;
        std::string fam_id;
        std::string grandmother_id;  // empty => no grandmother in pedigree
        std::string mother_id;
        std::string child_id;
        int  g_dp = 0,    g_ad_ref = 0, g_ad_alt = 0;
        double g_vaf = 0.0;
        int  m_dp = 0,    c_dp = 0;
        int  m_ad_ref = 0, m_ad_alt = 0;
        int  c_ad_ref = 0, c_ad_alt = 0;
        double m_vaf = 0.0, c_vaf = 0.0;
        int  has_g = 0;              // 1 iff grandmother columns are populated
        std::string qc;              // "PASS" or reason for failure
    };

    explicit TransmissionPrep(int argc, char* argv[]);
    explicit TransmissionPrep(Config config);          // for unit tests / library callers
    ~TransmissionPrep() = default;

    // Main entry point: stream the VCF, emit the pair TSV.
    // Returns the number of PASS rows written.
    long long run();

    const Config& config() const { return _config; }
    const MatchingStats& matching_stats() const { return _stats; }

    // -----------------------------------------------------------------
    // Pure helpers, exposed for unit tests.
    // -----------------------------------------------------------------

    // Parse a PLINK FAM file and resolve each trio's child / mother to
    // the VCF sample indices present in `hdr`.  Trios whose mother or
    // child is missing in the VCF are dropped (and counted in `stats`).
    static std::vector<Trio> parse_fam(const std::string& fam_path,
                                       const ngslib::VCFHeader& hdr,
                                       MatchingStats& stats);

    // Parse the optional GM FAM file (same format as parse_fam; here the
    // "child" is the mother of an MC pair, and the "mother" is her
    // grandmother).  Returned trios have `child_id` = the MC mother and
    // `mother_id` = the grandmother.  The MatchingStats counters prefixed
    // with `gm_` are populated.
    static std::vector<Trio> parse_gm_fam(const std::string& fam_path,
                                          const ngslib::VCFHeader& hdr,
                                          MatchingStats& stats);

    // For each MC trio in `trios`, if its `mother_id` appears as the
    // `child_id` of any GM trio in `gm_trios`, copy the grandmother's
    // id + VCF index into the MC trio (grandmother_id / grandmother_idx).
    // Updates stats.gm_matched_trios for every match.
    static void resolve_gm_for_trios(std::vector<Trio>& trios,
                                     const std::vector<Trio>& gm_trios,
                                     MatchingStats& stats);

    // TSV header line written to the output file.  The header is always
    // the full wide (HAS_G-aware) form; rows that have HAS_G = 0 fill the
    // grandmother columns with NA / 0.
    static const char* tsv_header();

    // Format a single PairRecord as one TSV line (no trailing newline).
    static std::string format_row(const PairRecord& r);

private:
    TransmissionPrep(const TransmissionPrep&)            = delete;
    TransmissionPrep& operator=(const TransmissionPrep&) = delete;

    static void usage();
    void _parse_args(int argc, char* argv[]);
    void _print_matching_report(std::ostream& os) const;

    // Validate that the VCF declares FORMAT/AD with Number=R and FORMAT/DP.
    static void _validate_vcf_format(const ngslib::VCFHeader& hdr);

    Config        _config;
    MatchingStats _stats;
    std::string   _cmdline_string;
};

#endif // _MT_TRANS_PREP_H_
