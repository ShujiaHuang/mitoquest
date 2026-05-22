// Author: Shujia Huang
// Date: 2026-05-22
//
// Unit tests for the `mitoquest copynum` subcommand (src/mt_copynum.{h,cpp}).
//
// Coverage:
//   * mt_utils helpers `is_autosomal` / `is_mitochondrial`
//   * Pure-math helpers `compute_gc_content`, `compute_statistics`,
//     `compute_normalized_ratios`
//   * End-to-end `run()` against tests/data/smp2_subset.bam (chrM + decoy contig)
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "mt_copynum.h"
#include "mt_utils.h"

// =====================================================================
//  mt_utils helpers
// =====================================================================

TEST(MtUtilsLookup, IsMitochondrial) {
    EXPECT_TRUE (is_mitochondrial("chrM"));
    EXPECT_TRUE (is_mitochondrial("chrMT"));
    EXPECT_TRUE (is_mitochondrial("MT"));
    EXPECT_TRUE (is_mitochondrial("M"));

    EXPECT_FALSE(is_mitochondrial("chr1"));
    EXPECT_FALSE(is_mitochondrial("Mt"));    // case-sensitive on purpose
    EXPECT_FALSE(is_mitochondrial("chrm"));
    EXPECT_FALSE(is_mitochondrial(""));
    EXPECT_FALSE(is_mitochondrial("NUMT_decoy"));
}

TEST(MtUtilsLookup, IsAutosomal) {
    EXPECT_TRUE (is_autosomal("1"));
    EXPECT_TRUE (is_autosomal("9"));
    EXPECT_TRUE (is_autosomal("22"));
    EXPECT_TRUE (is_autosomal("chr1"));
    EXPECT_TRUE (is_autosomal("chr22"));

    EXPECT_FALSE(is_autosomal("0"));
    EXPECT_FALSE(is_autosomal("23"));
    EXPECT_FALSE(is_autosomal("chrX"));
    EXPECT_FALSE(is_autosomal("chrY"));
    EXPECT_FALSE(is_autosomal("chrM"));
    EXPECT_FALSE(is_autosomal(""));
    EXPECT_FALSE(is_autosomal("chr1a"));
}

// =====================================================================
//  MtCopyNumber::compute_gc_content
// =====================================================================

TEST(MtCopyNumGcContent, AllInformativeBases) {
    EXPECT_DOUBLE_EQ(MtCopyNumber::compute_gc_content("GCGC"), 1.0);
    EXPECT_DOUBLE_EQ(MtCopyNumber::compute_gc_content("ATAT"), 0.0);
    EXPECT_DOUBLE_EQ(MtCopyNumber::compute_gc_content("AGCT"), 0.5);
}

TEST(MtCopyNumGcContent, IgnoresNs) {
    // Ns are excluded from BOTH the numerator and denominator.
    EXPECT_DOUBLE_EQ(MtCopyNumber::compute_gc_content("NNNN"),  0.0);
    EXPECT_DOUBLE_EQ(MtCopyNumber::compute_gc_content("AGCTNN"), 0.5);
    EXPECT_DOUBLE_EQ(MtCopyNumber::compute_gc_content("GNNN"),   1.0);
}

TEST(MtCopyNumGcContent, CaseInsensitive) {
    EXPECT_DOUBLE_EQ(MtCopyNumber::compute_gc_content("gcgc"),     1.0);
    EXPECT_DOUBLE_EQ(MtCopyNumber::compute_gc_content("aGcTnN"),   0.5);
    EXPECT_DOUBLE_EQ(MtCopyNumber::compute_gc_content("ATGCatgc"), 0.5);
}

TEST(MtCopyNumGcContent, EmptyAndAllN) {
    EXPECT_DOUBLE_EQ(MtCopyNumber::compute_gc_content(""),       0.0);
    EXPECT_DOUBLE_EQ(MtCopyNumber::compute_gc_content("NnNnNn"), 0.0);
}

// =====================================================================
//  MtCopyNumber::compute_statistics
// =====================================================================

TEST(MtCopyNumStatistics, EmptyVector) {
    auto s = MtCopyNumber::compute_statistics({});
    EXPECT_DOUBLE_EQ(s.mean,     0.0);
    EXPECT_DOUBLE_EQ(s.ci_lower, 0.0);
    EXPECT_DOUBLE_EQ(s.ci_upper, 0.0);
}

TEST(MtCopyNumStatistics, SingleValueCollapsesToMean) {
    auto s = MtCopyNumber::compute_statistics({4.2});
    EXPECT_DOUBLE_EQ(s.mean,     4.2);
    EXPECT_DOUBLE_EQ(s.ci_lower, 4.2);
    EXPECT_DOUBLE_EQ(s.ci_upper, 4.2);
}

TEST(MtCopyNumStatistics, ConstantValuesGiveZeroWidthInterval) {
    // All identical values => stddev = 0 => CI collapses to mean.
    auto s = MtCopyNumber::compute_statistics({5.0, 5.0, 5.0, 5.0});
    EXPECT_DOUBLE_EQ(s.mean,     5.0);
    EXPECT_NEAR     (s.ci_lower, 5.0, 1e-12);
    EXPECT_NEAR     (s.ci_upper, 5.0, 1e-12);
}

TEST(MtCopyNumStatistics, SymmetricIntervalAroundMean) {
    std::vector<double> v{1.0, 2.0, 3.0, 4.0, 5.0};
    auto s = MtCopyNumber::compute_statistics(v);
    EXPECT_DOUBLE_EQ(s.mean, 3.0);
    EXPECT_LT       (s.ci_lower, s.mean);
    EXPECT_GT       (s.ci_upper, s.mean);
    // Interval is centered on the mean (95% CI = mean +/- 1.96 * SE).
    EXPECT_NEAR     (s.mean - s.ci_lower, s.ci_upper - s.mean, 1e-12);
    // Hand-calculated: stddev = sqrt(2.5) ~= 1.5811, SE = stddev/sqrt(5).
    double se     = std::sqrt(2.5) / std::sqrt(5.0);
    double half   = 1.96 * se;
    EXPECT_NEAR(s.ci_lower, 3.0 - half, 1e-9);
    EXPECT_NEAR(s.ci_upper, 3.0 + half, 1e-9);
}

// =====================================================================
//  MtCopyNumber::compute_normalized_ratios
// =====================================================================

namespace {

// Build a minimal, deterministic genome: two autosomes + chrM.
std::vector<MtCopyNumber::ChromosomeData> make_diploid_with_mt() {
    using CD = MtCopyNumber::ChromosomeData;
    std::vector<CD> v;
    v.emplace_back("chr1",  /*length=*/100000);
    v.emplace_back("chr2",  /*length=*/100000);
    v.emplace_back("chrM",  /*length=*/16569);
    return v;
}

}  // namespace

TEST(MtCopyNumNormalize, EmptyInputThrows) {
    std::vector<MtCopyNumber::ChromosomeData> empty;
    EXPECT_THROW(MtCopyNumber::compute_normalized_ratios(empty), std::runtime_error);
}

TEST(MtCopyNumNormalize, ZeroFragmentsThrows) {
    auto v = make_diploid_with_mt();
    // All counts default to 0 -> total_fragments == 0 -> should throw.
    EXPECT_THROW(MtCopyNumber::compute_normalized_ratios(v), std::runtime_error);
}

TEST(MtCopyNumNormalize, BalancedAutosomesGiveCnEqualOne) {
    auto v = make_diploid_with_mt();
    // Same length-normalized fragment density on both autosomes; chrM is empty.
    v[0].count = 100000;
    v[1].count = 100000;
    v[2].count =      0;
    MtCopyNumber::compute_normalized_ratios(v);

    // Length-normalized ratio on autosomes should be equal -> CN ~= 1.
    EXPECT_NEAR(v[0].cn_stats.mean, 1.0, 1e-9);
    EXPECT_NEAR(v[1].cn_stats.mean, 1.0, 1e-9);
    // chrM has zero fragments -> CN mean is 0.
    EXPECT_DOUBLE_EQ(v[2].cn_stats.mean, 0.0);
}

TEST(MtCopyNumNormalize, MitochondrialWeightedByTwo) {
    auto v = make_diploid_with_mt();
    // Make autosomes balanced, give chrM the same length-normalized density
    // as the autosomes => chrM raw normalized_ratio == 1, and the reported
    // CN should be 2 * 1 / 1 = 2 ("copies per diploid cell" baseline).
    v[0].count = 100000;
    v[1].count = 100000;
    v[2].count =  16569;
    MtCopyNumber::compute_normalized_ratios(v);

    EXPECT_NEAR(v[2].normalized_ratio, 1.0, 1e-9);
    EXPECT_NEAR(v[2].cn_stats.mean,    2.0, 1e-9);
}

TEST(MtCopyNumNormalize, ImbalancedAutosomesProduceWiderCi) {
    auto v = make_diploid_with_mt();
    v[0].count =  50000;   // half the density of chr2
    v[1].count = 100000;
    v[2].count =  16569 * 100;  // ~100 mtDNA copies per diploid cell

    MtCopyNumber::compute_normalized_ratios(v);

    // chrM CN > 0 and CI lower < mean < CI upper.
    EXPECT_GT(v[2].cn_stats.mean,    0.0);
    EXPECT_LT(v[2].cn_stats.ci_lower, v[2].cn_stats.mean);
    EXPECT_GT(v[2].cn_stats.ci_upper, v[2].cn_stats.mean);

    // Decoy / non-autosomal contigs are *not* used as the baseline reference,
    // so adding one shouldn't bias chrM's CN value.
    auto v2 = v;
    v2.emplace_back("NUMT_decoy", 50000);
    v2.back().count = 0;
    MtCopyNumber::compute_normalized_ratios(v2);
    EXPECT_NEAR(v[2].cn_stats.mean, v2[2].cn_stats.mean, 1e-9);
}

TEST(MtCopyNumNormalize, NoAutosomesProducesZeroCn) {
    using CD = MtCopyNumber::ChromosomeData;
    std::vector<CD> v;
    v.emplace_back("chrM",       16569);
    v.emplace_back("NUMT_decoy", 64266);
    v[0].count = 1000;
    v[1].count =   10;
    // No autosome anywhere -> nothing to normalize against -> CN stays 0.
    MtCopyNumber::compute_normalized_ratios(v);
    EXPECT_DOUBLE_EQ(v[0].cn_stats.mean,     0.0);
    EXPECT_DOUBLE_EQ(v[0].cn_stats.ci_lower, 0.0);
    EXPECT_DOUBLE_EQ(v[0].cn_stats.ci_upper, 0.0);
    // But normalized_ratio (length-normalized fragment density) IS computed.
    EXPECT_GT(v[0].normalized_ratio, 0.0);
}

// =====================================================================
//  End-to-end run() against the bundled CRAM (chrM + NUMT decoy).
// =====================================================================

namespace {

// Helper: parse a TSV produced by MtCopyNumber and return the row for one
// chromosome (or std::nullopt-equivalent: an empty vector). Skips comment
// lines starting with '#'.
struct TsvRow {
    std::string chrom;
    int64_t     fragments      = 0;
    uint32_t    chrom_length   = 0;
    double      gc_content     = 0.0;
    double      normalized     = 0.0;
    double      copy_number    = 0.0;
    double      ci_lower       = 0.0;
    double      ci_upper       = 0.0;
};

std::vector<TsvRow> parse_tsv(const std::string &path) {
    std::vector<TsvRow> rows;
    std::ifstream ifs(path);
    if (!ifs.is_open()) return rows;
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        TsvRow r;
        ss >> r.chrom >> r.fragments >> r.chrom_length
           >> r.gc_content >> r.normalized
           >> r.copy_number >> r.ci_lower >> r.ci_upper;
        if (!ss.fail()) rows.push_back(std::move(r));
    }
    return rows;
}

}  // namespace

TEST(MtCopyNumRunE2E, SmallBamSampleProducesValidTsv) {
    MtCopyNumber::Config cfg;
    cfg.reference_file = "data/chrM_rCRS.decoy.fa.gz";
    cfg.input_bam      = "data/smp2_subset.bam";
    cfg.output_file    = "data/copynum_smp2.tsv";
    cfg.min_mapq       = 20;
    cfg.thread_count   = 2;
    cfg.seq_type       = SeqType::AUTO;

    // Wipe any stale output from a previous run so a no-op write is detectable.
    std::remove(cfg.output_file.c_str());

    MtCopyNumber tool(cfg);
    ASSERT_NO_THROW(tool.run());

    auto rows = parse_tsv(cfg.output_file);
    ASSERT_FALSE(rows.empty()) << "TSV output not produced or unreadable";

    // Both contigs from the test BAM should be present.
    auto find_row = [&rows](const std::string &name) {
        return std::find_if(rows.begin(), rows.end(),
                            [&](const TsvRow &r){ return r.chrom == name; });
    };
    auto it_mt    = find_row("chrM");
    auto it_decoy = find_row("NUMT_JoinedSequences_gaps1000N_decoy");
    ASSERT_NE(it_mt,    rows.end()) << "chrM row missing from TSV";
    ASSERT_NE(it_decoy, rows.end()) << "decoy row missing from TSV";

    // chrM expectations: subset has thousands of fragments, GC ~ 0.44, length 16569.
    EXPECT_EQ(it_mt->chrom_length, 16569u);
    EXPECT_GT(it_mt->fragments,    100);
    EXPECT_NEAR(it_mt->gc_content, 0.4436, 0.01);

    // The test reference has no autosomes, so CN should be zero on every row
    // but normalized_ratio must be positive for chrM (since it has fragments).
    EXPECT_GT(it_mt->normalized, 0.0);
    EXPECT_DOUBLE_EQ(it_mt->copy_number, 0.0);
    EXPECT_DOUBLE_EQ(it_mt->ci_lower,    0.0);
    EXPECT_DOUBLE_EQ(it_mt->ci_upper,    0.0);

    std::remove(cfg.output_file.c_str());
}

TEST(MtCopyNumRunE2E, RejectsMissingReference) {
    MtCopyNumber::Config cfg;
    cfg.reference_file = "data/__does_not_exist__.fa.gz";
    cfg.input_bam      = "data/smp2_subset.bam";
    cfg.output_file    = "data/copynum_should_not_exist.tsv";
    cfg.min_mapq       = 20;
    cfg.thread_count   = 1;
    cfg.seq_type       = SeqType::PE;

    std::remove(cfg.output_file.c_str());
    MtCopyNumber tool(cfg);
    EXPECT_ANY_THROW(tool.run());
    std::remove(cfg.output_file.c_str());
}

TEST(MtCopyNumRunE2E, RejectsMissingInputBam) {
    MtCopyNumber::Config cfg;
    cfg.reference_file = "data/chrM_rCRS.decoy.fa.gz";
    cfg.input_bam      = "data/__does_not_exist__.bam";
    cfg.output_file    = "data/copynum_should_not_exist.tsv";
    cfg.min_mapq       = 20;
    cfg.thread_count   = 1;
    cfg.seq_type       = SeqType::PE;

    std::remove(cfg.output_file.c_str());
    MtCopyNumber tool(cfg);
    EXPECT_ANY_THROW(tool.run());
    std::remove(cfg.output_file.c_str());
}

TEST(MtCopyNumRunE2E, ManualSeqTypeMatchesAutoDetection) {
    // Force PE explicitly and compare to AUTO. The test BAM is paired-end,
    // so the per-contig fragment counts must agree.
    auto run_with = [](SeqType st, const std::string &out_path) {
        MtCopyNumber::Config cfg;
        cfg.reference_file = "data/chrM_rCRS.decoy.fa.gz";
        cfg.input_bam      = "data/smp2_subset.bam";
        cfg.output_file    = out_path;
        cfg.min_mapq       = 20;
        cfg.thread_count   = 2;
        cfg.seq_type       = st;
        std::remove(out_path.c_str());
        MtCopyNumber tool(cfg);
        tool.run();
        return parse_tsv(out_path);
    };

    auto rows_auto = run_with(SeqType::AUTO, "data/copynum_auto.tsv");
    auto rows_pe   = run_with(SeqType::PE,   "data/copynum_pe.tsv");
    ASSERT_FALSE(rows_auto.empty());
    ASSERT_EQ(rows_auto.size(), rows_pe.size());
    for (size_t i = 0; i < rows_auto.size(); ++i) {
        EXPECT_EQ(rows_auto[i].chrom,        rows_pe[i].chrom);
        EXPECT_EQ(rows_auto[i].fragments,    rows_pe[i].fragments);
        EXPECT_EQ(rows_auto[i].chrom_length, rows_pe[i].chrom_length);
    }

    std::remove("data/copynum_auto.tsv");
    std::remove("data/copynum_pe.tsv");
}
