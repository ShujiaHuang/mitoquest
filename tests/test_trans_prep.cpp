// Author: Shujia Huang
// Date: 2026-05-28
//
// Unit tests for the `mitoquest trans-prep` subcommand
// (src/trans_prep.{h,cpp}).
//
// Coverage:
//   * FAM parsing + VCF-sample matching statistics
//   * End-to-end run() against a tiny in-test VCF + FAM fixture
//   * QC gating (FILTER=PASS, min-depth threshold, SNV-only filter)

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "io/vcf.h"
#include "io/vcf_header.h"
#include "trans_prep.h"

namespace {

// Small helper: write `content` to a fresh path and return the path.
std::string write_tmp(const std::string& path, const std::string& content) {
    std::ofstream f(path);
    f << content;
    f.close();
    return path;
}

// Minimal multi-sample VCF with three samples (mom_a, child_a, unrelated).
// FORMAT/AD uses the GT-aligned per-sample layout that `mitoquest caller`
// actually emits (AD[i] is the depth of the allele at GT position i).
const char* TEST_VCF = R"(##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allele depths in GT order">
##contig=<ID=chrM,length=16569>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	mom_a	child_a	unrelated
chrM	100	.	A	G	30	PASS	.	GT:DP:AD	0/1:1000:200,800	0/1:900:300,600	0:1100:1100
chrM	200	.	C	T	30	q10	.	GT:DP:AD	0/1:1000:100,900	0/1:900:50,850	0:1100:1100
chrM	300	.	G	A	30	PASS	.	GT:DP:AD	0/1:600:100,500	0/1:50:5,45	0:1100:1100
chrM	400	.	A	AT	30	PASS	.	GT:DP:AD	0/1:1000:200,800	0/1:900:300,600	0:1100:1100
)";

// Two valid trios in the FAM, but the second trio's child is missing
// from the VCF.  PLINK 6 columns: fam child father mother sex pheno.
const char* TEST_FAM =
    "FAM1 child_a 0 mom_a 2 -9\n"
    "FAM2 child_b 0 mom_b 2 -9\n";

}  // namespace

// ---------------------------------------------------------------------
// FAM parsing
// ---------------------------------------------------------------------

TEST(TransPrepParseFam, MatchesAndStats) {
    const std::string vcf_path = "tp_parsefam.vcf";
    const std::string fam_path = "tp_parsefam.fam";
    write_tmp(vcf_path, TEST_VCF);
    write_tmp(fam_path, TEST_FAM);

    ngslib::VCFFile reader(vcf_path);
    TransmissionPrep::MatchingStats stats;
    auto trios = TransmissionPrep::parse_fam(fam_path, reader.header(), stats);

    // FAM2 has child_b/mom_b; neither is in the VCF -> dropped.
    ASSERT_EQ(trios.size(), 1u);
    EXPECT_EQ(trios[0].fam_id,    "FAM1");
    EXPECT_EQ(trios[0].child_id,  "child_a");
    EXPECT_EQ(trios[0].mother_id, "mom_a");
    EXPECT_GE(trios[0].child_idx,  0);
    EXPECT_GE(trios[0].mother_idx, 0);

    EXPECT_EQ(stats.total_vcf_samples,        3);
    EXPECT_EQ(stats.total_fam_lines,          2);
    EXPECT_EQ(stats.valid_mother_child_pairs, 1);
    // FAM2 has both child and mother missing -> incremented in BOTH counts.
    EXPECT_EQ(stats.ignored_missing_child,    1);
    EXPECT_EQ(stats.ignored_missing_mother,   1);

    std::remove(vcf_path.c_str());
    std::remove(fam_path.c_str());
}

TEST(TransPrepParseFam, IgnoresNoMother) {
    const std::string vcf_path = "tp_nomother.vcf";
    const std::string fam_path = "tp_nomother.fam";
    write_tmp(vcf_path, TEST_VCF);
    write_tmp(fam_path,
              "FAM1 child_a dad_a 0 1 -9\n"   // mother encoded as 0 -> ignored
              "FAM2 child_a 0 mom_a 2 -9\n"); // valid

    ngslib::VCFFile reader(vcf_path);
    TransmissionPrep::MatchingStats stats;
    auto trios = TransmissionPrep::parse_fam(fam_path, reader.header(), stats);

    EXPECT_EQ(trios.size(), 1u);
    EXPECT_EQ(stats.ignored_no_mother_in_fam, 1);
    EXPECT_EQ(stats.valid_mother_child_pairs, 1);

    std::remove(vcf_path.c_str());
    std::remove(fam_path.c_str());
}

// ---------------------------------------------------------------------
// End-to-end run()
// ---------------------------------------------------------------------

namespace {

// Read a TSV file written by trans-prep.run() and skip provenance/header.
struct TsvRow {
    std::string chrom, alt, qc;
    std::string grandmother_id, mother_id, child_id;
    int pos, g_dp, g_ad_alt, m_dp, c_dp, m_ad_alt, c_ad_alt, has_g;
};

// Build a header index: column name -> 0-based index.
std::vector<std::string> split_tabs(const std::string& s) {
    std::vector<std::string> tk;
    std::string cur;
    for (char c : s) {
        if (c == '\t') { tk.push_back(cur); cur.clear(); }
        else cur.push_back(c);
    }
    tk.push_back(cur);
    return tk;
}

std::vector<TsvRow> read_pairs_tsv(const std::string& path) {
    std::ifstream f(path);
    std::string line;
    std::unordered_map<std::string, int> idx;
    std::vector<TsvRow> rows;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;          // provenance
        std::vector<std::string> tk = split_tabs(line);
        if (idx.empty()) {
            for (int i = 0; i < (int)tk.size(); ++i) idx[tk[i]] = i;
            continue;
        }
        TsvRow r;
        r.chrom          = tk[idx.at("CHROM")];
        r.pos            = std::stoi(tk[idx.at("POS")]);
        r.alt            = tk[idx.at("ALT")];
        r.grandmother_id = tk[idx.at("GRANDMOTHER_ID")];
        r.mother_id      = tk[idx.at("MOTHER_ID")];
        r.child_id       = tk[idx.at("CHILD_ID")];
        r.g_dp           = std::stoi(tk[idx.at("GRANDMOTHER_DP")]);
        r.g_ad_alt       = std::stoi(tk[idx.at("GRANDMOTHER_AD_ALT")]);
        r.has_g          = std::stoi(tk[idx.at("HAS_G")]);
        r.m_dp           = std::stoi(tk[idx.at("MOTHER_DP")]);
        r.m_ad_alt       = std::stoi(tk[idx.at("MOTHER_AD_ALT")]);
        r.c_dp           = std::stoi(tk[idx.at("CHILD_DP")]);
        r.c_ad_alt       = std::stoi(tk[idx.at("CHILD_AD_ALT")]);
        r.qc             = tk[idx.at("QC")];
        rows.push_back(r);
    }
    return rows;
}

}  // namespace

TEST(TransPrepRun, EndToEndPassFiltering) {
    const std::string vcf_path = "tp_e2e.vcf";
    const std::string fam_path = "tp_e2e.fam";
    const std::string out_path = "tp_e2e.tsv";
    write_tmp(vcf_path, TEST_VCF);
    write_tmp(fam_path, TEST_FAM);

    TransmissionPrep::Config cfg;
    cfg.vcf_path     = vcf_path;
    cfg.fam_path     = fam_path;
    cfg.output_file  = out_path;
    cfg.min_depth    = 500;
    cfg.require_pass = true;
    cfg.snv_only     = true;

    TransmissionPrep tp(cfg);
    long long pass_rows = tp.run();

    auto rows = read_pairs_tsv(out_path);

    // The four input variants:
    //   POS=100 PASS SNV  -> child DP=900>=500, PASS
    //   POS=200 q10       -> dropped by FILTER
    //   POS=300 PASS SNV  -> child DP=50<500    -> LOW_DEPTH (still emitted, not PASS)
    //   POS=400 PASS INS  -> dropped by --snv-only
    //
    // -> 1 PASS row, 1 LOW_DEPTH row.
    EXPECT_EQ(pass_rows, 1);
    ASSERT_EQ(rows.size(), 2u);

    // Locate the rows by position.
    const TsvRow* r100 = nullptr;
    const TsvRow* r300 = nullptr;
    for (const auto& r : rows) {
        if (r.pos == 100) r100 = &r;
        if (r.pos == 300) r300 = &r;
    }
    ASSERT_NE(r100, nullptr);
    ASSERT_NE(r300, nullptr);

    EXPECT_EQ(r100->qc,       "PASS");
    EXPECT_EQ(r100->alt,      "G");
    EXPECT_EQ(r100->m_dp,     1000);
    EXPECT_EQ(r100->m_ad_alt, 800);
    EXPECT_EQ(r100->c_dp,     900);
    EXPECT_EQ(r100->c_ad_alt, 600);

    EXPECT_EQ(r300->qc,       "LOW_DEPTH");

    std::remove(vcf_path.c_str());
    std::remove(fam_path.c_str());
    std::remove(out_path.c_str());
}

TEST(TransPrepRun, NoRequirePassKeepsFiltered) {
    const std::string vcf_path = "tp_nopass.vcf";
    const std::string fam_path = "tp_nopass.fam";
    const std::string out_path = "tp_nopass.tsv";
    write_tmp(vcf_path, TEST_VCF);
    write_tmp(fam_path, TEST_FAM);

    TransmissionPrep::Config cfg;
    cfg.vcf_path     = vcf_path;
    cfg.fam_path     = fam_path;
    cfg.output_file  = out_path;
    cfg.min_depth    = 0;        // accept everything depth-wise
    cfg.require_pass = false;    // keep q10 records too
    cfg.snv_only     = true;

    TransmissionPrep tp(cfg);
    tp.run();

    auto rows = read_pairs_tsv(out_path);
    // Three SNV rows: 100, 200, 300 (the indel at 400 is still excluded).
    EXPECT_EQ(rows.size(), 3u);

    std::remove(vcf_path.c_str());
    std::remove(fam_path.c_str());
    std::remove(out_path.c_str());
}

// ---------------------------------------------------------------------
// Multi-allelic site handling
// ---------------------------------------------------------------------

// Two records that exercise multi-allelic logic:
//   POS=500: pure tri-allelic SNV  (A>G,T)            -> 2 rows per trio
//   POS=600: mixed multi-allelic   (A>G,GT) with --snv-only on
//                                                    -> only the SNV ALT (G)
//                                                       survives, 1 row per trio
// AD is GT-aligned per sample (matches `mitoquest caller` semantics).
const char* TEST_VCF_MULTI = R"(##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allele depths in GT order">
##contig=<ID=chrM,length=16569>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	mom_a	child_a	unrelated
chrM	500	.	A	G,T	30	PASS	.	GT:DP:AD	0/1/2:1000:100,600,300	0/1/2:900:50,500,350	0:1100:1100
chrM	600	.	A	G,GT	30	PASS	.	GT:DP:AD	0/1/2:1000:200,700,100	0/1/2:900:150,650,100	0:1100:1100
)";

TEST(TransPrepRun, HandlesMultiAllelicSites) {
    const std::string vcf_path = "tp_multi.vcf";
    const std::string fam_path = "tp_multi.fam";
    const std::string out_path = "tp_multi.tsv";
    write_tmp(vcf_path, TEST_VCF_MULTI);
    write_tmp(fam_path, TEST_FAM);

    TransmissionPrep::Config cfg;
    cfg.vcf_path     = vcf_path;
    cfg.fam_path     = fam_path;
    cfg.output_file  = out_path;
    cfg.min_depth    = 0;
    cfg.require_pass = true;
    cfg.snv_only     = true;

    TransmissionPrep tp(cfg);
    long long pass_rows = tp.run();

    auto rows = read_pairs_tsv(out_path);

    // POS=500: 2 SNV alts (G,T) x 1 trio  = 2 rows
    // POS=600: 1 SNV alt (G) kept, indel (GT) skipped, x 1 trio = 1 row
    // -> total 3 PASS rows.
    EXPECT_EQ(pass_rows, 3);
    ASSERT_EQ(rows.size(), 3u);

    // Find the row for POS=500, ALT=T (the second alt of the tri-allelic).
    const TsvRow* r500_T = nullptr;
    const TsvRow* r500_G = nullptr;
    const TsvRow* r600_G = nullptr;
    for (const auto& r : rows) {
        if (r.pos == 500 && r.alt == "G") r500_G = &r;
        if (r.pos == 500 && r.alt == "T") r500_T = &r;
        if (r.pos == 600 && r.alt == "G") r600_G = &r;
    }
    ASSERT_NE(r500_G, nullptr);
    ASSERT_NE(r500_T, nullptr);
    ASSERT_NE(r600_G, nullptr);

    // POS=500 mom AD = 100,600,300 ; child AD = 50,500,350
    EXPECT_EQ(r500_G->m_ad_alt, 600);
    EXPECT_EQ(r500_G->c_ad_alt, 500);
    EXPECT_EQ(r500_T->m_ad_alt, 300);
    EXPECT_EQ(r500_T->c_ad_alt, 350);

    // POS=600 mom AD = 200,700,100 (kept ALT is G with index 1)
    EXPECT_EQ(r600_G->m_ad_alt, 700);
    EXPECT_EQ(r600_G->c_ad_alt, 650);

    // No row should be emitted for ALT=GT at POS=600.
    for (const auto& r : rows) {
        EXPECT_FALSE(r.pos == 600 && r.alt == "GT");
    }

    std::remove(vcf_path.c_str());
    std::remove(fam_path.c_str());
    std::remove(out_path.c_str());
}

// ---------------------------------------------------------------------
// GT-aligned AD: missing GT and truncated AD handling
// ---------------------------------------------------------------------

// VCF where some pairs have GT='.' at a site (the caller emits the
// truncated form `.:GQ:DP` under its full FORMAT; here, since our test
// FORMAT is just `GT:DP:AD`, we use `.:DP:.` — GT missing, DP positive,
// AD missing).  trans-prep cannot map AD positions to alleles when GT
// is missing, so the (mother, child) pair is dropped at that site and
// counted in `pair_site_dropped_gt_missing`.
const char* TEST_VCF_GT_MISSING = R"(##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allele depths in GT order">
##contig=<ID=chrM,length=16569>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	mom_a	child_a	unrelated
chrM	700	.	A	G	30	PASS	.	GT:DP:AD	.:1000:.	0/1:900:300,600	0:1100:1100
chrM	800	.	C	T	30	PASS	.	GT:DP:AD	0/1:1000:100,900	.:900:.	0:1100:1100
chrM	900	.	G	A	30	PASS	.	GT:DP:AD	0/1:1000:200,800	0/1:900:300,600	0:1100:1100
)";

TEST(TransPrepRun, GtMissingAlwaysDroppedAndCounted) {
    const std::string vcf_path = "tp_gt_missing.vcf";
    const std::string fam_path = "tp_gt_missing.fam";
    const std::string out_path = "tp_gt_missing.tsv";
    write_tmp(vcf_path, TEST_VCF_GT_MISSING);
    write_tmp(fam_path, TEST_FAM);

    TransmissionPrep::Config cfg;
    cfg.vcf_path     = vcf_path;
    cfg.fam_path     = fam_path;
    cfg.output_file  = out_path;
    cfg.min_depth    = 0;
    cfg.require_pass = true;
    cfg.snv_only     = true;

    TransmissionPrep tp(cfg);
    long long pass_rows = tp.run();

    auto rows = read_pairs_tsv(out_path);
    // POS=700 (mom GT='.') and POS=800 (child GT='.') are dropped; POS=900
    // is kept.  The dropped-pair counter must be 2.
    EXPECT_EQ(pass_rows, 1);
    ASSERT_EQ(rows.size(), 1u);
    EXPECT_EQ(rows[0].pos, 900);
    EXPECT_EQ(rows[0].m_ad_alt, 800);
    EXPECT_EQ(rows[0].c_ad_alt, 600);
    EXPECT_EQ(tp.matching_stats().pair_site_dropped_gt_missing, 2);

    std::remove(vcf_path.c_str());
    std::remove(fam_path.c_str());
    std::remove(out_path.c_str());
}

// VCF where the *child* has GT=0 (homoplasmic REF, AD has 1 value = REF
// depth, NO ALT depth recorded) and the mother is heteroplasmic with
// GT=0/1 and AD=REF,ALT.  trans-prep must report c_ad_alt = 0 (the ALT
// allele is absent from the child's GT) and c_ad_ref = the recorded REF
// depth.  This is the canonical "transmission-loss" event.
const char* TEST_VCF_HOMOPLASMIC_CHILD = R"(##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allele depths in GT order">
##contig=<ID=chrM,length=16569>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	mom_a	child_a	unrelated
chrM	1000	.	A	G	30	PASS	.	GT:DP:AD	0/1:1000:200,800	0:900:900	0:1100:1100
)";

TEST(TransPrepRun, ChildGtZeroYieldsAltCountZero) {
    const std::string vcf_path = "tp_homo_child.vcf";
    const std::string fam_path = "tp_homo_child.fam";
    const std::string out_path = "tp_homo_child.tsv";
    write_tmp(vcf_path, TEST_VCF_HOMOPLASMIC_CHILD);
    write_tmp(fam_path, TEST_FAM);

    TransmissionPrep::Config cfg;
    cfg.vcf_path     = vcf_path;
    cfg.fam_path     = fam_path;
    cfg.output_file  = out_path;
    cfg.min_depth    = 0;
    cfg.require_pass = true;
    cfg.snv_only     = true;

    TransmissionPrep tp(cfg);
    tp.run();

    auto rows = read_pairs_tsv(out_path);
    ASSERT_EQ(rows.size(), 1u);
    EXPECT_EQ(rows[0].pos,      1000);
    // Child's GT=0 ⇒ ALT not in GT ⇒ c_ad_alt must be 0.
    EXPECT_EQ(rows[0].c_ad_alt, 0);
    EXPECT_EQ(rows[0].c_dp,     900);
    // Mother heteroplasmic: GT=0/1, AD=200,800.
    EXPECT_EQ(rows[0].m_ad_alt, 800);

    std::remove(vcf_path.c_str());
    std::remove(fam_path.c_str());
    std::remove(out_path.c_str());
}

TEST(TransPrepFormatRow, Roundtrip) {
    TransmissionPrep::PairRecord r{};
    r.chrom = "chrM"; r.pos = 7;
    r.ref   = "A";   r.alt   = "G";
    r.fam_id = "F1"; r.mother_id = "M1"; r.child_id = "C1";
    r.m_dp = 1000; r.m_ad_ref = 200; r.m_ad_alt = 800; r.m_vaf = 0.8;
    r.c_dp =  900; r.c_ad_ref = 300; r.c_ad_alt = 600; r.c_vaf = 0.666666;
    r.qc   = "PASS";

    const std::string s = TransmissionPrep::format_row(r);
    EXPECT_NE(s.find("chrM"), std::string::npos);
    EXPECT_NE(s.find("F1"),   std::string::npos);
    EXPECT_NE(s.find("PASS"), std::string::npos);
    // Maternal VAF must serialise as "0.8" (no trailing zeros).
    EXPECT_NE(s.find("\t0.8\t"), std::string::npos);
    // When grandmother_id is empty, the row must write NA and HAS_G = 0.
    EXPECT_NE(s.find("\tNA\t"), std::string::npos);
    EXPECT_NE(s.find("\t0\t"), std::string::npos);   // HAS_G column
}

// ---------------------------------------------------------------------
// 3-generation (grandmother-mother-child) extension
// ---------------------------------------------------------------------

// VCF with 4 samples: gma_a (grandmother), mom_a (mother), child_a (child),
// and unrelated.  Two heterozygous SNV sites for all three generations.
const char* TEST_VCF_3GEN = R"(##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allele depths in GT order">
##contig=<ID=chrM,length=16569>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	gma_a	mom_a	child_a	unrelated
chrM	100	.	A	G	30	PASS	.	GT:DP:AD	0/1:1000:200,800	0/1:1000:300,700	0/1:900:400,500	0:1100:1100
chrM	200	.	C	T	30	PASS	.	GT:DP:AD	0/1:1000:250,750	0/1:1000:350,650	0/1:900:450,450	0:1100:1100
)";

// MC FAM: mom_a -> child_a.
const char* TEST_MC_FAM =
    "FAM1 child_a 0 mom_a 2 -9\n";

// GM FAM: gma_a -> mom_a.  The GM FAM's CHILD_ID is the MC-FAM's MOTHER_ID,
// so trans-prep should resolve mom_a as the MC-mother and gma_a as the
// grandmother of child_a.
const char* TEST_GM_FAM =
    "FAM1 mom_a 0 gma_a 2 -9\n";

TEST(TransPrepThreeGen, ResolvesGrandmotherAndEmitsHasGOne) {
    const std::string vcf_path = "tp_3gen.vcf";
    const std::string mc_path  = "tp_3gen_mc.fam";
    const std::string gm_path  = "tp_3gen_gm.fam";
    const std::string out_path = "tp_3gen.tsv";
    write_tmp(vcf_path, TEST_VCF_3GEN);
    write_tmp(mc_path,  TEST_MC_FAM);
    write_tmp(gm_path,  TEST_GM_FAM);

    TransmissionPrep::Config cfg;
    cfg.vcf_path     = vcf_path;
    cfg.fam_path     = mc_path;
    cfg.gm_fam_path  = gm_path;
    cfg.output_file  = out_path;
    cfg.min_depth    = 0;
    cfg.require_pass = true;
    cfg.snv_only     = true;

    TransmissionPrep tp(cfg);
    long long pass_rows = tp.run();

    auto rows = read_pairs_tsv(out_path);
    EXPECT_EQ(pass_rows, 2);
    ASSERT_EQ(rows.size(), 2u);

    for (const auto& r : rows) {
        EXPECT_EQ(r.has_g, 1);
        EXPECT_EQ(r.grandmother_id, "gma_a");
        EXPECT_EQ(r.mother_id,      "mom_a");
        EXPECT_EQ(r.child_id,       "child_a");
    }
    // POS=100: gma AD=(200,800), mom AD=(300,700), child AD=(400,500).
    const TsvRow* r100 = nullptr;
    const TsvRow* r200 = nullptr;
    for (const auto& r : rows) {
        if (r.pos == 100) r100 = &r;
        if (r.pos == 200) r200 = &r;
    }
    ASSERT_NE(r100, nullptr);
    ASSERT_NE(r200, nullptr);
    EXPECT_EQ(r100->g_dp,     1000);
    EXPECT_EQ(r100->g_ad_alt, 800);
    EXPECT_EQ(r100->m_ad_alt, 700);
    EXPECT_EQ(r100->c_ad_alt, 500);
    EXPECT_EQ(r200->g_ad_alt, 750);

    // GM matching stats: 1 MC pair matched to a grandmother.
    EXPECT_EQ(tp.matching_stats().gm_matched_trios,   1);
    EXPECT_EQ(tp.matching_stats().gm_total_fam_lines, 1);

    std::remove(vcf_path.c_str());
    std::remove(mc_path.c_str());
    std::remove(gm_path.c_str());
    std::remove(out_path.c_str());
}

TEST(TransPrepThreeGen, NoGmFamYieldsHasGZero) {
    // Without --gm-fam every row should still have the wide format (so the
    // header includes GRANDMOTHER_* and HAS_G columns) but HAS_G == 0 and
    // GRANDMOTHER_ID == "NA" for every row.
    const std::string vcf_path = "tp_no_gm.vcf";
    const std::string mc_path  = "tp_no_gm_mc.fam";
    const std::string out_path = "tp_no_gm.tsv";
    write_tmp(vcf_path, TEST_VCF_3GEN);
    write_tmp(mc_path,  TEST_MC_FAM);

    TransmissionPrep::Config cfg;
    cfg.vcf_path     = vcf_path;
    cfg.fam_path     = mc_path;
    cfg.gm_fam_path.clear();
    cfg.output_file  = out_path;
    cfg.min_depth    = 0;
    cfg.require_pass = true;
    cfg.snv_only     = true;

    TransmissionPrep tp(cfg);
    tp.run();

    auto rows = read_pairs_tsv(out_path);
    ASSERT_EQ(rows.size(), 2u);
    for (const auto& r : rows) {
        EXPECT_EQ(r.has_g, 0);
        EXPECT_EQ(r.grandmother_id, "NA");
        EXPECT_EQ(r.g_dp, 0);
        EXPECT_EQ(r.g_ad_alt, 0);
    }
    // GM counters must all be zero when --gm-fam was not given.
    EXPECT_EQ(tp.matching_stats().gm_total_fam_lines, 0);
    EXPECT_EQ(tp.matching_stats().gm_matched_trios,   0);

    std::remove(vcf_path.c_str());
    std::remove(mc_path.c_str());
    std::remove(out_path.c_str());
}

TEST(TransPrepThreeGen, GmMotherNotInVcfFallsBackToHasGZero) {
    // GM FAM's MOTHER_ID (the pedigree grandmother) is not in the VCF, so
    // the GM parse drops that line and the MC pair stays HAS_G = 0.
    const char* gm_missing =
        "FAM1 mom_a 0 gma_xxx 2 -9\n";   // gma_xxx not in VCF

    const std::string vcf_path = "tp_gm_miss.vcf";
    const std::string mc_path  = "tp_gm_miss_mc.fam";
    const std::string gm_path  = "tp_gm_miss_gm.fam";
    const std::string out_path = "tp_gm_miss.tsv";
    write_tmp(vcf_path, TEST_VCF_3GEN);
    write_tmp(mc_path,  TEST_MC_FAM);
    write_tmp(gm_path,  gm_missing);

    TransmissionPrep::Config cfg;
    cfg.vcf_path     = vcf_path;
    cfg.fam_path     = mc_path;
    cfg.gm_fam_path  = gm_path;
    cfg.output_file  = out_path;
    cfg.min_depth    = 0;
    cfg.require_pass = true;
    cfg.snv_only     = true;

    TransmissionPrep tp(cfg);
    tp.run();

    auto rows = read_pairs_tsv(out_path);
    ASSERT_EQ(rows.size(), 2u);
    for (const auto& r : rows) {
        EXPECT_EQ(r.has_g, 0);
        EXPECT_EQ(r.grandmother_id, "NA");
    }
    EXPECT_EQ(tp.matching_stats().gm_ignored_missing_mother, 1);
    EXPECT_EQ(tp.matching_stats().gm_matched_trios,          0);

    std::remove(vcf_path.c_str());
    std::remove(mc_path.c_str());
    std::remove(gm_path.c_str());
    std::remove(out_path.c_str());
}

