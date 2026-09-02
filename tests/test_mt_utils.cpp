#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

#include "mt_utils.h"

namespace {

VariantInfo make_variant_info(const std::vector<std::string>& alleles,
                              const std::vector<int>& depths,
                              size_t major_index,
                              int total_depth) {
    VariantInfo info("chrM", 100, total_depth, 60);
    info.major_allele_idx = major_index;
    info.ref_bases.assign(alleles.size(), "A");
    info.alt_bases = alleles;
    info.depths = depths;
    info.var_types.assign(alleles.size(), "SNV");
    info.ci.assign(alleles.size(), std::make_pair(0.0, 1.0));
    info.strand_bias.assign(alleles.size(), StrandBiasInfo{0, 0, 0.0, 0.0});
    return info;
}

}  // namespace

TEST(MtUtilsAlleleQuality, UsesExactBinomialUpperTailForEveryAllele) {
    VariantInfo info = make_variant_info({"A", "G"}, {10, 2}, 0, 10);
    const auto annotation = process_sample_variant(info, {"A", "G"}, 0.2);
    ASSERT_EQ(annotation.aq.size(), 2u);

    // P[X >= 10 | X~Binomial(10, 0.2)] = 0.2^10.
    EXPECT_EQ(annotation.aq[0], 69);
    // P[X >= 2 | X~Binomial(10, 0.2)] = 0.6241903616.
    EXPECT_EQ(annotation.aq[1], 2);
}

TEST(MtUtilsStrandBias, MatchesGatkSymmetricOddsRatioAndIsStrandInvariant) {
    std::vector<std::string> bases;
    std::vector<char> strands;
    auto add = [&](const std::string& base, char strand, int count) {
        for (int index = 0; index < count; ++index) {
            bases.push_back(base);
            strands.push_back(strand);
        }
    };
    add("A", '+', 40);
    add("A", '-', 45);
    add("G", '+', 13);
    add("G", '-', 2);

    const StrandBiasInfo forward = strand_bias("A", "G", bases, strands);
    const double major_forward = 41.0;
    const double major_reverse = 46.0;
    const double alt_forward = 14.0;
    const double alt_reverse = 3.0;
    const double odds = (major_forward / major_reverse) * (alt_reverse / alt_forward);
    const double expected = std::log(odds + 1.0 / odds)
        + std::log(std::min(major_forward, major_reverse) / std::max(major_forward, major_reverse))
        - std::log(std::min(alt_forward, alt_reverse) / std::max(alt_forward, alt_reverse));
    EXPECT_NEAR(forward.sor, expected, 1e-12);

    for (char& strand : strands) strand = (strand == '+') ? '-' : '+';
    const StrandBiasInfo reversed = strand_bias("A", "G", bases, strands);
    EXPECT_NEAR(reversed.sor, forward.sor, 1e-12);
}

TEST(MtUtilsStrandBias, MonomorphicBalanceUsesBinomialTestAndSymmetricLogRatio) {
    // 10 reads all on '+': absolute strand-balance test vs the 50:50
    // expectation.  p = 2 * P[X >= 10 | Bin(10, 0.5)] = 2 * 0.5^10.
    std::vector<std::string> bases(10, "A");
    std::vector<char> strands(10, '+');
    const StrandBiasInfo bias = strand_bias("A", "A", bases, strands);
    const double p_two = 2.0 * std::pow(0.5, 10);
    EXPECT_NEAR(bias.fs, -10.0 * std::log10(p_two), 1e-9);
    // x = (10+1)/(0+1) = 11 -> SOR = ln(11 + 1/11).
    EXPECT_NEAR(bias.sor, std::log(11.0 + 1.0 / 11.0), 1e-12);
    EXPECT_EQ(bias.fwd, 10);
    EXPECT_EQ(bias.rev, 0);

    // Strand-flip invariance: 0 forward / 10 reverse gives the same metrics.
    std::vector<char> rev_strands(10, '-');
    const StrandBiasInfo flipped = strand_bias("A", "A", bases, rev_strands);
    EXPECT_NEAR(flipped.fs, bias.fs, 1e-12);
    EXPECT_NEAR(flipped.sor, bias.sor, 1e-12);

    // Balanced monomorphic: p = 1 -> FS = 0, SOR = ln(2).
    std::vector<std::string> balanced_bases;
    std::vector<char> balanced_strands;
    for (int i = 0; i < 5; ++i) { balanced_bases.push_back("A"); balanced_strands.push_back('+'); }
    for (int i = 0; i < 5; ++i) { balanced_bases.push_back("A"); balanced_strands.push_back('-'); }
    const StrandBiasInfo balanced = strand_bias("A", "A", balanced_bases, balanced_strands);
    EXPECT_NEAR(balanced.fs, 0.0, 1e-12);
    EXPECT_NEAR(balanced.sor, std::log(2.0), 1e-12);

    // The output layer now carries values (no longer missing) for the sole
    // allele, while the boundary logit stays missing.
    VariantInfo info = make_variant_info({"A"}, {10}, 0, 10);
    info.strand_bias[0] = bias;
    const auto annotation = process_sample_variant(info, {"A"}, 0.2);
    ASSERT_EQ(annotation.fs_strings.size(), 1u);
    ASSERT_EQ(annotation.sor_strings.size(), 1u);
    EXPECT_EQ(annotation.fs_strings[0], "27.093");
    EXPECT_EQ(annotation.sor_strings[0], "2.406");
    ASSERT_EQ(annotation.logit_af.size(), 1u);
    EXPECT_EQ(annotation.logit_af[0], ".");
}