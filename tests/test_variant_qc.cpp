// Unit tests for the `mitoquest variant-qc` subcommand
// (src/variant_qc.{h,cpp}).
//
// Coverage:
//   * log_beta_fn:   closed-form spot checks.
//   * beta_pdf:      known Beta distribution values.
//   * binomial_ppf:  quantile function correctness.
//   * fit_beta_mle:  parameter recovery from synthetic VAFs.
//   * fit_betabinom_mle: parameter recovery from count data.
//   * bayesian_filter:   posterior boundary and edge cases.
//   * kl_divergence: divergence of known distributions.
//   * is_blacklisted: blacklist region containment.
//   * nelder_mead_2d: minimize a simple quadratic.

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <random>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <numeric>
#include <utility>

#include "variant_qc.h"
#include "log_factorial.h"
#include "io/iobgzf.h"
#include "io/vcf.h"

namespace {

void write_text_file(const std::string& path, const std::string& content) {
    std::ofstream output(path);
    ASSERT_TRUE(output.is_open());
    output << content;
}

std::string read_text_file(const std::string& path) {
    std::ifstream input(path);
    return std::string(std::istreambuf_iterator<char>(input),
                       std::istreambuf_iterator<char>());
}

VariantQC::Config make_variant_qc_config(const std::string& input,
                                         const std::string& output_vcf,
                                         const std::string& output_tsv) {
    VariantQC::Config config{};
    config.input_vcf = input;
    config.output_vcf = output_vcf;
    config.output_tsv = output_tsv;
    config.max_alt_alleles = 2;
    config.dp_threshold = 1;
    config.hq_threshold = 20;
    config.bins = 100;
    config.pi = 5e-8 * 16569;
    config.threshold = 0.9;
    config.max_iter = 1;
    config.convergence_eps = 0.001;
    return config;
}

}  // namespace

// =====================================================================
// log_beta_fn
// =====================================================================

TEST(VariantQCMath, LogBetaFnKnownValues) {
    // B(1, 1) = 1 -> log = 0
    EXPECT_NEAR(VariantQC::log_beta_fn(1.0, 1.0), 0.0, 1e-12);

    // B(2, 3) = Gamma(2)*Gamma(3)/Gamma(5) = 1*2/24 = 1/12
    EXPECT_NEAR(VariantQC::log_beta_fn(2.0, 3.0), std::log(1.0 / 12.0), 1e-9);

    // B(0.5, 0.5) = pi
    EXPECT_NEAR(VariantQC::log_beta_fn(0.5, 0.5), std::log(M_PI), 1e-9);

    // Symmetry: B(a, b) = B(b, a)
    EXPECT_NEAR(VariantQC::log_beta_fn(3.7, 2.1),
                VariantQC::log_beta_fn(2.1, 3.7), 1e-12);
}

// =====================================================================
// beta_pdf
// =====================================================================

TEST(VariantQCMath, BetaPdfUniformPrior) {
    // Beta(1, 1) is Uniform(0, 1): f(x) = 1 for x in (0, 1)
    for (double x : {0.1, 0.3, 0.5, 0.7, 0.9}) {
        EXPECT_NEAR(VariantQC::beta_pdf(x, 1.0, 1.0), 1.0, 1e-9);
    }
}

TEST(VariantQCMath, BetaPdfBoundaryValues) {
    // PDF should be 0 at x=0 and x=1
    EXPECT_NEAR(VariantQC::beta_pdf(0.0, 2.0, 3.0), 0.0, 1e-12);
    EXPECT_NEAR(VariantQC::beta_pdf(1.0, 2.0, 3.0), 0.0, 1e-12);

    // Invalid parameters
    EXPECT_NEAR(VariantQC::beta_pdf(0.5, -1.0, 2.0), 0.0, 1e-12);
    EXPECT_NEAR(VariantQC::beta_pdf(0.5, 2.0, 0.0), 0.0, 1e-12);
}

TEST(VariantQCMath, BetaPdfKnownPoint) {
    // Beta(2, 5) at x=0.3: f(0.3) = 0.3^1 * 0.7^4 / B(2,5)
    // B(2,5) = 1*24/720 = 1/30
    double expected = std::pow(0.3, 1.0) * std::pow(0.7, 4.0) / (1.0 / 30.0);
    EXPECT_NEAR(VariantQC::beta_pdf(0.3, 2.0, 5.0), expected, 1e-6);
}

// =====================================================================
// binomial_ppf
// =====================================================================

TEST(VariantQCMath, BinomialPpfBasic) {
    mitoquest::LogFactorial lf(1000);

    // PPF(0.5) for Binom(10, 0.5) should be 5 (median)
    int ppf50 = VariantQC::binomial_ppf(10, 0.5, 0.5, lf);
    EXPECT_EQ(ppf50, 5);

    // PPF(0.99) for Binom(100, 0.01): should be a small number
    int ppf99 = VariantQC::binomial_ppf(100, 0.01, 0.99, lf);
    EXPECT_GE(ppf99, 0);
    EXPECT_LE(ppf99, 5);  // very low error rate, 99th percentile is small
}

TEST(VariantQCMath, BinomialPpfBoundaries) {
    mitoquest::LogFactorial lf(100);

    // p=0: always 0
    EXPECT_EQ(VariantQC::binomial_ppf(10, 0.0, 0.99, lf), 0);

    // p=1: always n
    EXPECT_EQ(VariantQC::binomial_ppf(10, 1.0, 0.5, lf), 10);

    // q=1: should return n
    EXPECT_EQ(VariantQC::binomial_ppf(10, 0.5, 1.0, lf), 10);

    // n=0
    EXPECT_EQ(VariantQC::binomial_ppf(0, 0.5, 0.5, lf), 0);
}

TEST(VariantQCMath, BinomialPpfMonotonicity) {
    mitoquest::LogFactorial lf(1000);
    // PPF should be monotonically non-decreasing in q
    int prev = 0;
    for (double q : {0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99}) {
        int k = VariantQC::binomial_ppf(100, 0.1, q, lf);
        EXPECT_GE(k, prev) << "q=" << q;
        prev = k;
    }
}

// =====================================================================
// fit_beta_mle
// =====================================================================

TEST(VariantQCFitting, BetaMLERecoversParameters) {
    // Generate samples from Beta(2, 5) and check MLE recovery
    std::mt19937 rng(42);
    std::gamma_distribution<double> ga(2.0, 1.0);
    std::gamma_distribution<double> gb(5.0, 1.0);

    std::vector<double> vafs;
    for (int i = 0; i < 5000; ++i) {
        double a = ga(rng);
        double b = gb(rng);
        vafs.push_back(a / (a + b));
    }

    auto params = VariantQC::fit_beta_mle(vafs);
    // Allow 20% relative error
    EXPECT_NEAR(params.alpha, 2.0, 0.4);
    EXPECT_NEAR(params.beta, 5.0, 1.0);
}

TEST(VariantQCFitting, BetaMLEEmptyInput) {
    auto params = VariantQC::fit_beta_mle({});
    // Should return defaults
    EXPECT_GT(params.alpha, 0.0);
    EXPECT_GT(params.beta, 0.0);
}

TEST(VariantQCFitting, BetaMLEFewValidValues) {
    auto params = VariantQC::fit_beta_mle({0.5});
    EXPECT_GT(params.alpha, 0.0);
    EXPECT_GT(params.beta, 0.0);
}

// =====================================================================
// fit_betabinom_mle
// =====================================================================

TEST(VariantQCFitting, BetaBinomMLERecoversParameters) {
    // Generate data from BetaBinom(n, alpha=2, beta=5)
    std::mt19937 rng(123);
    double true_alpha = 2.0, true_beta = 5.0;

    std::vector<double> vafs;
    std::vector<int> k_obs, n_obs;

    std::gamma_distribution<double> ga(true_alpha, 1.0);
    std::gamma_distribution<double> gb(true_beta, 1.0);

    for (int i = 0; i < 2000; ++i) {
        double a = ga(rng);
        double b = gb(rng);
        double p = a / (a + b);
        int n = 500 + (rng() % 500); // n in [500, 1000)
        std::binomial_distribution<int> binom(n, p);
        int k = binom(rng);
        vafs.push_back(static_cast<double>(k) / n);
        k_obs.push_back(k);
        n_obs.push_back(n);
    }

    mitoquest::LogFactorial lf(2000);
    auto params = VariantQC::fit_betabinom_mle(vafs, k_obs, n_obs, lf);

    // Parameters should be in the right ballpark (30% tolerance)
    EXPECT_NEAR(params.alpha, true_alpha, true_alpha * 0.5);
    EXPECT_NEAR(params.beta, true_beta, true_beta * 0.5);
}

TEST(VariantQCFitting, BetaBinomMLEEmptyInput) {
    mitoquest::LogFactorial lf(100);
    auto params = VariantQC::fit_betabinom_mle({}, {}, {}, lf);
    EXPECT_GT(params.alpha, 0.0);
    EXPECT_GT(params.beta, 0.0);
}

TEST(VariantQCFitting, BetaBinomMLERejectsMismatchedOrInvalidObservations) {
    mitoquest::LogFactorial lf(100);
    EXPECT_THROW(VariantQC::fit_betabinom_mle({0.2}, {}, {}, lf),
                 std::invalid_argument);
    EXPECT_THROW(VariantQC::fit_betabinom_mle({0.2}, {2}, {1}, lf),
                 std::invalid_argument);
    EXPECT_THROW(VariantQC::fit_betabinom_mle(
                     {std::numeric_limits<double>::quiet_NaN()}, {1}, {10}, lf),
                 std::invalid_argument);
}

// =====================================================================
// bayesian_filter
// =====================================================================

TEST(VariantQCFilter, HighAltCountGivesHighPosterior) {
    // With 200 alt reads out of 500 total, this should be a strong mutation
    // signal under any reasonable background model.
    double pp = VariantQC::bayesian_filter(
        /*A=*/200, /*D=*/500,
        /*srf=*/0.9, /*hq=*/30.0, /*hq_threshold=*/20,
        /*q_alpha=*/1.0, /*q_beta=*/100.0,  // background noise: low error
        /*alpha_h1=*/2.0, /*beta_h1=*/3.0,   // true mutation: moderate VAF
        /*pi=*/0.001);

    EXPECT_GT(pp, 0.5);
}

TEST(VariantQCFilter, LowAltCountGivesLowPosterior) {
    // 2 alt reads out of 500 with a low-error background: likely noise
    double pp = VariantQC::bayesian_filter(
        /*A=*/2, /*D=*/500,
        /*srf=*/0.5, /*hq=*/30.0, /*hq_threshold=*/20,
        /*q_alpha=*/1.0, /*q_beta=*/100.0,
        /*alpha_h1=*/2.0, /*beta_h1=*/3.0,
        /*pi=*/0.001);

    EXPECT_LT(pp, 0.5);
}

TEST(VariantQCFilter, EdgeCasesReturnZero) {
    // D=0
    EXPECT_EQ(VariantQC::bayesian_filter(0, 0, 0.5, 30.0, 20, 1.0, 1.0, 1.0, 1.0, 0.001), 0.0);

    // A < 0
    EXPECT_EQ(VariantQC::bayesian_filter(-1, 100, 0.5, 30.0, 20, 1.0, 1.0, 1.0, 1.0, 0.001), 0.0);

    // A > D
    EXPECT_EQ(VariantQC::bayesian_filter(200, 100, 0.5, 30.0, 20, 1.0, 1.0, 1.0, 1.0, 0.001), 0.0);
}

TEST(VariantQCFilter, StrandBiasReducesPosterior) {
    // Same data, but different SRF values
    double pp_balanced = VariantQC::bayesian_filter(
        50, 500, /*srf=*/0.95, 30.0, 20,
        1.0, 50.0, 2.0, 3.0, 0.001);

    double pp_biased = VariantQC::bayesian_filter(
        50, 500, /*srf=*/0.05, 30.0, 20,
        1.0, 50.0, 2.0, 3.0, 0.001);

    // Balanced strand ratio should give higher posterior
    EXPECT_GT(pp_balanced, pp_biased);
}

TEST(VariantQCFilter, LowHqReducesPosterior) {
    double pp_high_hq = VariantQC::bayesian_filter(
        50, 500, 0.9, /*hq=*/30.0, 20,
        1.0, 50.0, 2.0, 3.0, 0.001);

    double pp_low_hq = VariantQC::bayesian_filter(
        50, 500, 0.9, /*hq=*/5.0, 20,
        1.0, 50.0, 2.0, 3.0, 0.001);

    EXPECT_GT(pp_high_hq, pp_low_hq);
}

TEST(VariantQCFilter, PosteriorBoundedIn01) {
    for (int A = 0; A <= 100; A += 10) {
        double pp = VariantQC::bayesian_filter(
            A, 100, 0.5, 25.0, 20,
            2.0, 5.0, 1.0, 2.0, 0.01);
        EXPECT_GE(pp, 0.0);
        EXPECT_LE(pp, 1.0);
    }
}

TEST(VariantQCQual, PhredMatchesPosteriorProbabilitiesAndExtremeLogOdds) {
    for (double p : {0.5, 0.9, 0.99}) {
        SCOPED_TRACE(p);
        const double z = std::log(p) - std::log1p(-p);
        EXPECT_NEAR(VariantQC::phred_from_log_odds(z), -10.0 * std::log10(1.0 - p), 1e-12);
    }
    EXPECT_NEAR(VariantQC::phred_from_log_odds(1000.0), 10000.0 / std::log(10.0), 1e-10);
    EXPECT_DOUBLE_EQ(VariantQC::phred_from_log_odds(-1000.0), 0.0);
    // Negative evidence is not rounded to zero before Phred conversion.
    EXPECT_GT(VariantQC::phred_from_log_odds(-100.0), 0.0);
}

TEST(VariantQCQual, SaturatedPosteriorStillHasFiniteDistinctUncappedQual) {
    double previous = 0.0;
    for (int count : {128, 256, 512}) {
        SCOPED_TRACE(count);
        const double pp = VariantQC::bayesian_filter(
            count, 1024, 1.0, 30.0, 20, 2.5, 1680.0, 0.5, 0.2, 5e-8 * 16569);
        const double z = VariantQC::bayesian_log_odds(
            count, 1024, 1.0, 30.0, 20, 2.5, 1680.0, 0.5, 0.2, 5e-8 * 16569);
        const double q = VariantQC::phred_from_log_odds(z);
        EXPECT_DOUBLE_EQ(pp, 1.0);
        EXPECT_TRUE(std::isfinite(z));
        EXPECT_TRUE(std::isfinite(q));
        EXPECT_GT(q, previous);
        previous = q;
    }
    EXPECT_GT(previous, 1000.0);
}

TEST(VariantQCQual, InvalidCountsAndNonfiniteEvidenceAreMissingRatherThanZeroQual) {
    for (const auto& counts : {std::make_pair(0, 0), std::make_pair(0, -1),
                               std::make_pair(-1, 100), std::make_pair(101, 100)}) {
        const double z = VariantQC::bayesian_log_odds(
            counts.first, counts.second, 1.0, 30.0, 20, 2.5, 1680.0, 0.5, 0.2, 0.001);
        EXPECT_TRUE(std::isnan(z));
        EXPECT_TRUE(std::isnan(VariantQC::phred_from_log_odds(z)));
        EXPECT_DOUBLE_EQ(VariantQC::bayesian_filter(
            counts.first, counts.second, 1.0, 30.0, 20, 2.5, 1680.0, 0.5, 0.2, 0.001), 0.0);
    }
    for (double z : {std::numeric_limits<double>::quiet_NaN(),
                     std::numeric_limits<double>::infinity(),
                     -std::numeric_limits<double>::infinity()}) {
        EXPECT_TRUE(std::isnan(VariantQC::phred_from_log_odds(z)));
    }
    for (double srf : {std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::infinity()}) {
        const double z = VariantQC::bayesian_log_odds(
            50, 100, srf, 30.0, 20, 2.5, 1680.0, 0.5, 0.2, 0.001);
        EXPECT_FALSE(std::isfinite(z));
        EXPECT_TRUE(std::isnan(VariantQC::phred_from_log_odds(z)));
        // Keep the legacy PP nonfinite fallback, even though QUAL is missing.
        EXPECT_DOUBLE_EQ(VariantQC::bayesian_filter(
            50, 100, srf, 30.0, 20, 2.5, 1680.0, 0.5, 0.2, 0.001), z > 0.0 ? 1.0 : 0.0);
    }
}

TEST(VariantQCQual, LogOddsExtractionPreservesLegacyPosteriorArithmetic) {
    for (int count : {0, 2, 50, 100}) {
        for (double hq : {19.0, 30.0}) {
            for (double pi : {0.0, 0.001, 1.0}) {
                SCOPED_TRACE(::testing::Message() << count << ", HQ=" << hq << ", pi=" << pi);
                const double prior = std::max(1e-12, std::min(pi, 1.0 - 1e-12));
                const double penalty = hq < 20 ? 1.0 / (1.0 + std::exp(10.0 * (20 - hq))) : 1.0;
                const double log_l0 = VariantQC::log_beta_fn(count + 2.5, 100 - count + 1680.0) -
                                      VariantQC::log_beta_fn(2.5, 1680.0);
                double log_l1 = VariantQC::log_beta_fn(count + 0.5, 100 - count + 0.2) -
                                VariantQC::log_beta_fn(0.5, 0.2);
                log_l1 += std::log(penalty * 0.5 + 1e-12);
                const double z = (std::log(prior) - std::log1p(-prior)) + log_l1 - log_l0;
                EXPECT_DOUBLE_EQ(VariantQC::bayesian_log_odds(
                    count, 100, 0.5, hq, 20, 2.5, 1680.0, 0.5, 0.2, pi), z);
                EXPECT_DOUBLE_EQ(VariantQC::bayesian_filter(
                    count, 100, 0.5, hq, 20, 2.5, 1680.0, 0.5, 0.2, pi),
                    1.0 / (1.0 + std::exp(-z)));
            }
        }
    }
}

// =====================================================================
// kl_divergence
// =====================================================================

TEST(VariantQCKL, ZeroForEmptyVAFs) {
    std::vector<double> bin_edges(101);
    for (int i = 0; i <= 100; ++i) bin_edges[i] = i / 100.0;
    EXPECT_NEAR(VariantQC::kl_divergence({}, 1.0, 1.0, bin_edges), 0.0, 1e-12);
}

TEST(VariantQCKL, NonNegative) {
    std::mt19937 rng(42);
    std::vector<double> vafs;
    for (int i = 0; i < 200; ++i) {
        vafs.push_back(static_cast<double>(rng() % 1000) / 1000.0);
    }
    std::vector<double> bin_edges(101);
    for (int i = 0; i <= 100; ++i) bin_edges[i] = i / 100.0;

    double kl = VariantQC::kl_divergence(vafs, 2.0, 5.0, bin_edges);
    EXPECT_GE(kl, 0.0);
}

TEST(VariantQCKL, HigherForMismatchedDistribution) {
    // Generate data from Beta(2, 2) and compare KL with Beta(2,2) vs Beta(10,10)
    std::mt19937 rng(42);
    std::gamma_distribution<double> ga(2.0, 1.0);
    std::gamma_distribution<double> gb(2.0, 1.0);

    std::vector<double> vafs;
    for (int i = 0; i < 1000; ++i) {
        double a = ga(rng), b = gb(rng);
        vafs.push_back(a / (a + b));
    }

    std::vector<double> bin_edges(101);
    for (int i = 0; i <= 100; ++i) bin_edges[i] = i / 100.0;

    double kl_match = VariantQC::kl_divergence(vafs, 2.0, 2.0, bin_edges);
    double kl_mismatch = VariantQC::kl_divergence(vafs, 10.0, 10.0, bin_edges);

    EXPECT_LT(kl_match, kl_mismatch);
}

// =====================================================================
// is_blacklisted
// =====================================================================

TEST(VariantQCBlacklist, KnownBlacklistedPositions) {
    VariantQC::init_blacklist();

    // Positions within known blacklisted regions
    EXPECT_TRUE(VariantQC::is_blacklisted(300));   // in (299, 317)
    EXPECT_TRUE(VariantQC::is_blacklisted(317));   // boundary
    EXPECT_TRUE(VariantQC::is_blacklisted(511));   // in (511, 525)
    EXPECT_TRUE(VariantQC::is_blacklisted(16184)); // in (16180, 16187)
}

TEST(VariantQCBlacklist, NonBlacklistedPositions) {
    VariantQC::init_blacklist();

    EXPECT_FALSE(VariantQC::is_blacklisted(1));
    EXPECT_FALSE(VariantQC::is_blacklisted(1000));
    EXPECT_FALSE(VariantQC::is_blacklisted(5000));
    EXPECT_FALSE(VariantQC::is_blacklisted(16569));
}

TEST(VariantQCConfig, RejectsNonPositiveBinCounts) {
    VariantQC::Config config{};
    config.bins = 0;
    EXPECT_THROW((void)VariantQC{config}, std::invalid_argument);
    config.bins = -1;
    EXPECT_THROW((void)VariantQC{config}, std::invalid_argument);
}

TEST(VariantQCConfig, RejectsNonFiniteOrOutOfRangeStatisticalOptions) {
    VariantQC::Config config = make_variant_qc_config("input.vcf", "output.vcf", "output.tsv");
    config.pi = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW((void)VariantQC{config}, std::invalid_argument);

    config.pi = 0.5;
    config.threshold = 1.1;
    EXPECT_THROW((void)VariantQC{config}, std::invalid_argument);

    config.threshold = 0.9;
    config.convergence_eps = std::numeric_limits<double>::infinity();
    EXPECT_THROW((void)VariantQC{config}, std::invalid_argument);
}

TEST(VariantQCPipeline, ReferenceOnlyCleanSitePasses) {
    const std::string input = "/tmp/mitoquest_variant_qc_reference_only.vcf";
    const std::string output_vcf = "/tmp/mitoquest_variant_qc_reference_only_out.vcf";
    const std::string output_tsv = "/tmp/mitoquest_variant_qc_reference_only_out.tsv";
    write_text_file(input, R"(##fileformat=VCFv4.2
##contig=<ID=chrM,length=16569>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Depths">
##FORMAT=<ID=AF,Number=.,Type=Float,Description="Frequencies">
##FORMAT=<ID=AQ,Number=.,Type=Integer,Description="Allele quality">
##FORMAT=<ID=FS,Number=.,Type=Float,Description="Fisher strand">
##FORMAT=<ID=SOR,Number=.,Type=Float,Description="Strand odds ratio">
##FORMAT=<ID=SB,Number=1,Type=String,Description="Strand counts">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ref_only
chrM	100	.	A	G	60	PASS	.	GT:DP:AD:AF:AQ:FS:SOR:SB	0:500:500:1:30:0:0:250,250
)");

    VariantQC quality_control(make_variant_qc_config(input, output_vcf, output_tsv));
    quality_control.run();

    const std::string output = read_text_file(output_vcf);
    // A clean reference-only sample is a confident call (posterior 1.0), so
    // the site passes; it has no real ALT allele, so no TSV row is emitted.
    EXPECT_NE(output.find("\tPASS\t"), std::string::npos);
    EXPECT_EQ(output.find("\tLOW_QUALITY\t"), std::string::npos);
    EXPECT_NE(output.find(":1:True"), std::string::npos);
    EXPECT_EQ(read_text_file(output_tsv).find("ref_only\t"), std::string::npos);

    std::remove(input.c_str());
    std::remove(output_vcf.c_str());
    std::remove(output_tsv.c_str());
}

TEST(VariantQCPipeline, BlacklistSkippedRecordsUsePosOnlyFilter) {
    const std::string input = "/tmp/mitoquest_variant_qc_blacklist_overlap.vcf";
    const std::string output_vcf = "/tmp/mitoquest_variant_qc_blacklist_overlap_out.vcf";
    const std::string output_tsv = "/tmp/mitoquest_variant_qc_blacklist_overlap_out.tsv";
    write_text_file(input, R"(##fileformat=VCFv4.2
##contig=<ID=chrM,length=16569>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Depths">
##FORMAT=<ID=AF,Number=.,Type=Float,Description="Frequencies">
##FORMAT=<ID=AQ,Number=.,Type=Integer,Description="Allele quality">
##FORMAT=<ID=FS,Number=.,Type=Float,Description="Fisher strand">
##FORMAT=<ID=SOR,Number=.,Type=Float,Description="Strand odds ratio">
##FORMAT=<ID=SB,Number=1,Type=String,Description="Strand counts">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample
chrM	298	.	AC	A	60	PASS	.	GT:DP:AD:AF:AQ:FS:SOR:SB	0/1:200:100,100:0.5,0.5:30,30:0,0:0,0:50,50;50,50
chrM	316	.	T	C	60	PASS	.	GT:DP:AD:AF:AQ:FS:SOR:SB	0/1:200:100,100:0.5,0.5:30,30:0,0:0,0:50,50;50,50
)");

    VariantQC quality_control(make_variant_qc_config(input, output_vcf, output_tsv));
    quality_control.run();

    const std::string output = read_text_file(output_vcf);
    // Both records are skipped by the QC loop (298 spans blacklisted 299;
    // 316 is itself blacklisted), so no PP/GOOD_CALL is added. The FILTER
    // uses the pos-only check: 298 -> LOW_QUALITY, 316 -> BLACKLISTED_SITE.
    EXPECT_NE(output.find("\tLOW_QUALITY\t"), std::string::npos);
    EXPECT_NE(output.find("\tBLACKLISTED_SITE\t"), std::string::npos);
    EXPECT_EQ(output.find(":PP"), std::string::npos); // FORMAT untouched
    EXPECT_EQ(output.find(":True"), std::string::npos);
    EXPECT_EQ(output.find(":False"), std::string::npos);
    // No QC results -> no TSV rows.
    EXPECT_EQ(read_text_file(output_tsv).find("\tchrM\t"), std::string::npos);

    std::remove(input.c_str());
    std::remove(output_vcf.c_str());
    std::remove(output_tsv.c_str());
}

TEST(VariantQCPipeline, MissingStrandBiasFieldSuppressesBayesianCall) {
    const std::string input = "/tmp/mitoquest_variant_qc_missing_sb.vcf";
    const std::string output_vcf = "/tmp/mitoquest_variant_qc_missing_sb_out.vcf";
    const std::string output_tsv = "/tmp/mitoquest_variant_qc_missing_sb_out.tsv";
    write_text_file(input, R"(##fileformat=VCFv4.2
##contig=<ID=chrM,length=16569>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Depths">
##FORMAT=<ID=AF,Number=.,Type=Float,Description="Frequencies">
##FORMAT=<ID=AQ,Number=.,Type=Integer,Description="Allele quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ref_sample	alt_sample
chrM	100	.	A	G	60	PASS	.	GT:DP:AD:AF:AQ	0:500:500:1:30	1:500:200:0.4:30
)");

    VariantQC quality_control(make_variant_qc_config(input, output_vcf, output_tsv));
    quality_control.run();

    const std::string vcf_output = read_text_file(output_vcf);
    // With SB absent from the FORMAT header, no sample gets an SRF value, so
    // the Bayesian call is suppressed (posterior 0) for every sample: the
    // site is LOW_QUALITY and no TSV row is emitted.
    EXPECT_NE(vcf_output.find("\tLOW_QUALITY\t"), std::string::npos);
    EXPECT_EQ(vcf_output.find("\tPASS\t"), std::string::npos);
    EXPECT_NE(vcf_output.find(":0:False"), std::string::npos);
    EXPECT_NE(vcf_output.find("\t.\tLOW_QUALITY\tKL="), std::string::npos);
    EXPECT_EQ(read_text_file(output_tsv).find("\tchrM\t"), std::string::npos);

    std::remove(input.c_str());
    std::remove(output_vcf.c_str());
    std::remove(output_tsv.c_str());
}

class VariantQCOutput : public ::testing::TestWithParam<std::string> {};

TEST_P(VariantQCOutput, SuffixControlsCompressionAndPreservesContent) {
    const std::string suffix = GetParam();
    const bool compressed = suffix == ".gz" || suffix == ".bgz";
    const std::string prefix = "/tmp/mitoquest_variant_qc_output" + suffix;
    const std::string input = prefix + "_input.vcf" + suffix;
    const std::string plain_vcf = prefix + "_plain.vcf";
    const std::string plain_tsv = prefix + "_plain.tsv";
    const std::string output_vcf = prefix + "_out.vcf" + suffix;
    const std::string output_tsv = prefix + "_out.tsv" + suffix;
    {
        ngslib::BGZFile file(input, compressed ? "w" : "wu");
        file.write(R"(##fileformat=VCFv4.2
##contig=<ID=chrM,length=16569>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Depths">
##FORMAT=<ID=AF,Number=.,Type=Float,Description="Frequencies">
##FORMAT=<ID=AQ,Number=.,Type=Integer,Description="Allele quality">
##FORMAT=<ID=FS,Number=.,Type=Float,Description="Fisher strand">
##FORMAT=<ID=SOR,Number=.,Type=Float,Description="Strand odds ratio">
##FORMAT=<ID=SB,Number=1,Type=String,Description="Strand counts">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ref_sample	alt_sample
chrM	100	.	A	G	60	PASS	.	GT:DP:AD:AF:AQ:FS:SOR:SB	0:500:500:1:30:0:0:250,250	1:500:200:0.4:30:0:0:100,100
chrM	316	.	T	C	60	PASS	.	GT:DP:AD:AF:AQ:FS:SOR:SB	0:500:500:1:30:0:0:250,250	1:500:200:0.4:30:0:0:100,100
)");
    }

    VariantQC(make_variant_qc_config(input, plain_vcf, plain_tsv)).run();
    VariantQC(make_variant_qc_config(input, output_vcf, output_tsv)).run();

    EXPECT_NE(read_text_file(plain_vcf).find(":PP:GOOD_CALL"), std::string::npos);
    EXPECT_NE(read_text_file(plain_vcf).find("\tBLACKLISTED_SITE\t"), std::string::npos);
    EXPECT_NE(read_text_file(plain_vcf).find("##INFO=<ID=KL,Number=1,Type=Float,"), std::string::npos);
    EXPECT_NE(read_text_file(plain_vcf).find("\t.\tBLACKLISTED_SITE\tKL=.\t"), std::string::npos);
    EXPECT_NE(read_text_file(plain_tsv).find("alt_sample\tchrM\t100\tA\tG\t"),
              std::string::npos);

    for (const auto& paths : {std::make_pair(output_vcf, plain_vcf),
                              std::make_pair(output_tsv, plain_tsv)}) {
        ngslib::BGZFile reader(paths.first, "r");
        EXPECT_EQ(bgzf_compression(reader.bgzf_handle()), compressed ? 2 : 0);
        if (compressed) EXPECT_EQ(bgzf_check_EOF(reader.bgzf_handle()), 1);
        std::string content, line;
        while (reader.readline(line)) content += line + "\n";
        EXPECT_EQ(content, read_text_file(paths.second));
    }

    if (compressed) {
        const int status = tbx_index_build(output_vcf.c_str(), 0, &tbx_conf_vcf);
        EXPECT_EQ(status, 0);
        if (status == 0) {
            ngslib::VCFFile reader(output_vcf, "r");
            reader.index_load();
            EXPECT_TRUE(reader.fetch("chrM:100-100"));
            ngslib::VCFRecord record;
            EXPECT_GE(reader.read(record), 0);
            EXPECT_EQ(record.pos(), 99);
            EXPECT_EQ(reader.read(record), -1);
        }
        std::remove((output_vcf + ".tbi").c_str());
    }

    for (const auto& path : {input, plain_vcf, plain_tsv, output_vcf, output_tsv}) {
        std::remove(path.c_str());
    }
}

TEST_P(VariantQCOutput, EmptyInputProducesValidHeaderOnlyOutputs) {
    const std::string suffix = GetParam();
    const bool compressed = suffix == ".gz" || suffix == ".bgz";
    const std::string prefix = "/tmp/mitoquest_variant_qc_empty" + suffix;
    const std::string input = prefix + "_input.vcf";
    const std::string output_vcf = prefix + "_out.vcf" + suffix;
    const std::string output_tsv = prefix + "_out.tsv" + suffix;
    write_text_file(input, "##fileformat=VCFv4.2\n"
                           "##contig=<ID=chrM,length=16569>\n"
                           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    VariantQC(make_variant_qc_config(input, output_vcf, output_tsv)).run();

    {
        ngslib::VCFFile reader(output_vcf, "r");
        ngslib::VCFRecord record;
        EXPECT_EQ(reader.read(record), -1);
        EXPECT_EQ(reader.header().n_samples(), 0);
    }
    for (const auto& path : {output_vcf, output_tsv}) {
        ngslib::BGZFile reader(path, "r");
        EXPECT_EQ(bgzf_compression(reader.bgzf_handle()), compressed ? 2 : 0);
        if (compressed) EXPECT_EQ(bgzf_check_EOF(reader.bgzf_handle()), 1);
    }
    {
        ngslib::BGZFile reader(output_tsv, "r");
        std::string line;
        EXPECT_TRUE(reader.readline(line));
        EXPECT_EQ(line, "sample\tchrom\tpos\tref\talt\tploidy\tHQ\tSRF\t"
                        "q_alpha\tq_beta\tp_error\tposterior\tis_mutation");
        EXPECT_FALSE(reader.readline(line));
    }
    for (const auto& path : {input, output_vcf, output_tsv}) {
        std::remove(path.c_str());
    }
}

INSTANTIATE_TEST_SUITE_P(Formats, VariantQCOutput,
                        ::testing::Values("", ".gz", ".bgz", ".txt"));

namespace {

const char* const h1_vcf_header = R"(##fileformat=VCFv4.2
##contig=<ID=chrM,length=16569>
##FILTER=<ID=REF_ONLY,Description="Reference-only input record">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="GT-aligned depths">
##FORMAT=<ID=AF,Number=.,Type=Float,Description="GT-aligned frequencies">
##FORMAT=<ID=AQ,Number=.,Type=Integer,Description="Allele quality">
##FORMAT=<ID=SB,Number=1,Type=String,Description="Strand counts">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	background	query
)";

std::vector<std::string> split_qc_fields(const std::string& text, char delimiter) {
    std::vector<std::string> fields;
    std::istringstream input(text);
    std::string field;
    while (std::getline(input, field, delimiter)) fields.push_back(field);
    return fields;
}

std::map<std::string, std::string> parse_qc_info(const std::string& text) {
    std::map<std::string, std::string> info;
    for (const auto& field : split_qc_fields(text, ';')) {
        if (field == ".") continue;
        const auto separator = field.find('=');
        const auto key = field.substr(0, separator);
        const auto value = separator == std::string::npos ? "" : field.substr(separator + 1);
        EXPECT_TRUE(info.emplace(key, value).second) << "Duplicate INFO/" << key;
    }
    return info;
}

struct H1OutputRecord {
    std::string qual;
    std::string filter;
    std::map<std::string, std::string> info;
    std::vector<std::map<std::string, std::string>> samples;
};

struct H1Iteration {
    int round;
    size_t selected;
    size_t changed;
    size_t evaluable;
    double fraction;
    double max_pp_change;
    double max_log_parameter_change;
};

struct H1Run {
    VariantQC::BetaParams beta;
    std::string log;
    std::string tsv;
    std::map<int, H1OutputRecord> records;
    std::vector<H1Iteration> iterations;
    VariantQC::Config config;
    std::string trace;
    std::vector<std::vector<std::string>> trace_rows;
};

// Restore stderr even if run() or output parsing throws.
class H1StderrCapture {
public:
    H1StderrCapture() { ::testing::internal::CaptureStderr(); }
    ~H1StderrCapture() {
        if (active_) (void)::testing::internal::GetCapturedStderr();
    }
    std::string finish() {
        active_ = false;
        return ::testing::internal::GetCapturedStderr();
    }
private:
    bool active_ = true;
};

class VariantQCH1Training : public ::testing::Test {
protected:
    void SetUp() override {
        char directory[] = "/tmp/mitoquest_variant_qc_h1_XXXXXX";
        const char* created = ::mkdtemp(directory);
        ASSERT_NE(created, nullptr);
        directory_ = created;
    }

    void TearDown() override {
        if (directory_.empty()) return;
        std::error_code error;
        std::filesystem::remove_all(directory_, error);
        EXPECT_FALSE(error) << error.message();
    }

    std::string path(const std::string& filename) const {
        return directory_ + "/" + filename;
    }

    H1Run run(const std::string& body, int max_iter,
              double threshold = 0.9, double epsilon = 0.001,
              int dp_threshold = 1024, const std::string& trace_suffix = "",
              const std::string& header = h1_vcf_header) {
        const std::string prefix = path("run_" + std::to_string(next_run_++));
        auto config = make_variant_qc_config(prefix + ".vcf", prefix + "_out.vcf",
                                             prefix + "_out.tsv");
        config.max_iter = max_iter;
        config.threshold = threshold;
        config.convergence_eps = epsilon;
        config.dp_threshold = dp_threshold;
        if (!trace_suffix.empty()) config.training_trace = prefix + "_training" + trace_suffix;
        write_text_file(config.input_vcf, header + body);
        H1StderrCapture capture;
        VariantQC(config).run();
        H1Run result{};
        result.config = config;
        result.log = capture.finish();
        result.tsv = read_text_file(config.output_tsv);
        if (!config.training_trace.empty()) {
            ngslib::BGZFile reader(config.training_trace, "r");
            std::string line;
            while (reader.readline(line)) {
                result.trace += line + "\n";
                if (!line.empty() && line.front() >= '0' && line.front() <= '9') {
                    result.trace_rows.push_back(split_qc_fields(line, '\t'));
                }
            }
        }

        const std::regex beta_pattern(R"(Final Beta: alpha=([^\s]+)\s+beta=([^\s]+))");
        std::smatch match;
        if (!std::regex_search(result.log, match, beta_pattern)) {
            throw std::runtime_error("Missing final H1 parameters:\n" + result.log);
        }
        result.beta = {std::stod(match[1].str()), std::stod(match[2].str())};
        EXPECT_TRUE(std::isfinite(result.beta.alpha));
        EXPECT_TRUE(std::isfinite(result.beta.beta));
        EXPECT_GT(result.beta.alpha, 0.0);
        EXPECT_GT(result.beta.beta, 0.0);

        const std::regex iteration_pattern(
            R"(H1 iteration ([0-9]+): selected=([0-9]+) changed=([0-9]+)/([0-9]+) \(([^)]+)\) max_pp_change=([^\s]+) max_log_parameter_change=([^\s]+))");
        for (std::sregex_iterator it(result.log.begin(), result.log.end(), iteration_pattern), end;
             it != end; ++it) {
            const auto& m = *it;
            result.iterations.push_back({std::stoi(m[1].str()),
                static_cast<size_t>(std::stoull(m[2].str())),
                static_cast<size_t>(std::stoull(m[3].str())),
                static_cast<size_t>(std::stoull(m[4].str())),
                std::stod(m[5].str()), std::stod(m[6].str()), std::stod(m[7].str())});
        }

        std::istringstream output(read_text_file(config.output_vcf));
        std::string line;
        while (std::getline(output, line)) {
            if (line.empty() || line.front() == '#') continue;
            const auto fields = split_qc_fields(line, '\t');
            if (fields.size() != 11) throw std::runtime_error("Unexpected VCF row: " + line);
            H1OutputRecord record{fields[5], fields[6], parse_qc_info(fields[7]), {}};
            const auto format = split_qc_fields(fields[8], ':');
            for (size_t sample = 9; sample < fields.size(); ++sample) {
                const auto values = split_qc_fields(fields[sample], ':');
                if (values.size() != format.size()) {
                    throw std::runtime_error("Incomplete output FORMAT: " + line);
                }
                std::map<std::string, std::string> data;
                for (size_t j = 0; j < format.size(); ++j) data[format[j]] = values[j];
                record.samples.push_back(std::move(data));
            }
            EXPECT_TRUE(result.records.emplace(std::stoi(fields[1]), std::move(record)).second);
        }
        return result;
    }

    // AF=1 in the sole qualifying background sample forces the Beta(2.5, 1680) fallback.
    static std::string row(int pos, const std::string& query,
                           const std::string& alt = "G",
                           const std::string& filter = "PASS",
                           const std::string& background = "0:2048:2048:1:30:1024,1024",
                           const std::string& info = ".") {
        return "chrM\t" + std::to_string(pos) + "\t.\tA\t" + alt + "\t60\t" + filter +
            "\t" + info + "\tGT:DP:AD:AF:AQ:SB\t" + background + "\t" + query + "\n";
    }

    static std::string sample(int depth, int count, bool with_ref) {
        // Even counts and power-of-two depths keep AF exact and strand ratios equal to 1.
        std::ostringstream output;
        output << std::setprecision(17) << (with_ref ? "0/1:" : "1:") << depth << ':';
        if (with_ref) output << depth - count << ',';
        output << count << ':';
        if (with_ref) output << static_cast<double>(depth - count) / depth << ',';
        output << static_cast<double>(count) / depth << ':';
        if (with_ref) output << "30,";
        output << "30:";
        if (with_ref) output << (depth - count) / 2 << ',' << (depth - count) / 2 << ';';
        output << count / 2 << ',' << count / 2;
        return output.str();
    }

    static std::vector<int> training_counts() { return {16, 32, 64, 128, 256, 384, 1000}; }

    static std::string training_rows(int probe_count = 6) {
        std::string body;
        int pos = 1100;
        for (int count : training_counts()) body += row(pos++, sample(1024, count, true));
        // Below the initial DP gate, but evaluable during refitting/final calling.
        return body + row(pos, sample(128, probe_count, false));
    }

    static std::string reference_rows(const std::string& filter, int count = 96) {
        std::string body;
        for (int i = 0; i < count; ++i) {
            // Mostly ALT=.; also exercise GT=0 at a site with a declared ALT.
            body += row(2000 + i, "0:512:512:1:30:256,256",
                        i < count * 3 / 4 ? "." : "G", filter);
        }
        return body;
    }

    static double posterior(int count, int depth, const VariantQC::BetaParams& beta) {
        return VariantQC::bayesian_filter(count, depth, 1.0, 30.0, 20,
                                          2.5, 1680.0, beta.alpha, beta.beta,
                                          5e-8 * 16569);
    }

    static double quality(int count, int depth, const VariantQC::BetaParams& beta) {
        return VariantQC::phred_from_log_odds(VariantQC::bayesian_log_odds(
            count, depth, 1.0, 30.0, 20, 2.5, 1680.0, beta.alpha, beta.beta, 5e-8 * 16569));
    }

    static double expected_kl(const std::vector<double>& vafs) {
        std::vector<double> edges(101);
        for (int i = 0; i <= 100; ++i) edges[i] = static_cast<double>(i) / 100;
        const double kl = VariantQC::kl_divergence(vafs, 2.5, 1680.0, edges);
        return std::isfinite(kl) ? kl : 1000.0;
    }

    static void expect_vcf_float(const std::string& text, double expected) {
        ASSERT_NE(text, ".");
        ASSERT_TRUE(std::isfinite(expected));
        const double actual = std::stod(text);
        EXPECT_TRUE(std::isfinite(actual));
        // Allow Float serialization error, but not integer rounding or score clipping.
        EXPECT_NEAR(actual, expected, std::max(1e-12, std::abs(expected) * 5e-6));
    }

    static VariantQC::BetaParams fitted(const std::vector<int>& counts) {
        std::vector<double> vafs;
        for (int count : counts) vafs.push_back(static_cast<double>(count) / 1024);
        mitoquest::LogFactorial lf(100000);
        return VariantQC::fit_betabinom_mle(vafs, counts,
                                           std::vector<int>(counts.size(), 1024), lf);
    }

    static void expect_beta(const VariantQC::BetaParams& actual,
                            const VariantQC::BetaParams& expected) {
        // The final-parameter log must retain enough digits to round-trip.
        EXPECT_DOUBLE_EQ(actual.alpha, expected.alpha);
        EXPECT_DOUBLE_EQ(actual.beta, expected.beta);
    }

    static void expect_sample(const H1OutputRecord& record, size_t index,
                              const std::string& gt, double pp, double threshold = 0.9) {
        const auto& output = record.samples.at(index);
        EXPECT_EQ(output.at("GT"), gt);
        // FORMAT/PP is Float and htslib prints a finite number of digits.
        EXPECT_NEAR(std::stod(output.at("PP")), pp, 1e-6);
        EXPECT_EQ(output.at("GOOD_CALL"), pp > threshold ? "True" : "False");
    }

    static void expect_final_training_calls(const H1Run& result, int probe_count = 6) {
        int pos = 1100;
        for (int count : training_counts()) {
            SCOPED_TRACE(pos);
            const auto& record = result.records.at(pos++);
            const double ref_pp = posterior(1024 - count, 1024, result.beta);
            const double alt_pp = posterior(count, 1024, result.beta);
            const std::string gt = std::string(ref_pp > 0.9 ? "0" : ".") + "/" +
                                   (alt_pp > 0.9 ? "1" : ".");
            expect_sample(record, 0, "0", posterior(2048, 2048, result.beta));
            expect_sample(record, 1, gt, std::max(ref_pp, alt_pp));
            expect_vcf_float(record.qual, quality(count, 1024, result.beta));
            EXPECT_EQ(record.filter, "PASS");
        }
        const double pp = posterior(probe_count, 128, result.beta);
        const auto& probe = result.records.at(pos);
        expect_sample(probe, 0, "0", posterior(2048, 2048, result.beta));
        expect_sample(probe, 1, pp > 0.9 ? "1" : ".", pp);
        expect_vcf_float(probe.qual, quality(probe_count, 128, result.beta));
        EXPECT_EQ(probe.filter, "PASS");
    }

    static std::string precise(double value) {
        std::ostringstream output;
        output << std::setprecision(17) << value;
        return output.str();
    }

    static std::vector<std::pair<int, int>> training_observations(int probe_count = 6) {
        std::vector<std::pair<int, int>> observations;
        for (int count : training_counts()) observations.emplace_back(count, 1024);
        observations.emplace_back(probe_count, 128);
        return observations;
    }

    static void expect_trace(const H1Run& result, size_t initial_training_count,
                             const std::vector<std::pair<int, int>>& observations,
                             const std::string& stop_reason) {
        const auto& config = result.config;
        ASSERT_FALSE(config.training_trace.empty());
        const auto lines = split_qc_fields(result.trace, '\n');
        const std::vector<std::string> prefix = {
            "#schema=mitoquest.variant_qc.training.v1",
            "#max_iter=" + std::to_string(config.max_iter),
            "#convergence_eps=" + precise(config.convergence_eps),
            "#threshold=" + precise(config.threshold),
            "#pi=" + precise(config.pi),
            "#dp_threshold=" + std::to_string(config.dp_threshold),
            "#hq_threshold=" + std::to_string(config.hq_threshold),
            "#max_alt_alleles=" + std::to_string(config.max_alt_alleles),
            "#bins=" + std::to_string(config.bins),
            "iteration\tstage\talpha\tbeta\ttraining_count\tevaluable_alt_count\t"
            "selected_alt_count\tchanged_alt_count\tcall_change_fraction\t"
            "max_pp_change\tmax_log_parameter_change"
        };
        ASSERT_GE(lines.size(), prefix.size() + 2);
        for (size_t i = 0; i < prefix.size(); ++i) EXPECT_EQ(lines[i], prefix[i]);
        EXPECT_EQ(lines.back(), "#stop_reason=" + stop_reason);
        EXPECT_NE(result.log.find("H1 stop: " + stop_reason + "\n"), std::string::npos);
        EXPECT_NE(result.log.find("Training: " + std::to_string(initial_training_count) +
                                  " VAFs."), std::string::npos);
        EXPECT_NE(result.log.find("H1 evaluable ALT observations: " +
                                  std::to_string(observations.size()) + "\n"),
                  std::string::npos);
        ASSERT_EQ(result.trace_rows.size(), lines.size() - prefix.size() - 1);
        ASSERT_EQ(result.trace_rows.size(), result.iterations.size() + 1);

        VariantQC::BetaParams previous{};
        std::vector<double> previous_posteriors;
        size_t previous_selected = 0;
        for (size_t i = 0; i < result.trace_rows.size(); ++i) {
            SCOPED_TRACE(i + 1);
            const auto& fields = result.trace_rows[i];
            ASSERT_EQ(fields.size(), 11u);
            EXPECT_EQ(split_qc_fields(lines[prefix.size() + i], '\t'), fields);
            EXPECT_EQ(fields[0], std::to_string(i + 1));
            EXPECT_EQ(fields[1], i == 0 ? "initial" : "refit");
            const VariantQC::BetaParams beta{std::stod(fields[2]), std::stod(fields[3])};
            EXPECT_TRUE(std::isfinite(beta.alpha));
            EXPECT_TRUE(std::isfinite(beta.beta));
            EXPECT_GT(beta.alpha, 0.0);
            EXPECT_GT(beta.beta, 0.0);
            EXPECT_EQ(fields[2], precise(beta.alpha));
            EXPECT_EQ(fields[3], precise(beta.beta));
            EXPECT_EQ(fields[4], std::to_string(i == 0 ? initial_training_count : previous_selected));
            EXPECT_EQ(fields[5], std::to_string(observations.size()));

            std::vector<double> posteriors;
            size_t selected = 0;
            for (const auto& [count, depth] : observations) {
                const double pp = posterior(count, depth, beta);
                posteriors.push_back(pp);
                if (pp > config.threshold) ++selected;
            }
            EXPECT_EQ(fields[6], std::to_string(selected));
            if (i == 0) {
                for (size_t field = 7; field < 11; ++field) EXPECT_EQ(fields[field], ".");
            } else {
                ASSERT_FALSE(observations.empty());
                size_t changed = 0;
                double max_pp_change = 0.0;
                for (size_t j = 0; j < posteriors.size(); ++j) {
                    if ((posteriors[j] > config.threshold) !=
                        (previous_posteriors[j] > config.threshold)) ++changed;
                    max_pp_change = std::max(max_pp_change,
                        std::abs(posteriors[j] - previous_posteriors[j]));
                }
                const double fraction = static_cast<double>(changed) / observations.size();
                const double max_log_change = std::max(
                    std::abs(std::log(beta.alpha) - std::log(previous.alpha)),
                    std::abs(std::log(beta.beta) - std::log(previous.beta)));
                EXPECT_EQ(fields[7], std::to_string(changed));
                EXPECT_EQ(fields[8], precise(fraction));
                EXPECT_EQ(fields[9], precise(max_pp_change));
                EXPECT_EQ(fields[10], precise(max_log_change));
                const auto& iteration = result.iterations[i - 1];
                EXPECT_EQ(iteration.round, static_cast<int>(i + 1));
                // The log's selected count describes this fit's input, not its new selections.
                EXPECT_EQ(iteration.selected, previous_selected);
                EXPECT_EQ(iteration.changed, changed);
                EXPECT_EQ(iteration.evaluable, observations.size());
                EXPECT_NEAR(iteration.fraction, fraction, 1e-6);
                EXPECT_NEAR(iteration.max_pp_change, max_pp_change, 1e-6);
                EXPECT_NEAR(iteration.max_log_parameter_change, max_log_change,
                            1e-5 * std::max(1.0, max_log_change));
            }
            previous = beta;
            previous_posteriors = std::move(posteriors);
            previous_selected = selected;
        }
        EXPECT_EQ(previous.alpha, result.beta.alpha);
        EXPECT_EQ(previous.beta, result.beta.beta);
    }

    static void expect_same_outputs(const H1Run& actual, const H1Run& expected) {
        EXPECT_EQ(actual.beta.alpha, expected.beta.alpha);
        EXPECT_EQ(actual.beta.beta, expected.beta.beta);
        EXPECT_EQ(actual.tsv, expected.tsv);
        ASSERT_EQ(actual.records.size(), expected.records.size());
        for (const auto& [pos, record] : expected.records) {
            SCOPED_TRACE(pos);
            const auto& other = actual.records.at(pos);
            EXPECT_EQ(other.qual, record.qual);
            EXPECT_EQ(other.filter, record.filter);
            EXPECT_EQ(other.info, record.info);
            EXPECT_EQ(other.samples, record.samples);
        }
    }

private:
    std::string directory_;
    int next_run_ = 0;
};

TEST_F(VariantQCH1Training, ReferenceOnlyRecordsCannotChangeTrainingOrSharedCalls) {
    const auto baseline = run(training_rows(), 4);
    ASSERT_FALSE(baseline.iterations.empty());
    for (const std::string filter : {"REF_ONLY", "PASS", "."}) {
        SCOPED_TRACE(filter);
        const auto expanded = run(training_rows() + reference_rows(filter), 4);
        expect_beta(expanded.beta, baseline.beta);
        ASSERT_EQ(expanded.iterations.size(), baseline.iterations.size());
        for (size_t i = 0; i < baseline.iterations.size(); ++i) {
            EXPECT_EQ(expanded.iterations[i].selected, baseline.iterations[i].selected);
            EXPECT_EQ(expanded.iterations[i].changed, baseline.iterations[i].changed);
            EXPECT_EQ(expanded.iterations[i].evaluable, baseline.iterations[i].evaluable);
            EXPECT_DOUBLE_EQ(expanded.iterations[i].fraction, baseline.iterations[i].fraction);
            EXPECT_DOUBLE_EQ(expanded.iterations[i].max_pp_change,
                             baseline.iterations[i].max_pp_change);
            EXPECT_DOUBLE_EQ(expanded.iterations[i].max_log_parameter_change,
                             baseline.iterations[i].max_log_parameter_change);
        }
        ASSERT_EQ(expanded.records.size(), baseline.records.size() + 96);
        for (const auto& [pos, record] : baseline.records) {
            SCOPED_TRACE(pos);
            const auto& other = expanded.records.at(pos);
            EXPECT_EQ(other.qual, record.qual);
            EXPECT_EQ(other.filter, record.filter);
            EXPECT_EQ(other.info, record.info);
            EXPECT_EQ(other.samples, record.samples); // includes GT, PP and GOOD_CALL
        }
        for (int pos = 2000; pos < 2096; ++pos) {
            const auto& record = expanded.records.at(pos);
            expect_sample(record, 0, "0", posterior(2048, 2048, expanded.beta));
            expect_sample(record, 1, "0", posterior(512, 512, expanded.beta));
            EXPECT_EQ(record.filter, "PASS");
        }
        EXPECT_FALSE(baseline.tsv.empty());
        EXPECT_EQ(expanded.tsv, baseline.tsv); // no added pure-REF TSV rows
    }
}

TEST_F(VariantQCH1Training, OneAndTwoFittingRoundsUseFrozenParametersForAllCalls) {
    const auto once = run(training_rows(), 1);
    const auto twice = run(training_rows(), 2);
    // Initial balancing keeps the first low AF and the sole AF>0.95.
    expect_beta(once.beta, fitted({16, 1000}));
    // The weak ALT at 16/1024 and the below-gate probe are initially rejected.
    expect_beta(twice.beta, fitted({32, 64, 128, 256, 384, 1000}));
    EXPECT_TRUE(once.iterations.empty());
    ASSERT_EQ(twice.iterations.size(), 1u);
    EXPECT_EQ(twice.iterations[0].selected, 6u);
    for (const auto* result : {&once, &twice}) {
        EXPECT_NE(result->log.find("H1 stop: maximum fitting rounds reached"), std::string::npos);
        expect_final_training_calls(*result);
    }
    // Intermediate rejection must not prevent a final call when PP crosses the threshold.
    ASSERT_LT(posterior(6, 128, once.beta), 0.9);
    ASSERT_GT(posterior(6, 128, twice.beta), 0.9);
    EXPECT_EQ(once.records.at(1107).samples[1].at("GT"), ".");
    EXPECT_EQ(twice.records.at(1107).samples[1].at("GT"), "1");
    EXPECT_LT(std::stod(once.records.at(1107).qual), 10.0);
    EXPECT_GT(std::stod(twice.records.at(1107).qual), 10.0);
    EXPECT_EQ(once.tsv.find("query\tchrM\t1107\t"), std::string::npos);
    EXPECT_NE(twice.tsv.find("query\tchrM\t1107\t"), std::string::npos);
}

TEST_F(VariantQCH1Training, ReferenceGtUsesFinalModelAfterEarlyRejection) {
    for (const std::string alt : {".", "G"}) {
        SCOPED_TRACE(alt);
        const std::string filter = alt == "." ? "REF_ONLY" : "PASS";
        const std::string body = training_rows(4) +
            row(2200, "0:128:6:0.046875:30:3,3", alt, filter) +
            row(2201, "0:128:2:0.015625:30:1,1", alt, filter);
        const auto initial = run(body, 1);
        const auto limited = run(body, 2);
        const auto converged = run(body, 4);
        expect_beta(initial.beta, fitted({16, 1000}));
        expect_beta(limited.beta, fitted({32, 64, 128, 256, 384, 1000}));
        expect_beta(converged.beta, limited.beta);
        const double initial_pp = posterior(6, 128, initial.beta);
        ASSERT_LT(initial_pp, 0.9);
        expect_sample(initial.records.at(2200), 1, ".", initial_pp);
        ASSERT_EQ(limited.iterations.size(), 1u);
        ASSERT_EQ(converged.iterations.size(), 2u);
        EXPECT_NE(limited.log.find("H1 stop: maximum fitting rounds reached"), std::string::npos);
        EXPECT_NE(converged.log.find("H1 stop: converged"), std::string::npos);
        for (const auto* result : {&limited, &converged}) {
            const double final_pp = posterior(6, 128, result->beta);
            ASSERT_GT(final_pp, 0.9);
            expect_sample(result->records.at(2200), 1, "0", final_pp);
            const double rejected_pp = posterior(2, 128, result->beta);
            ASSERT_LT(rejected_pp, 0.9);
            expect_sample(result->records.at(2201), 1, ".", rejected_pp);
            expect_final_training_calls(*result, 4);
        }
    }
}

TEST_F(VariantQCH1Training, OriginallyMissingGtIsNotInventedWhenFinalPosteriorPasses) {
    const std::string body = training_rows(4) +
        row(2200, ".:128:6:0.046875:30:3,3");
    const auto initial = run(body, 1);
    ASSERT_LT(posterior(6, 128, initial.beta), 0.9);
    expect_sample(initial.records.at(2200), 1, ".", posterior(6, 128, initial.beta));
    EXPECT_EQ(initial.records.at(2200).qual, ".");
    for (int max_iter : {2, 4}) {
        SCOPED_TRACE(max_iter);
        const auto result = run(body, max_iter);
        expect_beta(result.beta, fitted({32, 64, 128, 256, 384, 1000}));
        const double final_pp = posterior(6, 128, result.beta);
        ASSERT_GT(final_pp, 0.9);
        expect_sample(result.records.at(2200), 1, ".", final_pp);
        EXPECT_EQ(result.records.at(2200).qual, ".");
    }
}

TEST_F(VariantQCH1Training, ReferenceOnlyInputStopsWithoutEvaluableAlt) {
    const auto result = run(reference_rows("REF_ONLY"), 4);
    expect_beta(result.beta, {0.5, 0.2});
    EXPECT_TRUE(result.iterations.empty());
    EXPECT_NE(result.log.find("H1 stop: no evaluable ALT observations"), std::string::npos);
    ASSERT_EQ(result.records.size(), 96u);
    for (const auto& [pos, record] : result.records) {
        SCOPED_TRACE(pos);
        expect_sample(record, 0, "0", posterior(2048, 2048, result.beta));
        expect_sample(record, 1, "0", posterior(512, 512, result.beta));
        EXPECT_EQ(record.filter, "PASS");
        EXPECT_EQ(record.qual, ".");
        expect_vcf_float(record.info.at("KL"), 0.0);
    }
    EXPECT_EQ(std::count(result.tsv.begin(), result.tsv.end(), '\n'), 1);
}

TEST_F(VariantQCH1Training, OneOrTwoSelectedAltsKeepInitialParameters) {
    for (int count : {1, 2}) {
        SCOPED_TRACE(count);
        std::string body;
        for (int i = 0; i < count; ++i) body += row(1100 + i, sample(1024, 128 << i, false));
        // All query DPs are below the initial gate; initial H1 is the fallback.
        const auto result = run(body, 4, 0.9, 0.001, 2048);
        expect_beta(result.beta, {0.5, 0.2});
        EXPECT_TRUE(result.iterations.empty());
        EXPECT_NE(result.log.find("H1 stop: fewer than 3 selected ALT observations"),
                  std::string::npos);
        for (int i = 0; i < count; ++i) {
            const double pp = posterior(128 << i, 1024, result.beta);
            ASSERT_GT(pp, 0.9);
            expect_sample(result.records.at(1100 + i), 1, "1", pp);
        }
    }
}

TEST_F(VariantQCH1Training, ReferencePosteriorDoesNotSelectLowPosteriorAlt) {
    std::string body;
    for (int i = 0; i < 3; ++i) body += row(1100 + i, sample(1024, 2, true));
    const auto result = run(body, 4, 0.9, 0.001, 2048);
    expect_beta(result.beta, {0.5, 0.2});
    ASSERT_GT(posterior(1022, 1024, result.beta), 0.9);
    ASSERT_LT(posterior(2, 1024, result.beta), 0.9);
    EXPECT_TRUE(result.iterations.empty());
    EXPECT_NE(result.log.find("H1 stop: fewer than 3 selected ALT observations"),
              std::string::npos);
    for (const auto& [pos, record] : result.records) {
        SCOPED_TRACE(pos);
        expect_sample(record, 1, "0/.", posterior(1022, 1024, result.beta));
        EXPECT_DOUBLE_EQ(std::stod(record.samples[1].at("PP")), 1.0);
        EXPECT_EQ(record.filter, "PASS");
        expect_vcf_float(record.qual, quality(2, 1024, result.beta));
        EXPECT_LT(std::stod(record.qual), 1.0); // REF PP=1 cannot inflate ALT QUAL.
    }
}

TEST_F(VariantQCH1Training, StableAltDecisionsRequirePosteriorAndParameterConvergence) {
    for (double epsilon : {0.001, 0.1}) {
        SCOPED_TRACE(epsilon);
        const auto result = run(training_rows(/*probe_count=*/4), 4, 0.9, epsilon);
        ASSERT_EQ(result.iterations.size(), 2u) << result.log;
        const auto& first = result.iterations[0];
        EXPECT_EQ(first.selected, 6u);
        EXPECT_EQ(first.changed, 0u);
        EXPECT_EQ(first.evaluable, 8u);
        EXPECT_DOUBLE_EQ(first.fraction, 0.0);
        if (epsilon == 0.001) {
            EXPECT_GT(first.max_pp_change, epsilon);
        } else {
            // Stable calls AND a small PP change are still not sufficient.
            EXPECT_LT(first.max_pp_change, epsilon);
        }
        EXPECT_GT(first.max_log_parameter_change, epsilon);
        const auto& last = result.iterations[1];
        EXPECT_GT(last.round, first.round);
        EXPECT_EQ(last.selected, 6u);
        EXPECT_EQ(last.changed, 0u);
        EXPECT_LT(last.max_pp_change, epsilon);
        EXPECT_LT(last.max_log_parameter_change, epsilon);
        EXPECT_NE(result.log.find("H1 stop: converged"), std::string::npos);
        EXPECT_EQ(result.log.find("H1 stop: maximum fitting rounds reached"), std::string::npos);
        expect_beta(result.beta, fitted({32, 64, 128, 256, 384, 1000}));
        expect_final_training_calls(result, 4);
    }
}

TEST_F(VariantQCH1Training, OnlyEvaluableAltsContributeToIterationDiagnostics) {
    std::string extra = reference_rows("PASS");
    extra += row(2200, "0/1:128:96,32:.:30,30:48,48;16,16"); // prefiltered AF
    extra += row(2201, "0/.:128:96,32:.:30,30:48,48;16,16"); // normalized missing GT
    extra += row(2202, "0/1:.:96,32:0.75,0.25:30,30:48,48;16,16"); // missing DP
    extra += row(2203, "0/1:128:96,32:0.75,0.25:30:48,48;16,16"); // short AQ
    extra += row(2204, "0/1:128:96,32:0.75,0.25:30,30:48,48"); // short SRF
    extra += row(2205, "0/1:128:96:0.75,0.25:30,30:48,48;16,16"); // short AD
    extra += row(2206, "1:0:1:0.5:30:1,1"); // DP=0
    extra += row(2207, "1:128:0:0:30:0,0"); // AD=0
    extra += row(2208, "1:128:128:1:30:64,64"); // AD=DP
    extra += row(2209, "1:128:130:1:30:65,65"); // AD>DP
    extra += row(2210, ".:128:32:0.25:30:16,16"); // no original ALT GT
    extra += row(2211, "0/1:128:.:0.75,0.25:30,30:48,48;16,16"); // prefiltered AD
    const auto result = run(training_rows() + extra, 2);
    const auto initial = fitted({16, 1000});
    expect_beta(result.beta, fitted({32, 64, 128, 256, 384, 1000}));
    ASSERT_EQ(result.iterations.size(), 1u) << result.log;
    const auto& iteration = result.iterations[0];
    EXPECT_EQ(iteration.selected, 6u);
    EXPECT_EQ(iteration.changed, 1u);
    EXPECT_EQ(iteration.evaluable, 8u); // includes the below-initial-DP-gate probe
    EXPECT_DOUBLE_EQ(iteration.fraction, 1.0 / 8.0);
    double max_pp_change = std::abs(posterior(6, 128, initial) - posterior(6, 128, result.beta));
    for (int count : training_counts()) {
        max_pp_change = std::max(max_pp_change,
            std::abs(posterior(count, 1024, initial) - posterior(count, 1024, result.beta)));
    }
    EXPECT_NEAR(iteration.max_pp_change, max_pp_change, 1e-6);
    const double parameter_change = std::max(
        std::abs(std::log(result.beta.alpha) - std::log(initial.alpha)),
        std::abs(std::log(result.beta.beta) - std::log(initial.beta)));
    EXPECT_NEAR(iteration.max_log_parameter_change, parameter_change, 1e-5);
    expect_final_training_calls(result);
    for (int pos : {2200, 2202, 2211}) expect_sample(result.records.at(pos), 1, "0/1", 0.0);
    expect_sample(result.records.at(2201), 1, ".", 0.0);
    // Exclusion from H1 does not exclude evaluable observations from final calls.
    expect_sample(result.records.at(2208), 1, "1", posterior(128, 128, result.beta));
    expect_sample(result.records.at(2210), 1, ".", posterior(32, 128, result.beta));
    for (int pos : {2200, 2201, 2202, 2203, 2204, 2205, 2206, 2209, 2210, 2211}) {
        SCOPED_TRACE(pos);
        EXPECT_EQ(result.records.at(pos).qual, ".");
    }
    // AD=0 and AD=DP are scoreable for QUAL, unlike interior-only H1 refitting.
    expect_vcf_float(result.records.at(2207).qual, quality(0, 128, result.beta));
    expect_vcf_float(result.records.at(2208).qual, quality(128, 128, result.beta));
    EXPECT_GT(std::stod(result.records.at(2207).qual), 0.0);
    expect_vcf_float(result.records.at(2207).info.at("KL"), 0.0);
    expect_vcf_float(result.records.at(2208).info.at("KL"), 0.0);
}

TEST_F(VariantQCH1Training, PosteriorEqualToThresholdIsNeitherSelectedNorCalled) {
    std::string body;
    for (int i = 0; i < 3; ++i) body += row(1100 + i, sample(1024, 128 << i, false));
    const auto result = run(body, 4, /*threshold=*/1.0, 0.001, 2048);
    expect_beta(result.beta, {0.5, 0.2});
    EXPECT_TRUE(result.iterations.empty());
    EXPECT_NE(result.log.find("H1 stop: fewer than 3 selected ALT observations"),
              std::string::npos);
    for (int i = 0; i < 3; ++i) {
        ASSERT_DOUBLE_EQ(posterior(128 << i, 1024, result.beta), 1.0);
        const auto& record = result.records.at(1100 + i);
        expect_sample(record, 0, ".", 1.0, 1.0);
        expect_sample(record, 1, ".", 1.0, 1.0);
        EXPECT_EQ(record.filter, "LOW_QUALITY");
        expect_vcf_float(record.qual, quality(128 << i, 1024, result.beta));
        if (i > 0) EXPECT_GT(std::stod(record.qual), std::stod(result.records.at(1099 + i).qual));
    }
    EXPECT_GT(std::stod(result.records.at(1102).qual), 1000.0);
    EXPECT_EQ(std::count(result.tsv.begin(), result.tsv.end(), '\n'), 1);
}

TEST_F(VariantQCH1Training, EqualityAtInteriorThresholdDoesNotSupplyThirdSelectedAlt) {
    const VariantQC::BetaParams initial{0.5, 0.2};
    const double threshold = posterior(2, 1024, initial);
    ASSERT_GT(threshold, 0.0);
    ASSERT_LT(threshold, 1.0);
    const std::string body = row(1100, sample(1024, 2, false)) +
                             row(1101, sample(1024, 128, false)) +
                             row(1102, sample(1024, 256, false));
    const auto result = run(body, 4, threshold, 0.001, 2048);
    expect_beta(result.beta, initial);
    EXPECT_TRUE(result.iterations.empty());
    EXPECT_NE(result.log.find("H1 stop: fewer than 3 selected ALT observations"),
              std::string::npos);
    // Use the unrounded public posterior as the exact threshold, not VCF/PP.
    expect_sample(result.records.at(1100), 1, ".", threshold, threshold);
    for (int i = 0; i < 2; ++i) {
        const double pp = posterior(128 << i, 1024, initial);
        ASSERT_GT(pp, threshold);
        expect_sample(result.records.at(1101 + i), 1, "1", pp, threshold);
    }
    EXPECT_EQ(std::count(result.tsv.begin(), result.tsv.end(), '\n'), 3);
}

TEST_F(VariantQCH1Training, QualIsMaximumAcrossOriginalGtAlignedAltsAndSamplesNotSum) {
    const std::string two_alts = "1/2:128:32,96:0.25,0.75:30,30:16,16;48,48";
    const std::string body = row(2200, two_alts, "G,T") +
        row(2201, "2/1:128:96,32:0.75,0.25:30,30:48,48;16,16", "G,T") +
        row(2202, two_alts, "G,T", "PASS", "1:2048:2048:1:30:1024,1024") +
        // A surplus, stronger AD slot without an original GT must not contribute.
        row(2203, "1:128:32,128:0.25,1:30,30:16,16;64,64");
    const auto result = run(body, 1, 0.9, 0.001, 2048);
    const double low = quality(32, 128, result.beta);
    const double high = quality(96, 128, result.beta);
    ASSERT_GT(high, low);
    ASSERT_GT(low, 1.0); // Summing either within or between samples is detectably different.
    for (int pos : {2200, 2201}) expect_vcf_float(result.records.at(pos).qual, high);
    const double other_sample = quality(2048, 2048, result.beta);
    ASSERT_GT(other_sample, high);
    expect_vcf_float(result.records.at(2202).qual, other_sample);
    expect_vcf_float(result.records.at(2203).qual, low);
}

TEST_F(VariantQCH1Training, SkippedRecordsHaveMissingQualAndKlWithoutShiftingFollowingRecords) {
    std::string overlapping = row(298, sample(128, 32, false));
    overlapping.replace(overlapping.find("\tA\tG\t"), 5, "\tAC\tG\t");
    const auto result = run(overlapping + row(316, sample(128, 32, false)) +
        row(2200, sample(128, 32, false), "G,T,C") +
        row(2201, sample(128, 32, false), "G", "PASS", "0:512:512:1:30:256,256") +
        row(2202, sample(128, 32, false)), 1);
    for (int pos : {298, 316, 2200, 2201}) {
        SCOPED_TRACE(pos);
        const auto& record = result.records.at(pos);
        EXPECT_EQ(record.qual, ".");
        EXPECT_EQ(record.info.at("KL"), ".");
        EXPECT_EQ(record.filter, pos == 316 ? "BLACKLISTED_SITE" : "LOW_QUALITY");
        for (const auto& data : record.samples) {
            EXPECT_EQ(data.count("PP"), 0u);
            EXPECT_EQ(data.count("GOOD_CALL"), 0u);
        }
    }
    expect_vcf_float(result.records.at(2202).qual, quality(32, 128, result.beta));
    expect_vcf_float(result.records.at(2202).info.at("KL"), expected_kl({0.25}));
    EXPECT_NE(result.log.find("No background samples at chrM:2201"), std::string::npos);
}

TEST_F(VariantQCH1Training, InfoKlIsUnroundedAndReplacesStaleValuesAndHeaderDefinitions) {
    const std::string background = "0:2048:2048:1:30:1024,1024";
    const double kl = expected_kl({126.0 / 128, 2.0 / 128});
    ASSERT_GT(std::abs(kl - std::round(kl)), 1e-3);
    for (bool incompatible : {false, true}) {
        SCOPED_TRACE(incompatible);
        std::string header = h1_vcf_header;
        header.insert(header.find("#CHROM"),
            std::string("##INFO=<ID=KL,Number=") + (incompatible ? ".,Type=String" : "1,Type=Float") +
            ",Description=\"Stale input KL definition\">\n"
            "##INFO=<ID=KEEP,Number=1,Type=Integer,Description=\"Unrelated annotation\">\n");
        const std::string stale = incompatible ? "KL=999,998;KEEP=7" : "KL=999;KEEP=7";
        const auto result = run(
            row(316, sample(128, 2, true), "G", "PASS", background, stale) +
            row(2200, sample(128, 2, true), "G", "PASS", background, stale) +
            row(2201, "0:128:128:1:30:64,64", ".", "REF_ONLY", background, stale),
            1, 0.9, 0.001, 2048, "", header);
        expect_vcf_float(result.records.at(2200).info.at("KL"), kl);
        expect_vcf_float(result.records.at(2200).qual, quality(2, 128, result.beta));
        expect_vcf_float(result.records.at(2201).info.at("KL"), 0.0);
        EXPECT_EQ(result.records.at(2201).qual, ".");
        EXPECT_EQ(result.records.at(316).info.at("KL"), ".");
        EXPECT_EQ(result.records.at(316).qual, ".");
        for (const auto& [pos, record] : result.records) {
            EXPECT_EQ(record.info.at("KEEP"), "7") << pos;
            EXPECT_EQ(record.info.size(), 2u);
        }

        int kl_headers = 0;
        std::string qual_header, low_quality_header;
        std::istringstream output(read_text_file(result.config.output_vcf));
        std::string line;
        while (std::getline(output, line)) {
            if (line.rfind("##INFO=<ID=KL,", 0) == 0) {
                ++kl_headers;
                EXPECT_NE(line.find("Number=1,"), std::string::npos);
                EXPECT_NE(line.find("Type=Float,"), std::string::npos);
                EXPECT_EQ(line.find("Stale input KL definition"), std::string::npos);
            }
            if (line.rfind("##variant_qc_qual=", 0) == 0) qual_header = line;
            if (line.rfind("##FILTER=<ID=LOW_QUALITY,", 0) == 0) low_quality_header = line;
        }
        EXPECT_EQ(kl_headers, 1);
        ASSERT_FALSE(qual_header.empty());
        for (const std::string term : {"approx", "strongest.*ALT|ALT.*strongest", "posterior", "Phred",
                                       "not.*(joint|calibrated)"}) {
            EXPECT_TRUE(std::regex_search(qual_header, std::regex(term, std::regex::icase)))
                << qual_header;
        }
        ASSERT_FALSE(low_quality_header.empty());
        EXPECT_FALSE(std::regex_search(low_quality_header, std::regex("KL|divergence", std::regex::icase)));
    }
}

TEST_F(VariantQCH1Training, DuplicateCoordinatesKeepQualAndKlByInputRecordOrdinal) {
    auto config = make_variant_qc_config(path("duplicates.vcf"), path("duplicates_out.vcf"),
                                         path("duplicates_out.tsv"));
    config.dp_threshold = 2048; // No query enters initial training; both models use their fallbacks.
    const std::string body = row(2200, sample(128, 2, false), "G") +
        row(2200, sample(128, 32, false), "G,T,C") + // skipped at a processed coordinate
        row(2200, sample(128, 64, false), "T") +
        row(2200, "0:128:128:1:30:64,64", "C");
    write_text_file(config.input_vcf, std::string(h1_vcf_header) + body);
    {
        H1StderrCapture capture;
        VariantQC(config).run();
    }
    // The fixture's POS-keyed parser intentionally cannot represent duplicates.
    std::vector<std::vector<std::string>> records;
    std::istringstream output(read_text_file(config.output_vcf));
    std::string line;
    while (std::getline(output, line)) {
        if (!line.empty() && line.front() != '#') records.push_back(split_qc_fields(line, '\t'));
    }
    ASSERT_EQ(records.size(), 4u);
    const std::vector<std::string> alts = {"G", "G,T,C", "T", "C"};
    for (size_t i = 0; i < records.size(); ++i) {
        ASSERT_EQ(records[i].size(), 11u);
        EXPECT_EQ(records[i][1], "2200");
        EXPECT_EQ(records[i][4], alts[i]);
    }
    expect_vcf_float(records[0][5], quality(2, 128, {0.5, 0.2}));
    expect_vcf_float(records[2][5], quality(64, 128, {0.5, 0.2}));
    EXPECT_LT(std::stod(records[0][5]), std::stod(records[2][5]));
    EXPECT_EQ(records[1][5], ".");
    EXPECT_EQ(records[3][5], ".");
    expect_vcf_float(parse_qc_info(records[0][7]).at("KL"), expected_kl({2.0 / 128}));
    expect_vcf_float(parse_qc_info(records[2][7]).at("KL"), expected_kl({0.5}));
    EXPECT_EQ(parse_qc_info(records[1][7]).at("KL"), ".");
    expect_vcf_float(parse_qc_info(records[3][7]).at("KL"), 0.0);
    // Do not assert duplicate-coordinate GT/FILTER isolation: that legacy lookup is unchanged.
}

TEST_F(VariantQCH1Training, TrainingTraceInitialAndRefitCountsDescribeTheirOwnModels) {
    const auto initial = fitted({16, 1000});
    for (int max_iter : {1, 2}) {
        SCOPED_TRACE(max_iter);
        const auto baseline = run(training_rows(), max_iter);
        const auto result = run(training_rows(), max_iter, 0.9, 0.001, 1024, ".tsv");
        ASSERT_NO_FATAL_FAILURE(expect_trace(result, 2, training_observations(),
                                             "maximum fitting rounds reached"));
        expect_same_outputs(result, baseline);
        ASSERT_EQ(result.trace_rows.size(), static_cast<size_t>(max_iter));
        const auto& first = result.trace_rows.front();
        EXPECT_EQ(std::stod(first[2]), initial.alpha);
        EXPECT_EQ(std::stod(first[3]), initial.beta);
        EXPECT_EQ(first[4], "2");
        EXPECT_EQ(first[5], "8");
        EXPECT_EQ(first[6], "6");
        if (max_iter == 2) {
            const auto& refit = result.trace_rows.back();
            EXPECT_EQ(refit[4], "6");
            EXPECT_EQ(refit[5], "8");
            EXPECT_EQ(refit[6], "7");
            EXPECT_EQ(refit[7], "1");
            EXPECT_EQ(refit[8], "0.125");
        }
    }
}

TEST_F(VariantQCH1Training, TrainingTraceConvergenceIncludesEveryCompletedFit) {
    for (double epsilon : {0.001, 0.1}) {
        SCOPED_TRACE(epsilon);
        const auto result = run(training_rows(4), 4, 0.9, epsilon, 1024, ".tsv");
        ASSERT_NO_FATAL_FAILURE(expect_trace(result, 2, training_observations(4), "converged"));
        ASSERT_EQ(result.trace_rows.size(), 3u);
        EXPECT_EQ(result.trace_rows[1][7], "0");
        EXPECT_GT(std::stod(result.trace_rows[1][10]), epsilon);
        const auto& last = result.trace_rows.back();
        EXPECT_EQ(last[4], "6");
        EXPECT_EQ(last[6], "6");
        EXPECT_EQ(last[7], "0");
        EXPECT_LT(std::stod(last[8]), epsilon);
        EXPECT_LT(std::stod(last[9]), epsilon);
        EXPECT_LT(std::stod(last[10]), epsilon);
    }
}

TEST_F(VariantQCH1Training, TrainingTraceReferenceOnlyInputStillRecordsInitialFallback) {
    const auto result = run(reference_rows("REF_ONLY"), 4, 0.9, 0.001, 1024, ".tsv");
    ASSERT_NO_FATAL_FAILURE(expect_trace(result, 0, {}, "no evaluable ALT observations"));
    ASSERT_EQ(result.trace_rows.size(), 1u);
    EXPECT_EQ(result.beta.alpha, 0.5);
    EXPECT_EQ(result.beta.beta, 0.2);
    EXPECT_EQ(result.trace_rows[0][4], "0");
    EXPECT_EQ(result.trace_rows[0][5], "0");
    EXPECT_EQ(result.trace_rows[0][6], "0");
}

TEST_F(VariantQCH1Training, TrainingTraceTooFewSelectedAltsDoesNotInventARefit) {
    for (int count : {1, 2}) {
        SCOPED_TRACE(count);
        std::string body;
        std::vector<std::pair<int, int>> observations;
        for (int i = 0; i < count; ++i) {
            body += row(1100 + i, sample(1024, 128 << i, false));
            observations.emplace_back(128 << i, 1024);
        }
        const auto result = run(body, 4, 0.9, 0.001, 2048, ".tsv");
        ASSERT_NO_FATAL_FAILURE(expect_trace(result, 0, observations,
                                             "fewer than 3 selected ALT observations"));
        ASSERT_EQ(result.trace_rows.size(), 1u);
        EXPECT_EQ(result.beta.alpha, 0.5);
        EXPECT_EQ(result.beta.beta, 0.2);
        EXPECT_EQ(result.trace_rows[0][6], std::to_string(count));
    }
}

TEST_F(VariantQCH1Training, TrainingTraceCountsOnlyOwnAltPosteriorStrictlyAboveThreshold) {
    const VariantQC::BetaParams initial{0.5, 0.2};
    const double interior_threshold = posterior(2, 1024, initial);
    ASSERT_GT(interior_threshold, 0.0);
    ASSERT_LT(interior_threshold, 0.9);
    const std::string body = row(1100, sample(1024, 2, true)) +
                             row(1101, sample(1024, 128, true)) +
                             row(1102, sample(1024, 256, true));
    for (double threshold : {0.9, interior_threshold, 1.0}) {
        SCOPED_TRACE(threshold);
        const auto result = run(body, 4, threshold, 0.012345678901234567, 2048, ".tsv");
        ASSERT_NO_FATAL_FAILURE(expect_trace(result, 0, {{2, 1024}, {128, 1024}, {256, 1024}},
                                             "fewer than 3 selected ALT observations"));
        ASSERT_EQ(result.trace_rows.size(), 1u);
        EXPECT_EQ(result.trace_rows[0][6], threshold == 1.0 ? "0" : "2");
        if (threshold < 1.0) {
            ASSERT_GT(posterior(1022, 1024, initial), threshold);
            EXPECT_EQ(result.records.at(1100).samples[1].at("GT"), "0/.");
        }
    }
}

TEST_F(VariantQCH1Training, TrainingTraceCompressionPreservesContentAndAllOutputs) {
    const auto baseline = run(training_rows(4), 4);
    EXPECT_TRUE(baseline.config.training_trace.empty());
    EXPECT_TRUE(baseline.trace.empty());
    EXPECT_TRUE(baseline.trace_rows.empty());
    size_t files = 0;
    for (const auto& entry : std::filesystem::directory_iterator(path(""))) {
        EXPECT_TRUE(entry.is_regular_file());
        ++files;
    }
    EXPECT_EQ(files, 3u);

    std::string plain_trace;
    for (const std::string suffix : {".tsv", ".tsv.gz", ".tsv.bgz"}) {
        SCOPED_TRACE(suffix);
        const auto result = run(training_rows(4), 4, 0.9, 0.001, 1024, suffix);
        ASSERT_NO_FATAL_FAILURE(expect_trace(result, 2, training_observations(4), "converged"));
        expect_same_outputs(result, baseline);
        const bool compressed = suffix != ".tsv";
        ngslib::BGZFile reader(result.config.training_trace, "r");
        EXPECT_EQ(bgzf_compression(reader.bgzf_handle()), compressed ? 2 : 0);
        if (compressed) {
            EXPECT_EQ(bgzf_check_EOF(reader.bgzf_handle()), 1);
            EXPECT_EQ(result.trace, plain_trace);
        } else {
            plain_trace = result.trace;
            EXPECT_EQ(read_text_file(result.config.training_trace), plain_trace);
        }
    }
}

TEST_F(VariantQCH1Training, TrainingTraceReferenceOnlyExpansionPreservesAllModelRows) {
    const auto baseline = run(training_rows(4), 4, 0.9, 0.001, 1024, ".tsv");
    ASSERT_NO_FATAL_FAILURE(expect_trace(baseline, 2, training_observations(4), "converged"));
    for (const std::string filter : {"REF_ONLY", "PASS", "."}) {
        SCOPED_TRACE(filter);
        const auto expanded = run(training_rows(4) + reference_rows(filter),
                                  4, 0.9, 0.001, 1024, ".tsv");
        EXPECT_EQ(expanded.trace, baseline.trace);
        EXPECT_EQ(expanded.trace_rows, baseline.trace_rows);
        EXPECT_EQ(expanded.beta.alpha, baseline.beta.alpha);
        EXPECT_EQ(expanded.beta.beta, baseline.beta.beta);
        EXPECT_EQ(expanded.tsv, baseline.tsv);
        ASSERT_EQ(expanded.records.size(), baseline.records.size() + 96);
        for (const auto& [pos, record] : baseline.records) {
            SCOPED_TRACE(pos);
            const auto& other = expanded.records.at(pos);
            EXPECT_EQ(other.qual, record.qual);
            EXPECT_EQ(other.filter, record.filter);
            EXPECT_EQ(other.info, record.info);
            EXPECT_EQ(other.samples, record.samples);
        }
    }
}

TEST_F(VariantQCH1Training, TrainingTraceRejectsDirectNormalizedAndRelativePathCollisions) {
    const auto config = make_variant_qc_config(path("input.vcf"), path("output.vcf"),
                                               path("output.tsv"));
    const std::string input = std::string(h1_vcf_header) + training_rows();
    write_text_file(config.input_vcf, input);
    ASSERT_TRUE(std::filesystem::create_directory(path("alias")));
    for (const auto& target : {config.input_vcf, config.output_vcf, config.output_tsv}) {
        SCOPED_TRACE(target);
        const auto filename = std::filesystem::path(target).filename().string();
        const std::vector<std::string> aliases = {
            target, path("./" + filename), path("alias/../" + filename),
            std::filesystem::relative(target, std::filesystem::current_path()).string()
        };
        for (const auto& alias : aliases) {
            SCOPED_TRACE(alias);
            auto colliding = config;
            colliding.training_trace = alias;
            EXPECT_THROW((void)VariantQC{colliding}, std::invalid_argument);
        }
    }
    EXPECT_EQ(read_text_file(config.input_vcf), input);
    EXPECT_FALSE(std::filesystem::exists(config.output_vcf));
    EXPECT_FALSE(std::filesystem::exists(config.output_tsv));
}

TEST_F(VariantQCH1Training, TrainingTraceRejectsSymlinkAndHardLinkPathCollisions) {
    const auto config = make_variant_qc_config(path("input.vcf"), path("output.vcf"),
                                               path("output.tsv"));
    const std::map<std::string, std::string> originals = {
        {config.input_vcf, std::string(h1_vcf_header) + training_rows()},
        {config.output_vcf, "existing VCF output\n"},
        {config.output_tsv, "existing TSV output\n"}
    };
    for (const auto& [target, content] : originals) write_text_file(target, content);
    int link_number = 0;
    for (bool hard_link : {false, true}) {
        for (const auto& [target, content] : originals) {
            SCOPED_TRACE(target);
            SCOPED_TRACE(hard_link ? "hard link" : "symlink");
            auto colliding = config;
            colliding.training_trace = path("link_" + std::to_string(link_number++) + ".tsv");
            std::error_code error;
            if (hard_link) std::filesystem::create_hard_link(target, colliding.training_trace, error);
            else std::filesystem::create_symlink(target, colliding.training_trace, error);
            if (error) GTEST_SKIP() << "Filesystem does not support this link: " << error.message();
            EXPECT_THROW((void)VariantQC{colliding}, std::invalid_argument);
            EXPECT_EQ(read_text_file(target), content);
        }
    }
}

TEST_F(VariantQCH1Training, TrainingTraceMissingParentAndDirectoryDestinationsThrow) {
    auto config = make_variant_qc_config(path("input.vcf"), path("output.vcf"),
                                         path("output.tsv"));
    write_text_file(config.input_vcf, std::string(h1_vcf_header) + reference_rows("REF_ONLY"));
    for (const std::string suffix : {".tsv", ".tsv.gz", ".tsv.bgz"}) {
        SCOPED_TRACE(suffix);
        config.training_trace = path("missing/trace" + suffix);
        ASSERT_FALSE(std::filesystem::exists(path("missing")));
        {
            H1StderrCapture capture;
            EXPECT_THROW(VariantQC(config).run(), std::exception);
        }
        EXPECT_FALSE(std::filesystem::exists(path("missing")));
        config.training_trace = path("directory" + suffix);
        ASSERT_TRUE(std::filesystem::create_directory(config.training_trace));
        {
            H1StderrCapture capture;
            EXPECT_THROW(VariantQC(config).run(), std::exception);
        }
        EXPECT_TRUE(std::filesystem::is_directory(config.training_trace));
    }
}

TEST_F(VariantQCH1Training, TrainingTraceUnwritableParentThrows) {
    auto config = make_variant_qc_config(path("input.vcf"), path("output.vcf"),
                                         path("output.tsv"));
    write_text_file(config.input_vcf, std::string(h1_vcf_header) + reference_rows("REF_ONLY"));
    const auto locked = path("read_only");
    ASSERT_TRUE(std::filesystem::create_directory(locked));
    std::error_code error;
    std::filesystem::permissions(locked, std::filesystem::perms::owner_read |
                                std::filesystem::perms::owner_exec, error);
    ASSERT_FALSE(error) << error.message();
    bool can_write = false;
    {
        std::ofstream probe(path("read_only/probe"));
        can_write = probe.is_open();
    }
    if (!can_write) {
        for (const std::string suffix : {".tsv", ".tsv.gz", ".tsv.bgz"}) {
            SCOPED_TRACE(suffix);
            config.training_trace = path("read_only/trace" + suffix);
            H1StderrCapture capture;
            EXPECT_THROW(VariantQC(config).run(), std::exception);
            EXPECT_FALSE(std::filesystem::exists(config.training_trace));
        }
    }
    // Restore traversal/deletion permissions before fixture cleanup, including privileged runs.
    std::filesystem::permissions(locked, std::filesystem::perms::owner_all, error);
    EXPECT_FALSE(error) << error.message();
    if (can_write) GTEST_SKIP() << "Current privileges bypass the fixture's write permissions";
}

}  // namespace

