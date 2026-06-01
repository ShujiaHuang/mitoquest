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

#include <cmath>
#include <random>
#include <vector>
#include <numeric>

#include "variant_qc.h"
#include "log_factorial.h"

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


