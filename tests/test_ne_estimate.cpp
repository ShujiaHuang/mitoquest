// Author: Shujia Huang
// Date: 2026-05-28
//
// Unit tests for the `mitoquest ne-estimate` subcommand
// (src/ne_estimate.{h,cpp}).
//
// Coverage:
//   * LogFactorial: log_fact, log_comb, log_binomial_pmf, log_betabinom_pmf
//     (closed-form spot checks against analytic values).
//   * log_sum_exp_pair handles -INFINITY sentinel correctly.
//   * compute_ll_single agrees with brute-force enumeration at small Ne.
//   * required_cache_size grows with the deepest input pair / search ceiling.
//   * find_optimal_ne / estimate recover the true Ne on a deterministic
//     synthetic dataset (seeded simulation, true Ne = 10).

#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <vector>

#include "log_factorial.h"
#include "ne_estimate.h"

// =====================================================================
// LogFactorial closed-form spot checks
// =====================================================================

TEST(NeEstLogFactorial, LogFactSmallValues) {
    NeEstimator::LogFactorial lf(20);
    // 0! = 1 -> log = 0
    EXPECT_NEAR(lf.log_fact(0), 0.0, 1e-12);
    // 5! = 120 -> log(120)
    EXPECT_NEAR(lf.log_fact(5), std::log(120.0), 1e-9);
    // Spillover beyond cache: should fall back to lgamma.
    EXPECT_NEAR(lf.log_fact(50),
                std::lgamma(51.0),
                1e-9);
}

TEST(NeEstLogFactorial, LogCombBoundaries) {
    NeEstimator::LogFactorial lf(20);
    // C(10, 0) = C(10, 10) = 1 -> log = 0
    EXPECT_NEAR(lf.log_comb(10, 0),  0.0, 1e-12);
    EXPECT_NEAR(lf.log_comb(10, 10), 0.0, 1e-12);
    // C(10, 3) = 120
    EXPECT_NEAR(lf.log_comb(10, 3), std::log(120.0), 1e-9);
    // Out-of-range: -infinity
    EXPECT_TRUE(std::isinf(lf.log_comb(5, -1)));
    EXPECT_LT(lf.log_comb(5, -1), 0.0);
    EXPECT_TRUE(std::isinf(lf.log_comb(5, 6)));
    EXPECT_LT(lf.log_comb(5, 6), 0.0);
}

TEST(NeEstLogFactorial, BinomialPmfBoundaries) {
    NeEstimator::LogFactorial lf(20);
    // P(k=0 | n=10, p=0) = 1
    EXPECT_NEAR(lf.log_binomial_pmf(10, 0, 0.0), 0.0, 1e-12);
    // P(k=10 | n=10, p=1) = 1
    EXPECT_NEAR(lf.log_binomial_pmf(10, 10, 1.0), 0.0, 1e-12);
    // P(k=5 | n=10, p=0) = 0
    EXPECT_TRUE(std::isinf(lf.log_binomial_pmf(10, 5, 0.0)));
    EXPECT_LT(lf.log_binomial_pmf(10, 5, 0.0), 0.0);
    // P(k=2 | n=4, p=0.5) = 6 * 0.5^4 = 0.375
    EXPECT_NEAR(std::exp(lf.log_binomial_pmf(4, 2, 0.5)), 0.375, 1e-9);
}

// Analytic Beta-Binomial PMF for a few small values:
//   P(k|n=2, alpha=1, beta=1) = 1/(n+1) for all k in [0..n].
TEST(NeEstLogFactorial, BetaBinomUniformPrior) {
    NeEstimator::LogFactorial lf(20);
    // n=2, uniform prior alpha=beta=1 -> uniform over k.
    for (int k = 0; k <= 2; ++k) {
        EXPECT_NEAR(std::exp(lf.log_betabinom_pmf(2, k, 1.0, 1.0)),
                    1.0 / 3.0, 1e-9);
    }
    // PMF must sum to 1 for general (alpha, beta).
    {
        double s = 0.0;
        for (int k = 0; k <= 5; ++k) {
            s += std::exp(lf.log_betabinom_pmf(5, k, 2.0, 3.0));
        }
        EXPECT_NEAR(s, 1.0, 1e-9);
    }
}

// =====================================================================
// log-sum-exp sentinel handling
// =====================================================================

TEST(NeEstLogSumExp, NegInfinitySentinel) {
    const double neg_inf = -std::numeric_limits<double>::infinity();
    EXPECT_NEAR(NeEstimator::log_sum_exp_pair(neg_inf, 0.0), 0.0, 1e-12);
    EXPECT_NEAR(NeEstimator::log_sum_exp_pair(0.0, neg_inf), 0.0, 1e-12);
    // log(exp(0)+exp(0)) = log(2)
    EXPECT_NEAR(NeEstimator::log_sum_exp_pair(0.0, 0.0),
                std::log(2.0), 1e-12);
}

// =====================================================================
// compute_ll_single: brute-force agreement at small Ne
// =====================================================================
//
// Reimplement the joint PMF directly from definitions and compare.
namespace {
double brute_force_ll(const NeEstimator::PairData& pd, int ne,
                      const NeEstimator::LogFactorial& lf) {
    const double alpha = static_cast<double>(pd.m_ad_alt) + 1.0;
    const double beta  = static_cast<double>(pd.m_dp - pd.m_ad_alt) + 1.0;
    double total = 0.0;
    for (int k = 0; k <= ne; ++k) {
        const double bb  = std::exp(lf.log_betabinom_pmf(ne, k, alpha, beta));
        const double p1  = static_cast<double>(k) / static_cast<double>(ne);
        const double bin = std::exp(lf.log_binomial_pmf(pd.c_dp, pd.c_ad_alt, p1));
        total += bb * bin;
    }
    return std::log(total);
}
}  // namespace

TEST(NeEstComputeLLSingle, AgreesWithBruteForce) {
    NeEstimator::LogFactorial lf(50);

    // Pair: mother 100 reads, 40 ALT  -> alpha=41, beta=61
    //       child   80 reads, 30 ALT
    NeEstimator::PairData pd{100, 40, 80, 30};

    for (int ne : {2, 3, 5, 8}) {
        const double a = NeEstimator::compute_ll_single(pd, ne, lf);
        const double b = brute_force_ll(pd, ne, lf);
        EXPECT_NEAR(a, b, 1e-8) << "Ne=" << ne;
    }
}

TEST(NeEstComputeLLSingle, BoundaryPairs) {
    NeEstimator::LogFactorial lf(50);
    // Mother is fixed homozygous ALT (alpha large, beta=1), child homozygous ALT.
    // The likelihood should peak at any Ne >= 1.
    NeEstimator::PairData pd_homo{50, 50, 50, 50};
    EXPECT_GT(NeEstimator::compute_ll_single(pd_homo, 1, lf), -10.0);
    EXPECT_GT(NeEstimator::compute_ll_single(pd_homo, 5, lf), -10.0);

    // Mother ~50/50 heteroplasmic, child homozygous: Ne=1 is more likely
    // than Ne >> 1 because a 1-read bottleneck cleanly explains a child
    // hitting the homozygous ALT outcome with probability ~0.5.
    NeEstimator::PairData pd_split{1000, 500, 500, 500};
    const double ll1   = NeEstimator::compute_ll_single(pd_split, 1,   lf);
    const double ll100 = NeEstimator::compute_ll_single(pd_split, 100, lf);
    EXPECT_GT(ll1, ll100);
}

// =====================================================================
// required_cache_size
// =====================================================================

TEST(NeEstRequiredCacheSize, GrowsWithDpAndMaxNe) {
    std::vector<NeEstimator::PairData> data = {
        {1000, 200,  500, 100},
        { 800, 100, 1500, 300},   // child DP=1500 dominates
        { 600,  60,  900, 100},
    };
    EXPECT_EQ(NeEstimator::required_cache_size(data, 50),  1500);
    EXPECT_EQ(NeEstimator::required_cache_size(data, 5000), 5000);
    EXPECT_EQ(NeEstimator::required_cache_size({}, 200),    200);
    EXPECT_GE(NeEstimator::required_cache_size({}, 0),      1);
}

// =====================================================================
// find_optimal_ne / estimate: recover true Ne from a synthetic cohort
// =====================================================================

namespace {
// Simulate a cohort of (mother, child) pairs assuming
//   p0 ~ Uniform[low, high]
//   k  ~ Binomial(true_ne, p0)
//   m_ad_alt ~ Binomial(m_dp, p0)
//   c_ad_alt ~ Binomial(c_dp, k / true_ne)
std::vector<NeEstimator::PairData>
simulate_pairs(int true_ne, int n_pairs,
               int m_dp, int c_dp,
               double vaf_low, double vaf_high,
               unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> pdist(vaf_low, vaf_high);

    std::vector<NeEstimator::PairData> data;
    data.reserve(n_pairs);
    for (int i = 0; i < n_pairs; ++i) {
        const double p0 = pdist(rng);
        std::binomial_distribution<int> bk(true_ne, p0);
        const int k = bk(rng);

        std::binomial_distribution<int> bm(m_dp, p0);
        std::binomial_distribution<int> bc(c_dp,
                                           static_cast<double>(k) / true_ne);
        NeEstimator::PairData pd;
        pd.m_dp     = m_dp;
        pd.m_ad_alt = bm(rng);
        pd.c_dp     = c_dp;
        pd.c_ad_alt = bc(rng);
        data.push_back(pd);
    }
    return data;
}
}  // namespace

TEST(NeEstFindOptimal, RecoversTrueNeOnSyntheticCohort) {
    constexpr int true_ne = 10;
    auto data = simulate_pairs(true_ne,
                               /*n_pairs=*/200,
                               /*m_dp=*/   1000,
                               /*c_dp=*/   1000,
                               /*vaf_low=*/0.20,
                               /*vaf_high=*/0.80,
                               /*seed=*/   42u);

    const int cache = NeEstimator::required_cache_size(data, 80);
    NeEstimator::LogFactorial lf(cache);
    const int est = NeEstimator::find_optimal_ne(data, lf,
                                                 /*min_ne=*/1,
                                                 /*max_ne=*/80);
    // The Beta-Binomial MLE is consistent but for finite cohorts a
    // 30% relative band around the truth is the expected regime.
    EXPECT_GE(est, 7);
    EXPECT_LE(est, 14);
}

TEST(NeEstEstimate, ReportsCIBracketingTruth) {
    constexpr int true_ne = 10;
    auto data = simulate_pairs(true_ne,
                               /*n_pairs=*/300,
                               /*m_dp=*/   2000,
                               /*c_dp=*/   2000,
                               /*vaf_low=*/0.20,
                               /*vaf_high=*/0.80,
                               /*seed=*/   7u);

    NeEstimator::Result r = NeEstimator::estimate(data,
                                                  /*min_ne=*/1,
                                                  /*max_ne=*/100,
                                                  /*threads=*/1);
    EXPECT_EQ(r.n_pairs, data.size());
    EXPECT_GE(r.ne, 1);
    EXPECT_LE(r.ci_low,  r.ne);
    EXPECT_LE(r.ne,      r.ci_high);

    // 95% profile-likelihood CI should contain the true Ne in this regime.
    EXPECT_LE(r.ci_low,  true_ne);
    EXPECT_GE(r.ci_high, true_ne);

    // Boundary clip flags should not fire when truth lies safely inside
    // the [1, 100] search window.
    EXPECT_FALSE(r.ci_low_clipped);
    EXPECT_FALSE(r.ci_high_clipped);
}

TEST(NeEstEstimate, ThrowsOnEmptyInput) {
    EXPECT_THROW(NeEstimator::estimate({}, 1, 50, 1), std::runtime_error);
}

// =====================================================================
// Standalone mitoquest::LogFactorial (after refactor out of NeEstimator)
// =====================================================================

TEST(LogFactorialStandalone, IdenticalToNeEstimatorAlias) {
    // The NeEstimator::LogFactorial type-alias must point at exactly
    // the same class.
    static_assert(std::is_same<NeEstimator::LogFactorial,
                               mitoquest::LogFactorial>::value,
                  "NeEstimator::LogFactorial must alias mitoquest::LogFactorial");

    mitoquest::LogFactorial lf(20);
    EXPECT_NEAR(lf.log_fact(0), 0.0, 1e-12);
    EXPECT_NEAR(lf.log_fact(5), std::log(120.0), 1e-9);
    EXPECT_NEAR(std::exp(lf.log_betabinom_pmf(2, 1, 1.0, 1.0)), 1.0 / 3.0, 1e-9);
}

// =====================================================================
// Wonnapinij / Kimura cross-check
// =====================================================================

TEST(NeEstKimura, RecoversBOnLargeNeCohort) {
    // For a very large bottleneck Ne, the per-generation drift variance
    // V = p_m (1 - p_m) / Ne is tiny, so b ~ 1 and Ne_kimura should be
    // large.
    constexpr int true_ne = 50;
    auto data = simulate_pairs(true_ne,
                               /*n_pairs=*/600,
                               /*m_dp=*/   2000,
                               /*c_dp=*/   2000,
                               /*vaf_low=*/0.30,
                               /*vaf_high=*/0.70,
                               /*seed=*/   123u);
    auto k = NeEstimator::compute_kimura_check(data);
    EXPECT_TRUE(k.computed);
    EXPECT_GT(k.n_informative, 500u);
    EXPECT_GT(k.b, 0.5);                 // far from the boundary
    EXPECT_GT(k.ne_kimura, 5.0);         // reasonable single-generation Ne
}

TEST(NeEstKimura, ReportsSmallNeOnTightBottleneck) {
    // True Ne = 2 produces large drift variance; b should be small and
    // the implied Ne_kimura should be of order O(1)-O(10).
    constexpr int true_ne = 2;
    auto data = simulate_pairs(true_ne,
                               /*n_pairs=*/800,
                               /*m_dp=*/   2000,
                               /*c_dp=*/   2000,
                               /*vaf_low=*/0.30,
                               /*vaf_high=*/0.70,
                               /*seed=*/   321u);
    auto k = NeEstimator::compute_kimura_check(data);
    EXPECT_TRUE(k.computed);
    EXPECT_GT(k.n_informative, 600u);
    // Wonnapinij is approximate, but on synthetic data the implied Ne
    // should land in [1, 10] for true Ne = 2.
    EXPECT_GE(k.ne_kimura, 1.0);
    EXPECT_LE(k.ne_kimura, 10.0);
}

TEST(NeEstKimura, EmptyOrHomoplasmicReturnsInformativeFlag) {
    // Empty cohort: still computed=true, but n_informative=0.
    auto k_empty = NeEstimator::compute_kimura_check({});
    EXPECT_TRUE(k_empty.computed);
    EXPECT_EQ(k_empty.n_informative, 0u);
    EXPECT_FALSE(k_empty.note.empty());

    // All mothers homoplasmic ALT: weight w = p(1-p) = 0; nothing
    // contributes, n_informative = 0.
    std::vector<NeEstimator::PairData> homo = {
        {1000, 1000, 1000, 1000},
        { 500,  500,  500,  500},
    };
    auto k_homo = NeEstimator::compute_kimura_check(homo);
    EXPECT_TRUE(k_homo.computed);
    EXPECT_EQ(k_homo.n_informative, 0u);
}
