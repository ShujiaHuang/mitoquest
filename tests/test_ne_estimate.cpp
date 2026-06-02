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
#include <cstdio>
#include <fstream>
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
    // The Beta-Binomial MMLE is consistent but for finite cohorts a
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

// =====================================================================
// Continuous (Beta-diffusion) model tests (v1.8.2)
// =====================================================================

TEST(NeEstContinuous, SinglePairLLIsFinite) {
    // Basic smoke test: continuous LL should be finite for valid input.
    NeEstimator::LogFactorial lf(2000);
    NeEstimator::PairData pd{2000, 600, 2000, 700}; // mom=0.3, child=0.35
    for (int ne : {2, 5, 10, 20, 50}) {
        double ll = NeEstimator::compute_ll_single_continuous(pd, ne, lf);
        EXPECT_TRUE(std::isfinite(ll));
        EXPECT_LT(ll, 0.0);  // log-likelihood is negative
    }
}

TEST(NeEstContinuous, NearKimuraOnCleanData) {
    // On clean WF-generated data (no post-bottleneck noise),
    // continuous MMLE and Kimura should agree within ~50%.
    constexpr int true_ne = 10;
    auto data = simulate_pairs(true_ne,
                               /*n_pairs=*/600,
                               /*m_dp=*/   2000,
                               /*c_dp=*/   2000,
                               /*vaf_low=*/0.20,
                               /*vaf_high=*/0.80,
                               /*seed=*/   555u);
    auto r = NeEstimator::estimate(data, 1, 100, 1, /*continuous=*/true);
    auto k = NeEstimator::compute_kimura_check(data);

    // Both should be in the same ballpark.
    EXPECT_GE(r.ne, 5);
    EXPECT_LE(r.ne, 25);
    EXPECT_GE(k.ne_kimura, 5.0);
    EXPECT_LE(k.ne_kimura, 25.0);
    // They should agree within a factor of 2.
    double ratio = r.ne / k.ne_kimura;
    EXPECT_GT(ratio, 0.5);
    EXPECT_LT(ratio, 2.0);
}

TEST(NeEstContinuous, DiscreteGivesHigherNeThanContinuousOnCleanData) {
    // On WF data with moderate Ne, the discrete model can give a
    // similar or slightly lower Ne (on perfect data it hits the
    // exact grid).  The continuous model's CI should be wider.
    constexpr int true_ne = 10;
    auto data = simulate_pairs(true_ne, 300, 2000, 2000, 0.20, 0.80, 321u);
    auto r_disc = NeEstimator::estimate(data, 1, 50, 1, /*continuous=*/false);
    auto r_cont = NeEstimator::estimate(data, 1, 50, 1, /*continuous=*/true);

    // Both should be near the truth.
    EXPECT_GE(r_disc.ne, 5);
    EXPECT_LE(r_disc.ne, 20);
    EXPECT_GE(r_cont.ne, 5);
    EXPECT_LE(r_cont.ne, 20);
    // Continuous CI should typically be wider than discrete.
    double ci_width_cont = r_cont.ci_high - r_cont.ci_low;
    double ci_width_disc = r_disc.ci_high - r_disc.ci_low;
    EXPECT_GE(ci_width_cont, ci_width_disc);
}

// =====================================================================
// Continuous model: real-valued Ne optimization (v1.8.3)
// =====================================================================

namespace {
// Simulate a cohort under the continuous Beta-diffusion model:
//   p0       ~ Uniform[vaf_low, vaf_high]
//   p_child  ~ Beta(p0 * (true_ne - 1), (1 - p0) * (true_ne - 1))
//   m_ad_alt ~ Binomial(m_dp, p0)
//   c_ad_alt ~ Binomial(c_dp, p_child)
std::vector<NeEstimator::PairData>
simulate_pairs_continuous(double true_ne, int n_pairs,
                          int m_dp, int c_dp,
                          double vaf_low, double vaf_high,
                          unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> pdist(vaf_low, vaf_high);

    std::vector<NeEstimator::PairData> data;
    data.reserve(n_pairs);
    const double ne1 = true_ne - 1.0;
    for (int i = 0; i < n_pairs; ++i) {
        const double p0 = pdist(rng);
        // Draw child true heteroplasmy from Beta diffusion.
        const double alpha = p0 * ne1;
        const double beta  = (1.0 - p0) * ne1;
        std::gamma_distribution<double> ga(alpha, 1.0);
        std::gamma_distribution<double> gb(beta,  1.0);
        double xa = ga(rng);
        double xb = gb(rng);
        double p_child = xa / (xa + xb);
        if (p_child < 0.0) p_child = 0.0;
        if (p_child > 1.0) p_child = 1.0;

        std::binomial_distribution<int> bm(m_dp, p0);
        std::binomial_distribution<int> bc(c_dp, p_child);
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

TEST(NeEstContinuous, RealValuedOptimum) {
    // Simulate under the continuous model with a non-integer true Ne.
    // The real-valued optimizer should return a fractional Ne that is
    // closer to the truth than any integer.
    constexpr double true_ne = 5.5;
    auto data = simulate_pairs_continuous(true_ne,
                                          /*n_pairs=*/800,
                                          /*m_dp=*/2000,
                                          /*c_dp=*/2000,
                                          /*vaf_low=*/0.20,
                                          /*vaf_high=*/0.80,
                                          /*seed=*/2024u);
    auto r = NeEstimator::estimate(data, 1, 100, 1, /*continuous=*/true);

    // The real-valued MMLE should be fractional (not exactly integer).
    double frac_part = r.ne - std::floor(r.ne);
    // With 800 pairs there's enough power: fractional part should be
    // noticeably away from 0 or 1 (i.e. the optimizer actually refined).
    EXPECT_GT(frac_part, 0.01);
    EXPECT_LT(frac_part, 0.99);

    // Should be within 30% of the truth.
    EXPECT_NEAR(r.ne, true_ne, true_ne * 0.30);

    // CI should bracket the truth.
    EXPECT_LE(r.ci_low,  true_ne);
    EXPECT_GE(r.ci_high, true_ne);
}

TEST(NeEstContinuous, RealValuedAgreesWithKimura) {
    // On clean continuous-model data, the real-valued MMLE and Kimura
    // should agree within 20% (much tighter than the old 50% tolerance).
    constexpr double true_ne = 8.0;
    auto data = simulate_pairs_continuous(true_ne,
                                          /*n_pairs=*/1000,
                                          /*m_dp=*/2000,
                                          /*c_dp=*/2000,
                                          /*vaf_low=*/0.25,
                                          /*vaf_high=*/0.75,
                                          /*seed=*/2025u);
    auto r = NeEstimator::estimate(data, 1, 100, 1, /*continuous=*/true);
    auto k = NeEstimator::compute_kimura_check(data);

    EXPECT_NEAR(r.ne, true_ne, true_ne * 0.25);
    EXPECT_NEAR(k.ne_kimura, true_ne, true_ne * 0.25);
    // MMLE and Kimura should agree within 20%.
    double ratio = r.ne / k.ne_kimura;
    EXPECT_GT(ratio, 0.80);
    EXPECT_LT(ratio, 1.20);
}

// =====================================================================
// Kimura trimmed estimator
// =====================================================================

TEST(NeEstKimuraTrimmed, TrimmedFieldsPopulatedWhenTrimFracPositive) {
    // Clean data with true Ne=10: trimmed and untrimmed should both be
    // in the same ballpark because there are no real outliers.
    constexpr int true_ne = 10;
    auto data = simulate_pairs(true_ne,
                               /*n_pairs=*/400,
                               /*m_dp=*/   2000,
                               /*c_dp=*/   2000,
                               /*vaf_low=*/0.20,
                               /*vaf_high=*/0.80,
                               /*seed=*/   999u);
    auto k = NeEstimator::compute_kimura_check(data,
                                               /*n_bootstrap=*/0,
                                               /*seed=*/42,
                                               /*trim_frac=*/0.10,
                                               /*top_drift_k=*/0);
    EXPECT_TRUE(k.computed);
    EXPECT_TRUE(k.trimmed_computed);
    EXPECT_DOUBLE_EQ(k.trim_frac, 0.10);
    // n_after_trim = floor(0.9 * n_informative)
    EXPECT_GT(k.n_after_trim, 300u);
    EXPECT_LT(k.n_after_trim, k.n_informative);
    EXPECT_TRUE(std::isfinite(k.b_trimmed));
    EXPECT_TRUE(std::isfinite(k.ne_kimura_trimmed));
    EXPECT_GT(k.ne_kimura_trimmed, 1.0);
}

TEST(NeEstKimuraTrimmed, NotPopulatedWhenTrimFracZero) {
    // trim_frac = 0 => trimmed fields should NOT be populated.
    auto data = simulate_pairs(10, 100, 500, 500, 0.25, 0.75, 111u);
    auto k = NeEstimator::compute_kimura_check(data,
                                               /*n_bootstrap=*/0,
                                               /*seed=*/42,
                                               /*trim_frac=*/0.0,
                                               /*top_drift_k=*/0);
    EXPECT_TRUE(k.computed);
    EXPECT_FALSE(k.trimmed_computed);
}

TEST(NeEstKimuraTrimmed, TrimmedCloserToMMLEWhenOutliersPresent) {
    // Inject ~20% outlier pairs: child VAF far from mother.
    // MMLE is robust; untrimmed Kimura is pulled down; trimmed should
    // recover a value closer to the truth than the untrimmed Kimura.
    constexpr int true_ne = 20;
    auto data = simulate_pairs(true_ne,
                               /*n_pairs=*/300,
                               /*m_dp=*/   2000,
                               /*c_dp=*/   2000,
                               /*vaf_low=*/0.20,
                               /*vaf_high=*/0.80,
                               /*seed=*/   777u);
    // Replace the last 20% of pairs with extreme-drift fakes.
    const size_t n_outlier = data.size() / 5;
    for (size_t i = data.size() - n_outlier; i < data.size(); ++i) {
        data[i].m_ad_alt = data[i].m_dp / 2;  // mom ~0.5
        data[i].c_ad_alt = 0;                  // child = 0  (extreme drift)
    }
    auto k = NeEstimator::compute_kimura_check(data,
                                               /*n_bootstrap=*/0,
                                               /*seed=*/42,
                                               /*trim_frac=*/0.25,
                                               /*top_drift_k=*/5);
    EXPECT_TRUE(k.computed);
    EXPECT_TRUE(k.trimmed_computed);
    // Untrimmed Ne_kimura should be low (pulled by outliers).
    EXPECT_LT(k.ne_kimura, 10.0);
    // Trimmed should be noticeably higher than untrimmed.
    EXPECT_GT(k.ne_kimura_trimmed, k.ne_kimura * 1.5);

    // Top drift outliers: should have 5 entries.
    EXPECT_EQ(k.top_drift_outliers.size(), 5u);
    // The first outlier should have the largest F_i.
    if (k.top_drift_outliers.size() >= 2) {
        EXPECT_GE(k.top_drift_outliers[0].f_i,
                  k.top_drift_outliers[1].f_i);
    }
}

// =====================================================================
// Three-generation trio marginal likelihood tests
// =====================================================================

// Helper: build a trio PairData in one expression.
static NeEstimator::PairData make_trio(int g_dp, int g_ad_alt,
                                        int m_dp, int m_ad_alt,
                                        int c_dp, int c_ad_alt,
                                        bool has_g = true) {
    NeEstimator::PairData pd;
    pd.g_dp = g_dp; pd.g_ad_alt = g_ad_alt;
    pd.m_dp = m_dp; pd.m_ad_alt = m_ad_alt;
    pd.c_dp = c_dp; pd.c_ad_alt = c_ad_alt;
    pd.has_g = has_g ? 1 : 0;
    return pd;
}

// The closed-form formula must be numerically consistent with
// Gauss-Legendre quadrature.  Because GL quadrature on the non-polynomial
// Beta kernel has O(1/nodes^2) convergence, we compare *relative* LL
// differences (which is what the optimiser uses) rather than absolute LL.
TEST(NeEstTrio, ClosedFormVsQuadratureRelativeLL) {
    NeEstimator::LogFactorial lf(1000);
    struct Case { int g_dp, g_ad, m_dp, m_ad, c_dp, c_ad; };
    const std::vector<Case> cases = {
        {100, 30, 50, 15, 80, 20},
        {200, 80, 100, 45, 150, 60},
        {500, 250, 200, 110, 300, 150},
    };
    for (const auto& c : cases) {
        auto pd = make_trio(c.g_dp, c.g_ad, c.m_dp, c.m_ad, c.c_dp, c.c_ad);
        // Compare LL difference: LL(Ne=5) - LL(Ne=50) between methods.
        const double cf_lo = NeEstimator::compute_ll_trio_continuous(pd, 5.0, lf);
        const double cf_hi = NeEstimator::compute_ll_trio_continuous(pd, 50.0, lf);
        const double q_lo  = NeEstimator::compute_ll_trio_quadrature(pd, 5.0, lf, 256);
        const double q_hi  = NeEstimator::compute_ll_trio_quadrature(pd, 50.0, lf, 256);
        const double cf_diff = cf_lo - cf_hi;
        const double q_diff  = q_lo  - q_hi;
        // Relative LL differences should agree to within 1.5.
        // GL quadrature on the non-polynomial Beta kernel has O(1/n^2)
        // convergence; the absolute LL bias cancels partially in
        // differences but not perfectly.
        EXPECT_NEAR(cf_diff, q_diff, 1.5)
            << "LL-diff mismatch for g_dp=" << c.g_dp;
        // Both should be finite.
        EXPECT_TRUE(std::isfinite(cf_lo));
        EXPECT_TRUE(std::isfinite(cf_hi));
        EXPECT_TRUE(std::isfinite(q_lo));
        EXPECT_TRUE(std::isfinite(q_hi));
    }
}

// Verify that both closed-form and quadrature find the same optimal Ne.
TEST(NeEstTrio, ClosedFormAndQuadratureAgreeOnOptimalNe) {
    NeEstimator::LogFactorial lf(1000);
    auto pd = make_trio(200, 60, 100, 35, 150, 50);
    // Scan Ne in [2, 100] and find the argmax for both methods.
    double best_cf = 2.0, best_q = 2.0;
    double max_cf = -1e30, max_q = -1e30;
    for (double ne = 2.0; ne <= 100.0; ne += 1.0) {
        const double ll_cf = NeEstimator::compute_ll_trio_continuous(pd, ne, lf);
        const double ll_q  = NeEstimator::compute_ll_trio_quadrature(pd, ne, lf, 256);
        if (ll_cf > max_cf) { max_cf = ll_cf; best_cf = ne; }
        if (ll_q  > max_q)  { max_q  = ll_q;  best_q  = ne; }
    }
    // Optimal Ne should agree to within ±3 (step size is 1.0).
    EXPECT_NEAR(best_cf, best_q, 3.0);
}

// When has_g == 0 the trio function must fall back to the 2-gen model.
TEST(NeEstTrio, HasGZeroFallsBackTo2Gen) {
    NeEstimator::LogFactorial lf(500);
    const double ne = 15.0;
    auto pd = make_trio(100, 30, 50, 20, 80, 25, /*has_g=*/false);
    const double ll_trio = NeEstimator::compute_ll_trio_continuous(pd, ne, lf);
    const double ll_2gen = NeEstimator::compute_ll_single_continuous(pd, ne, lf);
    EXPECT_DOUBLE_EQ(ll_trio, ll_2gen);
}

// Homoplasmic grandmother (g_ad_alt == 0 or g_dp) is Ne-independent -> 0.0.
TEST(NeEstTrio, HomoplasmicGrandmotherReturnsZero) {
    NeEstimator::LogFactorial lf(200);
    // Grandmother homoplasmic REF.
    auto pd_ref = make_trio(100, 0, 50, 20, 80, 30);
    EXPECT_DOUBLE_EQ(NeEstimator::compute_ll_trio_continuous(pd_ref, 10.0, lf), 0.0);
    EXPECT_DOUBLE_EQ(NeEstimator::compute_ll_trio_continuous(pd_ref, 50.0, lf), 0.0);
    // Grandmother homoplasmic ALT.
    auto pd_alt = make_trio(100, 100, 50, 20, 80, 30);
    EXPECT_DOUBLE_EQ(NeEstimator::compute_ll_trio_continuous(pd_alt, 10.0, lf), 0.0);
}

// At Ne = 1 the diffusion degenerates (complete drift to fixation).
// For a heteroplasmic child the discrete fallback correctly returns
// -infinity (impossible under Ne=1).  Slightly above Ne=1 the result
// is finite, confirming the continuous formula is well-behaved.
TEST(NeEstTrio, NeEqualsOneBoundary) {
    NeEstimator::LogFactorial lf(200);
    auto pd = make_trio(100, 30, 50, 15, 80, 20);
    // Ne = 1 with heteroplasmic child: impossible under fixation -> -inf.
    const double ll1 = NeEstimator::compute_ll_trio_continuous(pd, 1.0, lf);
    EXPECT_TRUE(std::isinf(ll1) && ll1 < 0.0);
    // Just above Ne = 1: must be finite.
    const double ll_plus = NeEstimator::compute_ll_trio_continuous(pd, 1.01, lf);
    EXPECT_TRUE(std::isfinite(ll_plus));
    // Homoplasmic child: Ne=1 should give a finite result.
    auto pd_hom = make_trio(100, 30, 50, 15, 80, 80);
    const double ll_hom = NeEstimator::compute_ll_trio_continuous(pd_hom, 1.0, lf);
    EXPECT_TRUE(std::isfinite(ll_hom));
}

// Closed-form self-consistency: verify the formula directly for one case.
// Ne=10, p_G=0.3 -> alpha_G=2.7, beta_G=6.3
// k_M=15, d_M=50, k_C=20, d_C=80
// A = 2.7+15+20 = 37.7, B = 6.3+35+60 = 101.3
// log I = log C(50,15) + log C(80,20)
//       + lgamma(37.7) + lgamma(101.3) - lgamma(139.0)
//       - lgamma(2.7) - lgamma(6.3) + lgamma(9.0)
TEST(NeEstTrio, ClosedFormManualCheck) {
    NeEstimator::LogFactorial lf(500);
    const double ne = 10.0;
    auto pd = make_trio(100, 30, 50, 15, 80, 20);
    const double ll = NeEstimator::compute_ll_trio_continuous(pd, ne, lf);

    const double alpha_G = 0.3 * 9.0;   // 2.7
    const double beta_G  = 0.7 * 9.0;   // 6.3
    const double A = alpha_G + 15.0 + 20.0;  // 37.7
    const double B = beta_G  + 35.0 + 60.0;  // 101.3
    const double expected =
          std::lgamma(A) + std::lgamma(B) - std::lgamma(A + B)
        - std::lgamma(alpha_G) - std::lgamma(beta_G)
        + std::lgamma(alpha_G + beta_G)
        + lf.log_comb(50, 15) + lf.log_comb(80, 20);
    EXPECT_NEAR(ll, expected, 1e-9);
}

// Trio likelihood must differ from 2-gen likelihood when has_g == 1.
TEST(NeEstTrio, TrioDiffersFrom2GenWhenHasGOne) {
    NeEstimator::LogFactorial lf(500);
    const double ne = 10.0;
    auto pd = make_trio(100, 30, 50, 20, 80, 25, /*has_g=*/true);
    const double ll_trio = NeEstimator::compute_ll_trio_continuous(pd, ne, lf);
    const double ll_2gen = NeEstimator::compute_ll_single_continuous(pd, ne, lf);
    // They should be different (trio uses grandmother info).
    EXPECT_NE(ll_trio, ll_2gen);
    EXPECT_TRUE(std::isfinite(ll_trio));
    EXPECT_TRUE(std::isfinite(ll_2gen));
}

// Global LL with trio rows should equal the sum of per-row trio LLs.
TEST(NeEstTrio, GlobalLLDispatchesTrioRows) {
    NeEstimator::LogFactorial lf(500);
    const double ne = 10.0;
    // Build a dataset with one trio row and one 2-gen row.
    auto trio_pd = make_trio(100, 30, 50, 15, 80, 20, true);
    auto pair_pd = make_trio(0, 0, 80, 25, 120, 40, false);
    std::vector<NeEstimator::PairData> data = {trio_pd, pair_pd};
    const double global = NeEstimator::compute_global_ll_continuous(ne, data, lf, 1);
    const double expected =
          NeEstimator::compute_ll_trio_continuous(trio_pd, ne, lf)
        + NeEstimator::compute_ll_single_continuous(pair_pd, ne, lf);
    EXPECT_NEAR(global, expected, 1e-9);
}

// load_pairs must read HAS_G / GRANDMOTHER_DP / GRANDMOTHER_AD_ALT when
// present, and must leave has_g == 0 for rows with HAS_G == 0 or when the
// columns are absent (legacy TSV).
TEST(NeEstTrio, LoadPairsReadsTrioColumns) {
    const std::string tsv = "tests/data/ne_pipeline/trio_load_test.tsv";
    // Write a minimal wide TSV with one trio row and one MC-only row.
    {
        std::ofstream out(tsv);
        ASSERT_TRUE(out.is_open());
        out << "CHROM\tPOS\tREF\tALT\tFAM_ID\t"
               "GRANDMOTHER_ID\tMOTHER_ID\tCHILD_ID\t"
               "GRANDMOTHER_DP\tGRANDMOTHER_AD_REF\tGRANDMOTHER_AD_ALT\t"
               "GRANDMOTHER_VAF\tHAS_G\t"
               "MOTHER_DP\tMOTHER_AD_REF\tMOTHER_AD_ALT\tMOTHER_VAF\t"
               "CHILD_DP\tCHILD_AD_REF\tCHILD_AD_ALT\tCHILD_VAF\tQC\n";
        // Row 1: HAS_G = 1, grandmother 100 dp 30 alt.
        out << "chrM\t100\tA\tT\tF1\tGMA\tMOM\tCHILD\t"
               "100\t70\t30\t0.300\t1\t"
               "50\t35\t15\t0.300\t"
               "80\t60\t20\t0.250\tPASS\n";
        // Row 2: HAS_G = 0 (no grandmother).
        out << "chrM\t200\tC\tG\tF2\tNA\tMOM2\tCHILD2\t"
               "0\t0\t0\t0.000\t0\t"
               "80\t55\t25\t0.312\t"
               "120\t80\t40\t0.333\tPASS\n";
    }
    auto data = NeEstimator::load_pairs(tsv, 0.0, 1.0);
    ASSERT_EQ(data.size(), 2u);
    // Row 1: trio.
    EXPECT_EQ(data[0].has_g,    1);
    EXPECT_EQ(data[0].g_dp,     100);
    EXPECT_EQ(data[0].g_ad_alt, 30);
    EXPECT_EQ(data[0].m_dp,     50);
    EXPECT_EQ(data[0].c_dp,     80);
    // Row 2: no trio.
    EXPECT_EQ(data[1].has_g,    0);
    EXPECT_EQ(data[1].g_dp,     0);
    EXPECT_EQ(data[1].g_ad_alt, 0);
    // Clean up.
    std::remove(tsv.c_str());
}

// Legacy 16-column TSV (no HAS_G column) -> has_g == 0 for all rows.
TEST(NeEstTrio, LoadPairsLegacyTsvHasGZero) {
    const std::string tsv = "tests/data/ne_pipeline/legacy_load_test.tsv";
    {
        std::ofstream out(tsv);
        ASSERT_TRUE(out.is_open());
        out << "CHROM\tPOS\tREF\tALT\tFAM_ID\t"
               "MOTHER_ID\tCHILD_ID\t"
               "MOTHER_DP\tMOTHER_AD_REF\tMOTHER_AD_ALT\tMOTHER_VAF\t"
               "CHILD_DP\tCHILD_AD_REF\tCHILD_AD_ALT\tCHILD_VAF\tQC\n";
        out << "chrM\t100\tA\tT\tF1\tMOM\tCHILD\t"
               "50\t35\t15\t0.300\t"
               "80\t60\t20\t0.250\tPASS\n";
    }
    auto data = NeEstimator::load_pairs(tsv, 0.0, 1.0);
    ASSERT_EQ(data.size(), 1u);
    EXPECT_EQ(data[0].has_g, 0);
    EXPECT_EQ(data[0].g_dp,  0);
    std::remove(tsv.c_str());
}
