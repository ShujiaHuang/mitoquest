/**
 * @file ne_estimate.h
 * @brief Estimate the mtDNA bottleneck size (Ne) from mother-child
 *        transmission pairs via maximum likelihood, with an optional
 *        Wonnapinij/Kimura cross-check.
 *
 * ====================================================================
 * Statistical model (primary estimator: Beta-Binomial MLE)
 * ====================================================================
 *
 * For each independent mother-child (M-C) transmission pair we treat
 * the maternal heteroplasmy `p0` as latent and *fit it from the actual
 * read counts*, rather than substituting the VCF point estimate.
 *
 *   Maternal posterior :  p0 ~ Beta(alpha = m_alt + 1,
 *                                   beta  = m_ref + 1)
 *                         (uniform Beta(1, 1) prior on heteroplasmy.)
 *
 *   Single-generation
 *   bottleneck         :  k  ~ Binomial(Ne, p0)
 *   Marginalising p0   :  k  ~ BetaBinomial(Ne, alpha, beta)
 *
 *   Child read counts  :  c_alt ~ Binomial(c_dp, p1),   p1 = k / Ne.
 *
 *   Per-pair logL(Ne)  =  log Sigma_{k=0..Ne} [
 *                              logBetaBin(k | Ne, alpha, beta)
 *                            + logBin     (c_alt | c_dp, k/Ne) ]
 *
 *   Global logL(Ne)    =  Sigma over independent M-C pairs.
 *
 * Optimum is found by a *brute-force integer scan* over [min_ne, max_ne].
 * The discrete LL is NOT unimodal in Ne — when the true Ne is small,
 * integer multiples (2*Ne, 3*Ne, ...) of the truth align with the
 * observed VAF grid and create secondary local maxima.  Golden-section
 * search converges to those side maxima and misses the true peak; a
 * full scan is O(max_ne) global-LL evaluations and is plenty fast for
 * the typical operating range max_ne <= a few hundred.
 *
 * 95% confidence interval is the *contiguous* Ne range with
 *   logL(Ne) >= logL_max - chi2_{1, 0.95} / 2  =  logL_max - 1.92
 * (Wilks's theorem on the profile log-likelihood; see Cox & Hinkley
 *  1974, "Theoretical Statistics", Section 9.3).
 *
 * ====================================================================
 * Why DP/AD instead of the VCF AF field?
 * ====================================================================
 *
 *   * VCF AF is a *point estimate* and discards depth-dependent
 *     uncertainty.  A homoplasmic site at depth 50 (AD=50/50, AF=1.0)
 *     and at depth 5,000 (AD=5000/5000, AF=1.0) carry very different
 *     amounts of information about the true maternal heteroplasmy.
 *   * The Beta-Binomial posterior naturally weights each pair by its
 *     depth and is the *exact* discrete-Wright-Fisher likelihood for a
 *     single transmission.  Any approach that plugs in AF as if it
 *     were the truth biases Ne (downward at low VAF, upward near 0.5)
 *     and underestimates the CI width.
 *
 * ====================================================================
 * Cross-check estimator (Wonnapinij / Kimura), optional
 * ====================================================================
 *
 * As an independent sanity check we also compute the Wonnapinij
 * bottleneck parameter `b` using the sampling-error-corrected
 * method-of-moments estimator (Wonnapinij et al., 2008/2010; Helgason
 * et al., 2024 Cell), then convert to a single-generation Ne via
 *   Ne_kimura = 1 / (1 - b).
 * This is reported alongside the Beta-Binomial MLE for comparison
 * with the deCODE genetics 2024 Cell paper, but it is NOT the primary
 * estimator: for single-generation M-C data the Beta-Binomial MLE is
 * the exact form that the Kimura distribution approximates and it is
 * strictly more efficient (Cramer-Rao).
 *
 * @author Shujia Huang (hshujia@qq.com)
 * @date 2026-05-28
 */
#ifndef _MT_NE_ESTIMATE_H_
#define _MT_NE_ESTIMATE_H_

#include <getopt.h>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "log_factorial.h"
#include "version.h"

class NeEstimator {
public:
    /// Type alias to the shared LogFactorial utility (in src/log_factorial.h).
    using LogFactorial = ::mitoquest::LogFactorial;

    // One transmission pair, after upstream QC.
    struct PairData {
        int m_dp;      // mother total depth
        int m_ad_alt;  // mother ALT-supporting reads
        int c_dp;      // child  total depth
        int c_ad_alt;  // child  ALT-supporting reads
    };

    /// Optional Wonnapinij/Kimura cross-check output.
    struct KimuraCheck {
        bool   computed       = false; // populated only when --cross-check kimura is set
        double b              = 0.0;   // bottleneck parameter, range (0, 1]
        double ne_kimura      = 0.0;   // Ne implied by 1 / (1 - b) for one generation
        size_t n_informative  = 0;     // pairs that contributed to b (per Wonnapinij filter)
        std::string note;              // free-text caveat (e.g. "b clipped to [eps, 1-eps]")

        // ---- Optional non-parametric bootstrap CI on the same pair set ----
        // Populated only when n_bootstrap > 0 in compute_kimura_check().
        bool   ci_computed       = false;
        int    n_bootstrap       = 0;
        uint64_t bootstrap_seed  = 0;
        double b_ci_low          = 0.0;     // 2.5 percentile of bootstrap b
        double b_ci_high         = 0.0;     // 97.5 percentile of bootstrap b
        double ne_kimura_ci_low  = 0.0;     // 2.5 percentile of 1/(1-b)
        double ne_kimura_ci_high = 0.0;     // 97.5 percentile of 1/(1-b)
    };

    // Final estimate.
    struct Result {
        int    ne          = 0;
        int    ci_low      = 0;
        int    ci_high     = 0;
        double max_log_lik = 0.0;
        size_t n_pairs     = 0;
        bool   ci_low_clipped  = false;   // optimum is at the search-boundary
        bool   ci_high_clipped = false;
        KimuraCheck kimura;               // populated only when requested
    };

    // CLI configuration parsed from argv.
    struct Config {
        std::string input_tsv;       // input pairs TSV (from `trans-prep`)
        std::string output_file;     // JSON output (empty => stdout)
        double      min_vaf;         // maternal VAF lower bound (inclusive)
        double      max_vaf;         // maternal VAF upper bound (inclusive)
        int         min_ne;          // smallest Ne to consider (>= 1)
        int         max_ne;          // largest Ne to consider
        int         threads;         // worker threads for log-likelihood
        bool        kimura_check;    // run Wonnapinij cross-check?  (--cross-check kimura)
        int         kimura_bootstrap;// non-parametric bootstrap iterations for the Kimura CI;
                                     // 0 disables the CI computation.
        uint64_t    kimura_seed;     // RNG seed for the Kimura bootstrap
    };

    explicit NeEstimator(int argc, char* argv[]);
    explicit NeEstimator(Config config);                // for tests / library callers
    ~NeEstimator() = default;

    // Run the end-to-end pipeline (parse -> filter -> estimate -> write JSON).
    Result run();

    const Config& config() const { return _config; }

    // -----------------------------------------------------------------
    // Pure helpers, exposed for unit tests.
    // -----------------------------------------------------------------

    // Read the TSV produced by `mitoquest trans-prep`, applying the
    // QC == "PASS" gate and the maternal-VAF window.
    static std::vector<PairData> load_pairs(const std::string& tsv_path,
                                            double min_vaf, double max_vaf);

    // Per-pair log-likelihood.
    static double compute_ll_single(const PairData& pd, int ne,
                                    const LogFactorial& lf);

    // Global log-likelihood (single-threaded).
    static double compute_global_ll(int ne,
                                    const std::vector<PairData>& data,
                                    const LogFactorial& lf);

    // Global log-likelihood with optional thread-pool parallelism.
    // Falls back to sequential when `threads <= 1`.
    static double compute_global_ll_parallel(int ne,
                                             const std::vector<PairData>& data,
                                             const LogFactorial& lf,
                                             int threads);

    // Brute-force integer scan over [min_ne, max_ne]; see header comment
    // for why golden-section search is unsafe here.
    static int find_optimal_ne(const std::vector<PairData>& data,
                               const LogFactorial& lf,
                               int min_ne, int max_ne, int threads = 1);

    // Full estimate (point + 95% CI by profile likelihood, threshold = -1.92).
    static Result estimate(const std::vector<PairData>& data,
                           int min_ne = 1, int max_ne = 200,
                           int threads = 1);

    // Determine the LogFactorial cache size needed by `data`.
    // Cache must cover max(child DP) and the search ceiling for Ne.
    static int required_cache_size(const std::vector<PairData>& data, int max_ne);

    // Stable log-sum-exp for two terms: log(exp(a) + exp(b)).
    static double log_sum_exp_pair(double a, double b);

    // -----------------------------------------------------------------
    // Wonnapinij / Kimura cross-check helpers (optional).
    // -----------------------------------------------------------------

    /**
     * @brief Compute the Wonnapinij bottleneck parameter `b` from a set
     *        of M-C pairs, optionally with a non-parametric bootstrap CI.
     *
     * Uses the maternal point estimate p_m = m_alt / m_dp and child
     * point estimate p_c = c_alt / c_dp, then for each pair:
     *
     *     observed deviation        : d_i  = (p_c - p_m)^2
     *     sampling-noise correction : s_i  = p_m (1 - p_m) / m_dp
     *                                       + p_c (1 - p_c) / c_dp
     *     numerator                 : Sigma_i (d_i - s_i)
     *     denominator               : Sigma_i p_m_i (1 - p_m_i)
     *     normalised drift variance : V    = numerator / denominator
     *
     * Then 1 - b = V (Wonnapinij 2008, eq. 4 with sampling correction
     * from Wonnapinij 2010), so
     *
     *     b = 1 - V,    Ne_kimura = 1 / (1 - b)   for one generation.
     *
     * `b` is clipped to (eps, 1 - eps) to guard against finite-sample
     * outliers; the clipping event is reported in KimuraCheck::note.
     *
     * When `n_bootstrap > 0`, this method additionally performs a
     * non-parametric pair-level bootstrap (sample pairs with replacement,
     * recompute b on each resample) and returns the 2.5 / 97.5 percentile
     * confidence interval for both `b` and `Ne_kimura`.
     */
    static KimuraCheck compute_kimura_check(const std::vector<PairData>& data,
                                            int n_bootstrap = 0,
                                            uint64_t seed   = 42);

private:
    NeEstimator(const NeEstimator&)            = delete;
    NeEstimator& operator=(const NeEstimator&) = delete;

    static void usage();
    void _parse_args(int argc, char* argv[]);
    void _write_json(const Result& r, std::ostream& out) const;

    Config      _config;
    std::string _cmdline_string;
};

#endif // _MT_NE_ESTIMATE_H_
