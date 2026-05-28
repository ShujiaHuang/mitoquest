/**
 * @file ne_estimate.h
 * @brief Estimate the mtDNA bottleneck size (Ne) from mother-child
 *        transmission pairs via maximum likelihood, with an optional
 *        Wonnapinij/Kimura cross-check.
 *
 * ====================================================================
 * Statistical model (primary estimator: continuous Beta-diffusion MLE)
 * ====================================================================
 *
 * For each independent mother-child (M-C) transmission pair the child's
 * heteroplasmy is modelled as a Kimura-diffusion draw from the maternal
 * heteroplasmy after a single-generation bottleneck of Ne transmitting
 * units:
 *
 *   p_child | p_mother  ~  Beta( p_mother * (Ne - 1),
 *                                (1 - p_mother) * (Ne - 1) )
 *
 *   c_alt | p_child     ~  Binomial(c_dp, p_child)
 *
 * After marginalising out p_child:
 *
 *   c_alt | p_mother    ~  BetaBinomial(c_dp, p_mother * (Ne - 1),
 *                                             (1 - p_mother) * (Ne - 1))
 *
 * The maternal p_mother is taken as the point estimate m_alt / m_dp
 * (at typical mtDNA depths >= 100 this is virtually identical to the
 * fully marginalised integral over a Beta posterior).
 *
 * This continuous model is appropriate for mtDNA transmission because
 * the child's actual heteroplasmy is shaped not only by the initial
 * bottleneck sampling but also by post-bottleneck vegetative
 * segregation during cell division.  The Kimura diffusion captures
 * *all* these sources of variance as a single effective Ne.
 *
 * The previous discrete model (k ~ BetaBin, c ~ Bin(c_dp, k/Ne))
 * restricted child heteroplasmy to the coarse grid {0, 1/Ne, ..., 1}
 * and suffered from systematic upward bias at high sequencing depths:
 * the tight Binomial likelihood forced the MLE to inflate Ne to obtain
 * a fine enough grid, giving answers ~5-10x above the deCODE / Kimura
 * consensus.  See release_v1.8.2.md for details.
 *
 * The discrete model is still available via `--model discrete` for
 * specialised use cases (e.g. virus-passage bottleneck experiments
 * where the physical inoculum count is the target).
 *
 * 95% confidence interval is the *contiguous* Ne range with
 *   logL(Ne) >= logL_max - chi2_{1, 0.95} / 2  =  logL_max - 1.92
 * (Wilks's theorem on the profile log-likelihood; see Cox & Hinkley
 *  1974, "Theoretical Statistics", Section 9.3).
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
 * On well-behaved data the continuous MLE and Wonnapinij b should
 * agree closely.  Residual discrepancies indicate either heavy-tailed
 * outliers (use `--kimura-trim`) or model misspecification.
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

        // ---- Robust trimmed Kimura (drops top trim_frac of high-drift pairs) ----
        // Populated only when trim_frac > 0 in compute_kimura_check().
        // Standard Wonnapinij b is variance-of-moments and is *not* robust to
        // outliers (NUMTs / sequencing errors / mixed populations).  A handful
        // of high-drift pairs can collapse Ne_kimura by an order of magnitude
        // even when the bulk of pairs are well-behaved.  The trimmed estimator
        // ranks pairs by their per-pair contribution F_i = (d_i - s_i) / w_i
        // and drops the top `trim_frac` of pairs (highest-drift) before
        // computing b.  When the gap between trimmed and untrimmed is large,
        // the data contain outliers that the standard Kimura is over-fitting.
        bool   trimmed_computed     = false;
        double trim_frac            = 0.0;  // user-specified --kimura-trim
        size_t n_after_trim         = 0;    // informative pairs surviving the trim
        double b_trimmed            = 0.0;
        double ne_kimura_trimmed    = 0.0;

        // ---- Per-pair drift outlier diagnostic (top-K) -----------------
        // Populated only when top_drift_k > 0 in compute_kimura_check().
        // Each entry identifies a pair by its 0-based index in the input
        // PairData vector and reports the per-pair Wonnapinij contribution
        // F_i = (d_i - s_i) / w_i, sorted in descending order.
        struct DriftOutlier {
            size_t pair_index;   // 0-based index into the input PairData vector
            int    m_dp;
            int    m_ad_alt;
            int    c_dp;
            int    c_ad_alt;
            double m_vaf;
            double c_vaf;
            double f_i;          // (d_i - s_i) / w_i, the Wonnapinij per-pair F
        };
        std::vector<DriftOutlier> top_drift_outliers;
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
        std::string model;           // "continuous" (default) or "discrete"
        double      min_vaf;         // maternal VAF lower bound (inclusive)
        double      max_vaf;         // maternal VAF upper bound (inclusive)
        int         min_ne;          // smallest Ne to consider (>= 1)
        int         max_ne;          // largest Ne to consider
        int         threads;         // worker threads for log-likelihood
        bool        kimura_check;    // run Wonnapinij cross-check?  (--cross-check kimura)
        int         kimura_bootstrap;// non-parametric bootstrap iterations for the Kimura CI;
                                     // 0 disables the CI computation.
        uint64_t    kimura_seed;     // RNG seed for the Kimura bootstrap
        double      kimura_trim;     // fraction of high-drift pairs to drop from
                                     // the trimmed Kimura cross-check;
                                     // 0.0 disables (default).
        int         top_drift_k;     // emit the top-K drift outlier pairs in the
                                     // JSON output for diagnostic inspection;
                                     // 0 disables (default).
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

    // Per-pair log-likelihood (discrete bottleneck model).
    static double compute_ll_single(const PairData& pd, int ne,
                                    const LogFactorial& lf);

    // Per-pair log-likelihood (continuous Beta-diffusion model).
    // Models: p_child | p_m ~ Beta(p_m*(Ne-1), (1-p_m)*(Ne-1)),
    // then c_alt | p_child ~ Bin(c_dp, p_child), marginalized to:
    //   c_alt ~ BetaBin(c_dp, p_m*(Ne-1), (1-p_m)*(Ne-1)).
    static double compute_ll_single_continuous(const PairData& pd, int ne,
                                               const LogFactorial& lf);

    // Global log-likelihood (single-threaded).
    static double compute_global_ll(int ne,
                                    const std::vector<PairData>& data,
                                    const LogFactorial& lf,
                                    bool continuous = false);

    // Global log-likelihood with optional thread-pool parallelism.
    // Falls back to sequential when `threads <= 1`.
    static double compute_global_ll_parallel(int ne,
                                             const std::vector<PairData>& data,
                                             const LogFactorial& lf,
                                             int threads,
                                             bool continuous = false);

    // Brute-force integer scan over [min_ne, max_ne]; see header comment
    // for why golden-section search is unsafe here.
    static int find_optimal_ne(const std::vector<PairData>& data,
                               const LogFactorial& lf,
                               int min_ne, int max_ne, int threads = 1,
                               bool continuous = false);

    // Full estimate (point + 95% CI by profile likelihood, threshold = -1.92).
    static Result estimate(const std::vector<PairData>& data,
                           int min_ne = 1, int max_ne = 200,
                           int threads = 1, bool continuous = false);

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
                                            int      n_bootstrap = 0,
                                            uint64_t seed        = 42,
                                            double   trim_frac   = 0.0,
                                            int      top_drift_k = 0);

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
