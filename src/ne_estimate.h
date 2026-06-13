/**
 * @file ne_estimate.h
 * @brief Estimate the mtDNA bottleneck size (Ne) from mother-child
 *        transmission pairs via Maximum Marginal Likelihood Estimation
 *        (MMLE), with an optional Wonnapinij/Kimura cross-check.
 *
 * ====================================================================
 * Statistical model (primary estimator: continuous Beta-diffusion MMLE)
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
 * After analytically marginalising out the latent p_child:
 *
 *   c_alt | p_mother    ~  BetaBinomial(c_dp, p_mother * (Ne - 1),
 *                                             (1 - p_mother) * (Ne - 1))
 *
 * yielding the per-pair *marginal* likelihood that depends only on Ne
 * and the observed counts.  Because the M-C pairs are treated as
 * independent, the global objective is the sum of per-pair marginal
 * log-likelihoods (a composite / pseudo-likelihood that is consistent
 * for Ne under the standard Wright-Fisher / Kimura assumptions).
 * Maximising this composite marginal likelihood over Ne yields the
 * Maximum Marginal Likelihood Estimator (MMLE) reported by this tool.
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
 * the tight Binomial likelihood forced the MMLE to inflate Ne to obtain
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
 * On well-behaved data the continuous MMLE and Wonnapinij b should
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

    // One transmission pair, after upstream QC.  For two-generation rows
    // (produced without a GM FAM) has_g == 0 and the g_* fields are zero.
    // For three-generation G-M-C trios (HAS_G = 1 in the TSV) g_dp /
    // g_ad_alt carry the grandmother's counts and the row-level
    // log-likelihood switches to the trio marginal formula in the
    // continuous model (see compute_ll_trio_continuous()).
    struct PairData {
        int m_dp;      // mother total depth
        int m_ad_alt;  // mother ALT-supporting reads
        int c_dp;      // child  total depth
        int c_ad_alt;  // child  ALT-supporting reads
        int g_dp     = 0;   // grandmother total depth (HAS_G=1 rows only)
        int g_ad_alt = 0;   // grandmother ALT reads  (HAS_G=1 rows only)
        int has_g    = 0;   // 1 iff this row is a G-M-C trio; 0 if MC pair
        // Family identifiers (read from TSV FAM_ID/MOTHER_ID/CHILD_ID columns;
        // used by per-family mode; empty strings for legacy TSVs without these).
        std::string fam_id;
        std::string mother_id;
        std::string child_id;
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

    // One family's worth of transmission data (grouped by FAM_ID + MOTHER_ID).
    struct FamilyData {
        std::string fam_id;
        std::string mother_id;
        std::vector<std::string> child_ids;  // may have >1 children
        std::vector<PairData>    pairs;      // all variant sites for this family
    };

    // Per-family estimation result.
    struct FamilyResult {
        std::string fam_id;
        std::string mother_id;
        size_t      n_children    = 0;
        size_t      n_pairs       = 0;      // total variant sites used
        size_t      n_informative = 0;      // sites with 0 < p_M < 1
        double      ne            = 0.0;    // per-family MMLE Ne
        double      ci_low        = 0.0;
        double      ci_high       = 0.0;
        double      max_log_lik   = 0.0;
        bool        ci_low_clipped  = false;
        bool        ci_high_clipped = false;
        double      mean_mother_dp  = 0.0;  // mean depth for quality assessment
        double      mean_child_dp   = 0.0;
        bool        skipped       = false;  // true when n_informative < min_family_sites
        KimuraCheck kimura;                 // per-family Kimura cross-check
        std::string warning;                // e.g. "small sample" or "Ne near boundary"
    };

    // Final estimate.
    struct Result {
        double ne          = 0.0;
        double ci_low      = 0.0;
        double ci_high     = 0.0;
        double max_log_lik = 0.0;
        size_t n_pairs     = 0;
        bool   ci_low_clipped  = false;   // optimum is at the search-boundary
        bool   ci_high_clipped = false;
        KimuraCheck kimura;               // populated only when requested
        // Per-family results (populated only when --per-family is set).
        std::vector<FamilyResult> family_results;
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
        std::string bin_simulation_file; // when non-empty, write a per-bin
                                         // observed-drift summary TSV (per
                                         // maternal-VAF bin: observed mean
                                         // drift vs analytical Kimura
                                         // prediction p(1-p)/Ne at the
                                         // fitted Ne and its 95% CI).
        int         bin_simulation_n_bins; // number of equal-width bins to use
                                           // when --bin-simulation is set.
        std::string ne_profile_file;     // when non-empty, write an Ne-profile
                                         // TSV that scores every candidate Ne
                                         // under both the MMLE and Kimura
                                         // models (dual-objective Ne scan).
        double      ne_profile_step;     // grid step on the Ne axis for
                                         // --ne-profile (default 0.1).
        // NEW: per-family options
        bool        per_family         = false;   // --per-family
        int         min_family_sites   = 3;       // --min-family-sites
        std::string per_family_output_file;       // --per-family-output FILE (TSV)
    };

    /// One row of the per-bin observed-vs-theoretical drift summary.
    /// All quantities are computed on the same pair set the MMLE /
    /// Kimura cross-check was fit on.
    struct BinSimulationRow {
        int    bin_idx       = 0;
        double bin_low       = 0.0;   // maternal-VAF bin lower edge (inclusive)
        double bin_high      = 0.0;   // maternal-VAF bin upper edge (exclusive,
                                      // except the top bin which is inclusive)
        double bin_center    = 0.0;   // (bin_low + bin_high) / 2
        size_t n_pairs       = 0;     // pairs with p_m in this bin AND informative
        double mean_pm       = 0.0;   // empirical mean p_mother in bin
        double mean_pc       = 0.0;   // empirical mean p_child  in bin
        double obs_var       = 0.0;   // empirical mean of (p_c - p_m)^2 (raw drift)
        double obs_var_corr  = 0.0;   // empirical mean of (d_i - s_i)
                                      //   = sampling-corrected drift squared
        double obs_F         = 0.0;   // empirical mean of F_i = (d_i - s_i) / w_i
                                      //   (= bin estimate of 1 - b)
        double obs_F_se      = 0.0;   // standard error of mean F_i in bin
    };

    /// One row of the Ne-profile scoring scan.  For each candidate Ne we
    /// report two independent goodness-of-fit metrics so the user can see
    /// which Ne each of the two estimators in the program prefers:
    ///
    ///   * `mmle_log_lik` -- global marginal log-likelihood under the configured
    ///                       model (continuous Beta-diffusion or discrete
    ///                       Beta-Binomial).  Maximised at the fitted
    ///                       Ne_MMLE.
    ///   * `kimura_ssr`   -- sum-of-squared per-pair residuals under the
    ///                       one-generation Wright-Fisher prediction
    ///                           E[d_i - s_i]  =  p_m_i (1 - p_m_i) / Ne.
    ///                       Minimised at the analytic Kimura-SSR best
    ///                       Ne =  Sigma w_i^2  /  Sigma r_i w_i.
    ///
    /// Both metrics are normalised in post-processing so the plotting
    /// script can render comparable curves on linear or log axes.
    struct NeProfileRow {
        double ne_candidate    = 0.0;
        double mmle_log_lik    = 0.0;
        double mmle_delta_2ll  = 0.0;   // -2 (LL - LL_max); 0 at fitted Ne_MMLE
        double kimura_ssr      = 0.0;   // Sigma_i (r_i - w_i / Ne)^2
        double kimura_norm_ssr = 0.0;   // ssr / ssr_min (>= 1, =1 at best fit)
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

    // Per-pair log-likelihood (continuous Beta-diffusion model, integer Ne).
    // Models: p_child | p_m ~ Beta(p_m*(Ne-1), (1-p_m)*(Ne-1)),
    // then c_alt | p_child ~ Bin(c_dp, p_child), marginalized to:
    //   c_alt ~ BetaBin(c_dp, p_m*(Ne-1), (1-p_m)*(Ne-1)).
    static double compute_ll_single_continuous(const PairData& pd, int ne,
                                               const LogFactorial& lf);

    // Per-pair log-likelihood (continuous model, real-valued Ne).
    static double compute_ll_single_continuous(const PairData& pd, double ne,
                                               const LogFactorial& lf);

    // Per-row marginal log-likelihood for a three-generation G-M-C trio
    // under the continuous Beta-diffusion model.  The mother's latent
    // heteroplasmy p_M is integrated out analytically:
    //
    //   I(Ne) = int_0^1 Beta(p_M | alpha_G, beta_G)
    //                  * Bin(k_M | d_M, p_M)
    //                  * BetaBin(k_C | d_C, p_M (Ne-1), (1-p_M)(Ne-1))
    //           dp_M
    //
    // where alpha_G = p_hat_G * (Ne - 1) and beta_G = (1 - p_hat_G) * (Ne - 1)
    // with p_hat_G = g_ad_alt / g_dp.  The closed form is:
    //
    //   I(Ne) = C(d_M, k_M) C(d_C, k_C)
    //           * B(alpha_G + k_M + k_C,
    //               beta_G  + (d_M - k_M) + (d_C - k_C))
    //           / B(alpha_G, beta_G)
    //
    // When pd.has_g == 0 (i.e. the row is a two-generation MC pair), this
    // function falls back to the standard two-generation continuous
    // likelihood (compute_ll_single_continuous).  When the grandmother is
    // homoplasmic (p_hat_G in {0, 1}) the row is non-informative for Ne
    // and the function returns 0.0 (a Ne-independent constant).
    static double compute_ll_trio_continuous(const PairData& pd, double ne,
                                             const LogFactorial& lf);

    // Gauss-Legendre quadrature fallback for the trio marginal likelihood.
    // Used by the unit tests to validate the closed-form formula; not the
    // default path in the optimiser.
    static double compute_ll_trio_quadrature(const PairData& pd, double ne,
                                             const LogFactorial& lf,
                                             int n_nodes = 64);

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

    // Global log-likelihood for the continuous model with real-valued Ne.
    static double compute_global_ll_continuous(double ne,
                                              const std::vector<PairData>& data,
                                              const LogFactorial& lf,
                                              int threads = 1);

    // Real-valued Ne optimizer for the continuous model.
    // Phase 1: coarse integer scan; Phase 2: golden-section refinement.
    static double find_optimal_ne_continuous(const std::vector<PairData>& data,
                                            const LogFactorial& lf,
                                            int min_ne, int max_ne,
                                            int threads = 1);

    // Full estimate with real-valued Ne for the continuous model.
    static Result estimate_continuous(const std::vector<PairData>& data,
                                     int min_ne = 1, int max_ne = 200,
                                     int threads = 1);

    // Full estimate (point + 95% CI by profile likelihood, threshold = -1.92).
    // For discrete model only (integer Ne scan).
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

    /**
     * @brief Aggregate per-pair drift statistics into equal-width
     *        maternal-VAF bins for the per-bin drift summary plot.
     *
     * For each bin the row reports:
     *   * `obs_var`      = mean (p_c - p_m)^2          (raw observed drift)
     *   * `obs_var_corr` = mean (d_i - s_i)            (sampling-corrected drift)
     *   * `obs_F`        = mean (d_i - s_i) / [p_m(1-p_m)]
     *                                                  (bin estimate of 1-b)
     *   * `obs_F_se`     = SE of mean F_i in the bin
     *
     * The downstream theoretical curves p_m(1-p_m) / Ne and 1/Ne are
     * computed by the caller (or the plotting script) from the fitted Ne.
     *
     * Pairs with p_m outside [vaf_low, vaf_high] or with p_m equal to 0
     * or 1 (mother homoplasmic) are skipped.  `n_bins` is clamped to >= 1.
     */
    static std::vector<BinSimulationRow>
    compute_bin_simulation(const std::vector<PairData>& data,
                           double vaf_low, double vaf_high, int n_bins);

    /**
     * @brief Score every candidate Ne in [min_ne, max_ne] (step `step`)
     *        under both the MMLE marginal log-likelihood and the Kimura
     *        sum-of-squared-residuals metric.
     *
     * The MMLE column uses the continuous Beta-diffusion marginal log-
     * likelihood when `continuous == true` and the discrete Beta-Binomial
     * marginal log-likelihood otherwise.  The Kimura column is independent
     * of the model selection: for every informative pair we compute
     *
     *     residual_i(Ne)  =  (d_i - s_i)  -  p_m_i (1 - p_m_i) / Ne
     *     ssr(Ne)         =  Sigma_i residual_i(Ne)^2
     *
     * After the scan we normalise: `mmle_delta_2ll = -2 (LL - LL_max)` and
     * `kimura_norm_ssr = ssr / ssr_min`.  The grid is clamped so that
     * `step > 0` and `min_ne >= 1`.
     */
    static std::vector<NeProfileRow>
    compute_ne_profile(const std::vector<PairData>& data,
                       const LogFactorial& lf,
                       double min_ne, double max_ne, double step,
                       int threads, bool continuous);

    // -----------------------------------------------------------------
    // Per-family estimation helpers.
    // -----------------------------------------------------------------

    // Group loaded pairs by (fam_id, mother_id) into FamilyData structs.
    // Pairs with empty fam_id (legacy TSVs) are grouped into a single
    // family with fam_id = "ALL".
    static std::vector<FamilyData>
    group_into_families(const std::vector<PairData>& data);

    // Estimate Ne for a single family using the continuous MMLE.
    // Returns a FamilyResult with skipped=true when n_informative < min_family_sites.
    static FamilyResult estimate_family(const FamilyData& fam,
                                         int min_ne, int max_ne,
                                         int min_family_sites);

    // Estimate Ne for all families (embarrassingly parallel across families).
    static std::vector<FamilyResult>
    estimate_all_families(const std::vector<FamilyData>& families,
                          int min_ne, int max_ne,
                          int min_family_sites,
                          int threads);

    // Per-family Kimura cross-check (site-level bootstrap, not pair-level).
    static KimuraCheck compute_family_kimura_check(
        const FamilyData& fam,
        int      n_bootstrap = 0,
        uint64_t seed        = 42,
        double   trim_frac   = 0.0);

    // Write per-family results as a TSV file (one row per family).
    void _write_family_tsv(const std::vector<FamilyResult>& results,
                           std::ostream& out) const;

    /// Closed-form Kimura-SSR best-fit Ne under the one-generation
    /// Wright-Fisher prediction:  Ne_best  =  Sigma w_i^2 / Sigma r_i w_i.
    /// Returns NaN when the denominator is non-positive (corrected drift
    /// summed across pairs goes the "wrong" direction; i.e. sampling
    /// correction overshoots and the data are statistically
    /// indistinguishable from no drift).
    static double kimura_ssr_best_ne(const std::vector<PairData>& data);

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
