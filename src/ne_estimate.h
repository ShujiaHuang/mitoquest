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
 * heteroplasmy is modelled with a working Beta transition approximation
 * matched to the mean and one-generation variance of neutral drift:
 *
 *   p_child | p_mother  ~  Beta(p_mother * (Ne - 1), (1 - p_mother) * (Ne - 1))
 *
 *   c_alt | p_child     ~  Binomial(c_dp, p_child)
 *
 * After analytically marginalising out the latent p_child:
 *
 *   c_alt | p_mother    ~  BetaBinomial(c_dp, p_mother * (Ne - 1), (1 - p_mother) * (Ne - 1))
 *
 * yielding the per-pair *marginal* likelihood that depends only on Ne
 * and the observed counts.  Because the M-C pairs are treated as
 * independent, the global objective is the sum of per-pair marginal
 * log-likelihoods (a composite / pseudo-likelihood that is consistent
 * for Ne under the standard Wright-Fisher / Kimura assumptions).
 * Maximising this composite marginal likelihood over Ne yields the
 * Maximum Marginal Likelihood Estimator (MMLE) reported by this tool.
 *
 * The maternal frequency p_mother is integrated out by DEFAULT
 * (ModelOptions::marginalize_maternal).  The CLI exposes only the opt-out
 * `--no-maternal-marginalization`, which restores the legacy plug-in -- there
 * is no positive form, because asking for the default would be a no-op that
 * only looks like a guarantee.  The plug-in
 * alternative takes p_mother = m_alt / m_dp, so the mother's own binomial read
 * noise is NOT modelled and is misread as drift.  The plug-in therefore
 * estimates a depth-dependent pseudo-true value rather than Ne:
 *
 *   1 / Ne_app  =  1 / Ne + 1 / d_M        <=>   Ne_app = Ne d_M / (Ne + d_M)
 *
 * which moves with sequencing depth, so plug-in estimates are not comparable
 * across cohorts or families of different coverage.  Replicated Monte-Carlo
 * (R=8 cohorts of ~1550 gated pairs each, continuous generator) confirms this
 * quantitatively: the plug-in mean lands on Ne d_M/(Ne + d_M) to within 3%,
 * with bias -2.7% at (Ne=5, d=100), -36.4% at (Ne=30, d=50), -21.3% at
 * (Ne=30, d=100) and -48.4% at (Ne=100, d=100).
 *
 * The default replaces the plug-in by an exact integration of p_M over its
 * Beta(k_M+1, d_M-k_M+1) posterior, which is the same maternal prior the
 * legacy discrete model uses.  On the same replicates the residual bias drops
 * to +0.7% / -2.7% / +1.1% / +3.9% and the RMSE improves 1.3x-4.7x, matching
 * the Wonnapinij/Kimura cross-check to within 0.5% in all four
 * configurations.  It does NOT remove model-form error: under a truly integer
 * Wright-Fisher bottleneck at Ne=5 both the plug-in and the marginalised
 * continuous estimator stay ~-50% biased (Kimura: -1.9%), because a continuous
 * Beta transition cannot represent an integer bottleneck -- use
 * `--model discrete` there.
 *
 * Cross-depth comparability, measured (true Ne=30, 20 families x 12 sites =
 * 240 pairs per cohort, d_C=2000, 20 replicate cohorts per depth, d_M swept
 * over 30..5000).  Sweeping d_M moves the plug-in fit over [16.5, 30.3], a
 * spread of 46% of the true Ne and closely tracking -d_M/(Ne+d_M); the
 * marginalised fit moves over [29.4, 34.2], a spread of 16%, i.e. 2.9x
 * tighter.  Excluding the shallowest depth the marginalised fit spans only
 * [29.4, 30.7] over d_M in [60, 5000] -- 4.5% -- so for d_M >= 60 the reported
 * Ne is genuinely comparable across coverage.  The residual is the Beta(1,1)
 * maternal prior being over-dispersed relative to a cohort whose true maternal
 * frequencies are truncated away from 0 and 1: +14% at d_M=30, <= +5% for
 * d_M >= 60, and within ~1 standard error (<= 2.5%) for d_M >= 500.  Use
 * --min-depth to gate very shallow maternal rows out rather than interpreting
 * them at face value.
 *
 * The maternal VAF window is NOT the source of that residual: whether a row
 * passes `k_M/d_M in [min_vaf, max_vaf]` depends on p_M and the prior but not
 * on Ne, so the window contributes an Ne-independent factor to the likelihood
 * and cannot shift the MLE.  Measured with p_M drawn near the window edge
 * (U[0.10,0.30], only 90% of rows passing at d_M=30) the marginalised bias is
 * +13.0% with the default window against +11.6% with `--min-vaf 0
 * --max-vaf 1`, a 1.4 point difference inside the 3.2 point standard error;
 * every depth agrees to < 1.5 points.
 *
 * Trio rows (has_g == 1) already marginalise p_M over the grandmother-informed
 * posterior (compute_ll_trio_gauss_jacobi), so the switch is a no-op for them;
 * the residual plug-in there is p_hat_G, i.e. the bias moves up a generation.
 *
 * Cost: marginalising replaces one Beta-Binomial evaluation per row by an
 * adaptive Gauss-Jacobi rule.  The rule construction (a Golub-Welsch
 * tridiagonal eigendecomposition) is 97.5% of that cost at d_M=d_C=2000, and
 * its key (k_M+1, d_M-k_M+1, order) carries no Ne, so estimate_continuous()
 * builds it once per row and reuses it across the whole Ne scan.  Measured on
 * the 236-pair shipped cohort at --max-ne 200, single-threaded: 19 s without
 * that reuse, 1.69 s with it, and both produce byte-identical JSON.
 *
 * This continuous model treats Ne as an apparent, variance-matching
 * effective bottleneck. It may absorb bottleneck sampling, vegetative
 * segregation, and unmodelled biological or technical variability; it is
 * not a literal molecule count or the exact transient density of a
 * Wright-Fisher diffusion.
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

class NeEstimator {
public:
    /// Type alias to the shared LogFactorial utility (in src/log_factorial.h).
    using LogFactorial = ::mitoquest::LogFactorial;

    /// Model-form switches threaded through the likelihood / optimizer stack.
    /// Passed by value with a trailing default so existing positional callers
    /// (including the unit tests) keep compiling unchanged.
    ///
    /// `marginalize_maternal` is TRUE by default and affects ONLY the
    /// continuous model's two-generation rows (has_g == 0):
    ///   * trio rows (has_g == 1) already integrate p_M over the
    ///     grandmother-informed Beta posterior -> no-op;
    ///   * the discrete model (`--model discrete`) already integrates p_M over
    ///     Beta(m_alt+1, m_dp-m_alt+1) in compute_ll_single() -> no-op, and
    ///     the CLI warns only when the user asks for the switch explicitly
    ///     alongside `--model discrete` (a default-on switch must not make
    ///     every discrete run look misconfigured).
    /// Set it to false for the legacy plug-in p_M = m_alt / m_dp.
    ///
    /// The constructors are user-provided on purpose: a default member
    /// initializer here would be odr-used by the `ModelOptions opts = {}`
    /// default arguments below, i.e. inside the still-incomplete enclosing
    /// class, which clang rejects ("default member initializer ... needed
    /// within definition of enclosing class").  Initialising in the ctor body
    /// keeps this header portable across gcc and clang.
    struct ModelOptions {
        bool marginalize_maternal;

        ModelOptions() : marginalize_maternal(true) {}
        explicit ModelOptions(bool marginalize_maternal_)
            : marginalize_maternal(marginalize_maternal_) {}
    };

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
        // Stable CHROM/POS/ALT key when the input TSV supplies all three;
        // used to resample all children of one locus together in a family.
        std::string locus_id;
    };

    /// Loading diagnostics collected by load_pairs.  Trio rows whose
    /// grandmother is homoplasmic are absorbing G->M->C states and are
    /// dropped at load time: compatible all-homoplasmic descendants
    /// contribute a Ne-independent constant (log 1 = 0), while
    /// incompatible segregating descendants are impossible (-inf) under
    /// the model and would abort the whole fit if passed through.
    struct LoadStats {
        size_t trio_founder_hom_skipped      = 0; // G hom, descendants all-hom
        size_t trio_founder_mismatch_skipped = 0; // G hom, descendants segregating
        size_t min_depth_skipped             = 0; // m_dp or c_dp below --min-depth
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
        double      mean_mother_dp  = 0.0;  // arithmetic mean depth, quality assessment
        double      mean_child_dp   = 0.0;
        // Harmonic means, 1 / mean(1/d).  Depth enters the plug-in bias as
        // 1/d_M, so this is the scale on which limited coverage actually
        // aggregates; the arithmetic mean is dominated by the deepest rows,
        // which are the least biased ones.  Neither is a back-correction
        // handle -- see the note on Result::mean_mother_dp.
        double      harmonic_mother_dp = 0.0;
        double      harmonic_child_dp  = 0.0;
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

        // Depth descriptors of the fitted pair set, computed on the SAME rows
        // that entered the likelihood (i.e. after the QC / VAF / min-depth
        // gates).  Reported for transparency only.
        //
        // They must NOT be advertised as a plug-in back-correction handle:
        // 1/Ne_app = 1/Ne + 1/d_M is a *single-depth* identity, and on a
        // depth-heterogeneous cohort the plug-in estimand is a per-row mixture
        // that no scalar mean of d_M reproduces.  Measured on a cohort with
        // d_M in {30,60,100,300,1000,3000} and true Ne=30, the plug-in fit
        // (~25.9) is predicted as 28.9 by the arithmetic mean (753), 22.7 by
        // the harmonic mean (93) and 22.7 by the w-weighted harmonic mean (93);
        // the implied effective depth (~190) lies between the two.  Both scales
        // are therefore emitted and neither is a correction recipe.
        double mean_mother_dp      = 0.0;  // arithmetic mean of m_dp
        double mean_child_dp       = 0.0;  // arithmetic mean of c_dp
        double harmonic_mother_dp  = 0.0;  // 1 / mean(1 / m_dp); the scale on
                                           // which 1/d_M effects aggregate
        double harmonic_child_dp   = 0.0;

        // Echoes of the model-form switches actually applied to this fit, so a
        // JSON reader can tell which likelihood the reported LL belongs to.
        // `maternal_marginalization` is true for a default continuous-model
        // run and stays false when the plug-in was selected
        // (--no-maternal-marginalization) or when the selected model ignores
        // the switch (--model discrete already integrates p_M; the CLI reports
        // that on stderr when the user asked for it explicitly).
        //
        // Max_Marginal_LogLik is NOT comparable across this flag.  Switching it
        // on replaces the plug-in child term log BB(k_C; d_C, p_hat_M s,
        // (1-p_hat_M) s) by log E_{Beta(k_M+1, d_M-k_M+1)}[BB(k_C; d_C, p_M s,
        // (1-p_M) s)], which is not the old term plus a constant: the fitted Ne
        // moves as well.  Only the joint read-sampling normalisation
        // log C(d_M,k_M) + log B(k_M+1, d_M-k_M+1) = -log(d_M+1) is
        // Ne-independent, and hence MLE-neutral.
        bool   maternal_marginalization = false;
        int    min_depth                = 0;

        KimuraCheck kimura;               // populated only when requested
        // Per-family results (populated only when --per-family is set).
        std::vector<FamilyResult> family_results;

        // Loading diagnostics (trio founder-homoplasmy filters etc.).
        LoadStats load_stats;
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
        // NEW: maternal-noise handling.  Marginalisation is ON by default; it
        // is the only continuous-model treatment of maternal read noise that
        // yields coverage-comparable Ne (see the class comment for the
        // measured cross-depth spread).  --no-maternal-marginalization selects
        // the legacy plug-in p_M = m_alt / m_dp.
        bool        maternal_marginalization = true;  // cleared only by
                                                      // --no-maternal-marginalization;
                                                      // see ModelOptions
        int         min_depth          = 0;       // --min-depth; drop rows with
                                                  // m_dp < N or c_dp < N (0 = off).
                                                  // Complements marginalization
                                                  // rather than replacing it: at
                                                  // true Ne=30 on a
                                                  // heterogeneous-depth cohort,
                                                  // --min-depth 500 reaches -1.3%
                                                  // plug-in bias only after
                                                  // discarding 2/3 of the pairs
                                                  // (1569 -> 528), while
                                                  // marginalization reaches +2.1%
                                                  // keeping all 1569.  Its role
                                                  // under the new default is to
                                                  // gate out the very shallow
                                                  // maternal rows (d_M ~ 30)
                                                  // where the residual Beta(1,1)
                                                  // prior mismatch is still
                                                  // ~+14%.
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
    // QC == "PASS" gate and the maternal-VAF window.  Trio rows with a
    // homoplasmic grandmother are dropped here (see LoadStats); pass
    // `stats` to receive the skip counters for logging / JSON output.
    // `min_depth > 0` additionally drops rows with m_dp < min_depth or
    // c_dp < min_depth (applied to BOTH sides: a shallow child gives an
    // uninformative k_C just as a shallow mother gives an unreliable
    // p_hat_M).  The VAF window is evaluated on m_ad_alt / m_dp, so
    // min_depth interacts with it only through which rows survive.
    static std::vector<PairData> load_pairs(const std::string& tsv_path,
                                            double min_vaf, double max_vaf,
                                            LoadStats* stats = nullptr,
                                            int min_depth = 0);

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

    // Per-row marginal log-likelihood for a TWO-generation M-C pair under the
    // continuous model with the maternal heteroplasmy integrated out instead of
    // plugged in.  With s = Ne - 1:
    //
    //   logL = log C(d_M, k_M) + log B(k_M+1, d_M-k_M+1)
    //          + log E_{p_M ~ Beta(k_M+1, d_M-k_M+1)}
    //                  [ BetaBin(k_C | d_C, s p_M, s (1 - p_M)) ]
    //
    // The Beta(k_M+1, d_M-k_M+1) weight is the exact posterior of p_M under a
    // Beta(1,1) prior -- the SAME maternal prior compute_ll_single() uses for
    // the discrete model, so this option aligns the continuous model's maternal
    // treatment with the discrete one rather than introducing a new prior.
    //
    // Two exact identities pin the normalisation and tie this to shipped code:
    //   * log C(d_M,k_M) + log B(k_M+1, d_M-k_M+1) = -log(d_M+1), a constant
    //     independent of Ne, so the MLE is unaffected by keeping it;
    //   * at Ne=3 with p_hat_G=0.5 the trio prior Beta(p_hat_G s, (1-p_hat_G) s)
    //     is exactly Beta(1,1), hence compute_ll_trio_continuous() evaluated on
    //     a row with g_dp=2, g_ad_alt=1 and Ne=3 returns precisely this value
    //     (verified against the production binary to all printed decimals).
    //
    // The child Beta-Binomial factor is a degree-d_C polynomial in p_M, so the
    // expectation is computed with the same adaptive Gauss-Jacobi driver as the
    // trio path.  Its rule parameters (k_M+1, d_M-k_M+1) do not depend on Ne,
    // so a single TrioQuadratureCache serves the whole Ne scan.
    //
    // Ne == 1 is the degenerate complete-drift limit and is handled in closed
    // form: E[p_M] = (k_M+1)/(d_M+2), giving log((1-E[p_M])/(d_M+1)) for k_C=0,
    // log(E[p_M]/(d_M+1)) for k_C=d_C and -inf otherwise -- the true Ne->1+
    // limit of the integral (endpoint rows converge to it, interior rows
    // diverge to -inf from both sides, so the composite LL is continuous).
    static double compute_ll_single_marginalized(const PairData& pd, double ne,
                                                 const LogFactorial& lf);

    // Per-row marginal log-likelihood for a three-generation G-M-C trio
    // under the continuous Beta-diffusion model.  The mother's latent
    // heteroplasmy p_M is integrated with adaptive Gauss-Jacobi quadrature:
    //
    //   I(Ne) = int_0^1 Beta(p_M | alpha_G, beta_G)
    //                  * Bin(k_M | d_M, p_M)
    //                  * BetaBin(k_C | d_C, p_M (Ne-1), (1-p_M)(Ne-1))
    //           dp_M
    //
    // where alpha_G = p_hat_G * (Ne - 1) and beta_G = (1 - p_hat_G) * (Ne - 1)
    // with p_hat_G = g_ad_alt / g_dp.  The child Beta-Binomial factor is not
    // conjugate in p_M, so this integral cannot be reduced to one Beta ratio.
    // The Gauss-Jacobi rule absorbs the maternal Beta posterior weight and is
    // refined until stable, with a polynomial-exact order cap.
    //
    // When pd.has_g == 0 (i.e. the row is a two-generation MC pair), this
    // function falls back to the standard two-generation continuous
    // likelihood (compute_ll_single_continuous).  At Ne == 1 and for a
    // homoplasmic grandmother, the exact degenerate G-M-C probabilities are
    // evaluated instead of discarding the row.
    static double compute_ll_trio_continuous(const PairData& pd, double ne,
                                             const LogFactorial& lf);

    // Diagnostic Gauss-Legendre evaluation of the trio marginal likelihood.
    // It evaluates the correct integrand but does not handle endpoint-singular
    // Beta weights as robustly as the production Gauss-Jacobi path.
    static double compute_ll_trio_quadrature(const PairData& pd, double ne,
                                             const LogFactorial& lf,
                                             int n_nodes = 64);

    // Global log-likelihood (single-threaded). When `continuous` is true,
    // HAS_G rows dispatch to the G-M-C trio likelihood.
    static double compute_global_ll(int ne,
                                    const std::vector<PairData>& data,
                                    const LogFactorial& lf,
                                    bool continuous = false);

    // Global log-likelihood with optional thread-pool parallelism. When
    // `continuous` is true, HAS_G rows dispatch to the G-M-C trio likelihood.
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
    // `opts.marginalize_maternal` (default true) selects
    // compute_ll_single_marginalized() for the has_g == 0 rows; passing
    // ModelOptions{false} selects the legacy plug-in.  Trio rows are
    // unaffected either way.
    //
    // Each call builds its own maternal Gauss-Jacobi rules.  A caller that
    // sweeps Ne itself should go through estimate_continuous(), which keeps one
    // rule cache alive for the whole scan; see the class comment for the
    // measured difference (19 s vs 1.69 s on the 236-pair shipped cohort).
    static double compute_global_ll_continuous(double ne,
                                               const std::vector<PairData>& data,
                                               const LogFactorial& lf,
                                               int threads = 1,
                                               ModelOptions opts = {});

    // Real-valued Ne optimizer for the continuous model.
    // Phase 1: coarse integer scan; Phase 2: golden-section refinement.
    // Same per-call rule caching caveat as compute_global_ll_continuous().
    static double find_optimal_ne_continuous(const std::vector<PairData>& data,
                                             const LogFactorial& lf,
                                             int min_ne, int max_ne,
                                             int threads = 1,
                                             ModelOptions opts = {});

    // Full estimate with real-valued Ne for the continuous model.  Owns the
    // scan-lifetime maternal rule cache shared by the coarse sweep, the
    // refinement and both CI walks.
    static Result estimate_continuous(const std::vector<PairData>& data,
                                      int min_ne = 1, int max_ne = 200,
                                      int threads = 1,
                                      ModelOptions opts = {});

    // Full estimate (point + 95% CI by profile likelihood, threshold = -1.92).
    // For discrete model only (integer Ne scan); `opts` is forwarded to
    // estimate_continuous() when `continuous` is true and silently ignored
    // otherwise, because the discrete model already marginalises the mother.
    // The CLI, not this function, reports that no-op: with the switch on by
    // default a warning here would fire on every `--model discrete` run.
    static Result estimate(const std::vector<PairData>& data,
                           int min_ne = 1, int max_ne = 200,
                           int threads = 1, bool continuous = false,
                           ModelOptions opts = {});

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
     *     sampling-noise correction : s_i  = p_m (1 - p_m) / (m_dp - 1)
     *                                       + p_c (1 - p_c) / (c_dp - 1)
     *     numerator                 : Sigma_i (d_i - s_i)
     *     denominator               : Sigma_i p_m_i (1 - p_m_i)
     *     normalised drift variance : V    = numerator / denominator
     *
     * The (d-1) denominators are deliberate: p_hat (1 - p_hat) / (d - 1) is
     * unbiased for p (1 - p) / d under binomial sampling, so s_i estimates the
     * read-noise contribution to d_i without the 1/d shrinkage a plain /d would
     * leave in.  Rows with m_dp <= 1 or c_dp <= 1 are excluded from the
     * aggregation because the plug-in variance is not identifiable at depth 1.
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
    * @brief Score every candidate Ne under both the MMLE marginal log-
    *        likelihood and the Kimura sum-of-squared-residuals metric.
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
    * `kimura_norm_ssr = ssr / ssr_min`. The continuous model uses `step`;
    * the discrete model evaluates and reports each integer Ne only. The
    * grid is clamped so that `step > 0` and `min_ne >= 1`.
     */
    static std::vector<NeProfileRow>
    compute_ne_profile(const std::vector<PairData>& data,
                       const LogFactorial& lf,
                       double min_ne, double max_ne, double step,
                       int threads, bool continuous,
                       ModelOptions opts = {});

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
                                         int min_family_sites,
                                         ModelOptions opts = {});

    // Estimate Ne for all families (embarrassingly parallel across families).
    static std::vector<FamilyResult>
    estimate_all_families(const std::vector<FamilyData>& families,
                          int min_ne, int max_ne,
                          int min_family_sites,
                          int threads,
                          ModelOptions opts = {});

    // Per-family Kimura cross-check. When locus_id is available, bootstrap
    // resampling uses whole locus clusters so all child observations at that
    // locus remain together; missing locus IDs fall back to row clusters.
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
