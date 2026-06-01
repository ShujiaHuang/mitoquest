/**
 * @file variant_qc.h
 * @brief Bayesian quality control for mtDNA variants from VCF files.
 *
 * Distinguishes true heteroplasmic mutations from sequencing artifacts
 * using a Beta-Binomial hypothesis testing framework:
 *
 *   H0: observed alt counts ~ BetaBinom(D, q_alpha, q_beta)   (background noise)
 *   H1: observed alt counts ~ BetaBinom(D, alpha_h1, beta_h1) (true mutation)
 *
 * Per-site background noise parameters (q_alpha, q_beta) are estimated from
 * homozygous-reference samples via MLE.  Global true-mutation parameters
 * (alpha_h1, beta_h1) are estimated from high-confidence variant observations.
 *
 * The posterior probability P(H1 | data) is computed via a log-likelihood
 * ratio test with quality penalties (strand bias, allele quality).
 *
 * @author Shujia Huang (hshujia@qq.com)
 * @date 2025-07-10 (C++ port: 2026-06-01)
 */
#ifndef _VARIANT_QC_H_
#define _VARIANT_QC_H_

#include <getopt.h>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <limits>

#include "log_factorial.h"
#include "version.h"

// Forward declarations for ngslib VCF I/O types (defined in io/*.h)
namespace ngslib {
    class VCFFile;
    class VCFHeader;
    class VCFRecord;
}

class VariantQC {
public:
    /// Type alias to the shared LogFactorial utility.
    using LogFactorial = ::mitoquest::LogFactorial;

    // =================================================================
    // Configuration
    // =================================================================
    struct Config {
        std::string input_vcf;         // input indexed VCF/BCF file
        std::string output_vcf;        // output filtered VCF file
        std::string output_tsv;        // output tabular report (TSV)

        int    max_alt_alleles;        // max ALT alleles per site (default 2)
        int    dp_threshold;           // min DP for per-sample pre-filter (default 100)
        int    hq_threshold;           // min allele quality (AQ) threshold (default 20)
        int    bins;                   // histogram bins for KL divergence (default 100)
        double pi;                     // prior probability of mutation (default 5e-8 * 16569)
        double threshold;              // posterior threshold for calling mutations (default 0.9)
        int    max_iter;               // max iterations for EM-style refinement (default 20)
        double convergence_eps;        // convergence threshold for call-change fraction (default 0.001)
    };

    // =================================================================
    // Per-sample per-site data collected from VCF FORMAT fields
    // =================================================================
    struct SampleInfo {
        std::string sample_name;
        int    ploidy;                 // determined from GT length
        std::vector<int>    gt;        // genotype allele indices
        std::vector<double> af;        // allele frequencies (AF)
        std::vector<int>    ad;        // allele depths (AD)
        std::vector<int>    aq;        // allele quality scores (AQ)
        int    dp;                     // total depth
        std::vector<double> fs;        // Fisher strand values
        std::vector<double> sor;       // strand odds ratio values
        std::vector<double> srf;       // strand ratio factor per allele (computed from SB)
        bool   pre_filtered;           // true if this sample fails pre-filtering
    };

    // =================================================================
    // One result row per sample per site
    // =================================================================
    struct ResultRecord {
        std::string sample_name;
        std::string chrom;
        uint32_t    pos;
        std::string ref;
        std::vector<std::string> alt;  // alleles corresponding to GT
        std::vector<int>    gt;        // original or updated GT
        std::vector<int>    original_gt; // GT as read from VCF (never modified during iteration)
        int    ploidy;
        double q_alpha;                // background noise Beta alpha at this site
        double q_beta;                 // background noise Beta beta at this site
        double p_error;                // mean background error rate at this site
        double posterior;              // max posterior across ALT alleles
        bool   is_mutation;            // posterior > threshold
        bool   pre_filtered;           // carried from SampleInfo
        // Per-allele data (kept for Bayesian filter re-computation)
        std::vector<int>    allele_depths;  // AD values
        std::vector<double> allele_srf;     // SRF values
        std::vector<int>    allele_aq;      // AQ values
    };

    // =================================================================
    // Beta distribution parameters (result of MLE fitting)
    // =================================================================
    struct BetaParams {
        double alpha;
        double beta;
    };

    // =================================================================
    // Public interface
    // =================================================================

    explicit VariantQC(int argc, char* argv[]);
    explicit VariantQC(Config config);   // for tests / library callers
    ~VariantQC() = default;

    /// Run the full QC pipeline: parse -> fit -> filter -> write.
    void run();

    const Config& config() const { return _config; }

    // -----------------------------------------------------------------
    // Pure statistical helpers, exposed for unit tests.
    // -----------------------------------------------------------------

    /// Log of the Beta function: log B(a, b) = lgamma(a) + lgamma(b) - lgamma(a+b).
    static double log_beta_fn(double a, double b);

    /// Beta PDF: f(x | a, b) = x^(a-1) * (1-x)^(b-1) / B(a, b).
    static double beta_pdf(double x, double a, double b);

    /// Fit Beta distribution parameters via MLE using Nelder-Mead on
    /// log-transformed parameters.  Returns (alpha, beta).
    /// Falls back to method-of-moments if MLE fails.
    static BetaParams fit_beta_mle(const std::vector<double>& vafs);

    /// Fit Beta-Binomial parameters via MLE.
    /// Uses VAFs for method-of-moments initialisation, then optimises
    /// the negative log-likelihood using the LogFactorial PMF.
    static BetaParams fit_betabinom_mle(const std::vector<double>& vafs,
                                        const std::vector<int>& k_obs,
                                        const std::vector<int>& n_obs,
                                        const LogFactorial& lf);

    /// Bayesian filter: compute posterior P(H1 | A, D, ...).
    /// Returns a value in [0, 1].
    static double bayesian_filter(
        int A, int D,
        double srf, double hq, int hq_threshold,
        double q_alpha, double q_beta,
        double alpha_h1, double beta_h1,
        double pi);

    /// Multi-sample KL divergence: D_KL(P_empirical || Beta(q_alpha, q_beta)).
    /// Uses a histogram approximation with the given bin edges.
    static double kl_divergence(const std::vector<double>& vafs,
                                double q_alpha, double q_beta,
                                const std::vector<double>& bin_edges);

    /// Binomial PPF (quantile function): smallest k such that
    /// P(X <= k) >= q under Binom(n, p).
    static int binomial_ppf(int n, double p, double q, const LogFactorial& lf);

    /// Check whether a 1-based position falls in a blacklisted region.
    static bool is_blacklisted(uint32_t pos);

    /// Initialise the global blacklist set (called once at startup).
    static void init_blacklist();

private:
    VariantQC(const VariantQC&)            = delete;
    VariantQC& operator=(const VariantQC&) = delete;

    static void usage();
    void _parse_args(int argc, char* argv[]);

    // Pipeline stages
    void _process();
    void _collect_samples_info(const ngslib::VCFRecord& rec,
                               const ngslib::VCFHeader& hdr,
                               const std::vector<std::string>& sample_names,
                               std::vector<SampleInfo>& out);

    std::vector<double> _fetch_background_error_rates(
        const std::vector<SampleInfo>& samples_info);

    void _collect_mutation_vafs(
        const std::vector<SampleInfo>& samples_info,
        double p_error,
        const LogFactorial& lf,
        std::vector<double>& vafs,
        std::vector<int>& k_obs,
        std::vector<int>& n_obs);

    void _iterative_fit_and_call(
        std::vector<ResultRecord>& results,
        const std::vector<double>& vafs,
        const std::vector<int>& k_obs,
        const std::vector<int>& n_obs,
        const LogFactorial& lf);

    void _write_vcf(const std::vector<ResultRecord>& results,
                    const std::map<std::pair<std::string, uint32_t>, double>& pos_kl_div);

    void _write_tsv(const std::vector<ResultRecord>& results);

    // Nelder-Mead simplex optimizer for 2D unconstrained problems.
    // Minimises f(x) where x is a 2-element vector.
    // Returns the best (x1, x2) found.
    struct NMResult { double x1, x2, fval; bool success; };
    static NMResult nelder_mead_2d(
        double (*objective)(double, double, const void*),
        const void* user_data,
        double x0_1, double x0_2,
        int max_iter = 500, double tol = 1e-8);

    Config      _config;
    std::string _cmdline_string;

    // Global blacklist of 1-based positions in known problematic mtDNA regions.
    static std::unordered_set<uint32_t> _blacklist;
    static bool _blacklist_initialized;
};

#endif // _VARIANT_QC_H_
