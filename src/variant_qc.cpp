/**
 * @file variant_qc.cpp
 * @brief Implementation of Bayesian quality control for mtDNA variants.
 *
 * @author Shujia Huang (hshujia@qq.com)
 * @date 2026-06-01
 */
#include "variant_qc.h"
#include "io/vcf.h"
#include "mt_utils.h"

#include <fstream>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <stdexcept>
#include <cstring>

// =====================================================================
// Static member definitions
// =====================================================================
std::unordered_set<uint32_t> VariantQC::_blacklist;
bool VariantQC::_blacklist_initialized = false;

// =====================================================================
// Blacklist initialisation
// =====================================================================
void VariantQC::init_blacklist() {
    if (_blacklist_initialized) return;

    // Known problematic mtDNA regions (1-based, inclusive)
    static const std::pair<uint32_t, uint32_t> regions[] = {
        {299, 317}, {511, 525}, {564, 571}, {952, 955},
        {3105, 3108}, {5895, 5899}, {8268, 8279},
        {13645, 13650}, {16180, 16187}
    };
    for (const auto& [start, end] : regions) {
        for (uint32_t pos = start; pos <= end; ++pos) {
            _blacklist.insert(pos);
        }
    }
    _blacklist_initialized = true;
}

bool VariantQC::is_blacklisted(uint32_t pos) {
    if (!_blacklist_initialized) init_blacklist();
    return _blacklist.count(pos) > 0;
}

// =====================================================================
// Usage
// =====================================================================
void VariantQC::usage() {
    std::cerr
        << "Usage: mitoquest variant-qc [options]\n\n"
        << "  Bayesian quality control for mtDNA variants from VCF files.\n\n"
        << "Required:\n"
        << "  -i, --input-vcf FILE       Input indexed VCF/BCF file\n"
        << "  -o, --output-vcf FILE      Output filtered VCF file\n"
        << "  -t, --output-tsv FILE      Output tabular report (TSV)\n\n"
        << "Optional:\n"
        << "  --max-alt-alleles INT      Max ALT alleles per site [2]\n"
        << "  --dp-threshold INT         Min DP for pre-filter [100]\n"
        << "  --hq-threshold INT         Min allele quality [20]\n"
        << "  --bins INT                 Histogram bins for KL div [100]\n"
        << "  --pi FLOAT                 Prior prob of mutation [5e-8*16569]\n"
        << "  --threshold FLOAT          Posterior threshold [0.9]\n"
        << "  --max-iter INT             Max EM iterations [20]\n"
        << "  --convergence-eps FLOAT    Convergence threshold [0.001]\n"
        << "  -h, --help                 Show this help\n";
}

// =====================================================================
// Argument parsing
// =====================================================================
void VariantQC::_parse_args(int argc, char* argv[]) {
    static struct option long_options[] = {
        {"input-vcf",        required_argument, nullptr, 'i'},
        {"output-vcf",       required_argument, nullptr, 'o'},
        {"output-tsv",       required_argument, nullptr, 't'},
        {"max-alt-alleles",  required_argument, nullptr, 1001},
        {"dp-threshold",     required_argument, nullptr, 1002},
        {"hq-threshold",     required_argument, nullptr, 1003},
        {"bins",             required_argument, nullptr, 1004},
        {"pi",               required_argument, nullptr, 1005},
        {"threshold",        required_argument, nullptr, 1006},
        {"max-iter",         required_argument, nullptr, 1007},
        {"convergence-eps",  required_argument, nullptr, 1008},
        {"help",             no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    // Set defaults
    _config.max_alt_alleles  = 2;
    _config.dp_threshold     = 100;
    _config.hq_threshold     = 20;
    _config.bins             = 100;
    _config.pi               = 5e-8 * 16569;
    _config.threshold        = 0.9;
    _config.max_iter         = 20;
    _config.convergence_eps  = 0.001;

    int opt;
    while ((opt = getopt_long(argc, argv, "i:o:t:h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'i': _config.input_vcf  = optarg; break;
            case 'o': _config.output_vcf = optarg; break;
            case 't': _config.output_tsv = optarg; break;
            case 1001: _config.max_alt_alleles = std::atoi(optarg); break;
            case 1002: _config.dp_threshold    = std::atoi(optarg); break;
            case 1003: _config.hq_threshold    = std::atoi(optarg); break;
            case 1004: _config.bins            = std::atoi(optarg); break;
            case 1005: _config.pi              = std::atof(optarg); break;
            case 1006: _config.threshold       = std::atof(optarg); break;
            case 1007: _config.max_iter        = std::atoi(optarg); break;
            case 1008: _config.convergence_eps = std::atof(optarg); break;
            case 'h': usage(); std::exit(0);
            default:  usage(); throw std::runtime_error("Unknown option.");
        }
    }
    if (_config.input_vcf.empty()) {
        std::cerr << "Error: Input VCF (-i/--input-vcf) is required.\n";
        usage();
        std::exit(EXIT_FAILURE);
    }
}

// =====================================================================
// Constructors
// =====================================================================
VariantQC::VariantQC(int argc, char* argv[]) {
    init_blacklist();
    // Build command-line string for VCF header recording
    _cmdline_string = "mitoquest variant-qc";
    for (int i = 0; i < argc; ++i) {
        _cmdline_string += " ";
        _cmdline_string += argv[i];
    }
    _parse_args(argc, argv);
}

VariantQC::VariantQC(Config config) : _config(std::move(config)) {
    init_blacklist();
    _cmdline_string = "mitoquest variant-qc (programmatic)";
}

// =====================================================================
// Statistical helpers
// =====================================================================

double VariantQC::log_beta_fn(double a, double b) {
    return std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
}

double VariantQC::beta_pdf(double x, double a, double b) {
    if (x <= 0.0 || x >= 1.0) return 0.0;
    if (a <= 0.0 || b <= 0.0) return 0.0;
    double log_pdf = (a - 1.0) * std::log(x) + (b - 1.0) * std::log(1.0 - x) - log_beta_fn(a, b);
    return std::exp(log_pdf);
}

int VariantQC::binomial_ppf(int n, double p, double q, const LogFactorial& lf) {
    if (n <= 0) return 0;
    if (p <= 0.0) return 0;
    if (p >= 1.0) return n;
    if (q <= 0.0) return 0;
    if (q >= 1.0) return n;

    // Binary search for smallest k such that P(X <= k) >= q
    int lo = 0, hi = n;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        // Compute CDF at mid: P(X <= mid) using log-sum-exp
        double log_cdf = -std::numeric_limits<double>::infinity();
        for (int j = 0; j <= mid; ++j) {
            double log_p = lf.log_binomial_pmf(n, j, p);
            double mx = std::max(log_cdf, log_p);
            if (mx == -std::numeric_limits<double>::infinity()) {
                log_cdf = mx;
            } else {
                log_cdf = mx + std::log(std::exp(log_cdf - mx) + std::exp(log_p - mx));
            }
        }
        if (log_cdf >= std::log(q)) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    return lo;
}

// =====================================================================
// Nelder-Mead 2D optimizer
// =====================================================================
VariantQC::NMResult VariantQC::nelder_mead_2d(
    double (*objective)(double, double, const void*),
    const void* user_data,
    double x0_1, double x0_2,
    int max_iter, double tol)
{
    // Simplex: 3 vertices in 2D
    struct Vertex { double x1, x2, f; };

    double step = 0.5;
    Vertex v[3];
    v[0] = {x0_1, x0_2, objective(x0_1, x0_2, user_data)};
    v[1] = {x0_1 + step, x0_2, objective(x0_1 + step, x0_2, user_data)};
    v[2] = {x0_1, x0_2 + step, objective(x0_1, x0_2 + step, user_data)};

    auto sort_simplex = [&]() {
        if (v[0].f > v[1].f) std::swap(v[0], v[1]);
        if (v[0].f > v[2].f) std::swap(v[0], v[2]);
        if (v[1].f > v[2].f) std::swap(v[1], v[2]);
    };

    for (int iter = 0; iter < max_iter; ++iter) {
        sort_simplex();

        // Convergence check
        if (std::abs(v[2].f - v[0].f) < tol) break;

        // Centroid of best two
        double c1 = (v[0].x1 + v[1].x1) / 2.0;
        double c2 = (v[0].x2 + v[1].x2) / 2.0;

        // Reflection
        double r1 = 2.0 * c1 - v[2].x1;
        double r2 = 2.0 * c2 - v[2].x2;
        double fr = objective(r1, r2, user_data);

        if (fr < v[0].f) {
            // Expansion
            double e1 = 3.0 * c1 - 2.0 * v[2].x1;
            double e2 = 3.0 * c2 - 2.0 * v[2].x2;
            double fe = objective(e1, e2, user_data);
            if (fe < fr) { v[2] = {e1, e2, fe}; }
            else         { v[2] = {r1, r2, fr}; }
        } else if (fr < v[1].f) {
            v[2] = {r1, r2, fr};
        } else {
            // Contraction
            double cc1 = (c1 + v[2].x1) / 2.0;
            double cc2 = (c2 + v[2].x2) / 2.0;
            double fc = objective(cc1, cc2, user_data);
            if (fc < v[2].f) {
                v[2] = {cc1, cc2, fc};
            } else {
                // Shrink towards best
                v[1] = {(v[0].x1 + v[1].x1) / 2.0, (v[0].x2 + v[1].x2) / 2.0,
                         objective((v[0].x1 + v[1].x1) / 2.0, (v[0].x2 + v[1].x2) / 2.0, user_data)};
                v[2] = {(v[0].x1 + v[2].x1) / 2.0, (v[0].x2 + v[2].x2) / 2.0,
                         objective((v[0].x1 + v[2].x1) / 2.0, (v[0].x2 + v[2].x2) / 2.0, user_data)};
            }
        }
    }

    sort_simplex();
    NMResult res;
    res.x1 = v[0].x1;
    res.x2 = v[0].x2;
    res.fval = v[0].f;
    res.success = std::isfinite(v[0].f);
    return res;
}

// =====================================================================
// Beta MLE fitting
// =====================================================================

namespace {
struct BetaMLEData {
    double sum_log_x;
    double sum_log_1mx;
    int n;
};

double beta_nll_objective(double log_a, double log_b, const void* data) {
    const auto* d = static_cast<const BetaMLEData*>(data);
    double a = std::exp(log_a);
    double b = std::exp(log_b);
    if (a <= 0.0 || b <= 0.0) return 1e15;
    double nll = -static_cast<double>(d->n) * VariantQC::log_beta_fn(a, b)
                 + (a - 1.0) * d->sum_log_x
                 + (b - 1.0) * d->sum_log_1mx;
    // Negate: we MINIMISE negative log-likelihood
    return -nll;
}
} // anonymous namespace

VariantQC::BetaParams VariantQC::fit_beta_mle(const std::vector<double>& vafs) {
    // Filter valid VAFs
    std::vector<double> valid;
    valid.reserve(vafs.size());
    for (double v : vafs) {
        if (v > 1e-10 && v < 1.0 - 1e-10) valid.push_back(v);
    }
    if (valid.size() < 2) {
        return {2.5, 1680.0}; // default from empirical data
    }

    // Method of moments initialisation
    double sum = 0.0, sum_sq = 0.0;
    for (double v : valid) { sum += v; sum_sq += v * v; }
    double n = static_cast<double>(valid.size());
    double mu = sum / n;
    double var = sum_sq / n - mu * mu;
    if (var <= 0.0) var = 1e-10;
    double max_var = mu * (1.0 - mu);
    if (var >= max_var) var = max_var * 0.99;
    double common = (mu * (1.0 - mu) / var) - 1.0;
    double a_init = std::max(mu * common, 0.01);
    double b_init = std::max((1.0 - mu) * common, 0.01);

    // Prepare data for NLL
    BetaMLEData data{};
    data.n = static_cast<int>(valid.size());
    data.sum_log_x = 0.0;
    data.sum_log_1mx = 0.0;
    for (double v : valid) {
        data.sum_log_x   += std::log(v);
        data.sum_log_1mx += std::log(1.0 - v);
    }

    NMResult res = nelder_mead_2d(beta_nll_objective, &data,
                                   std::log(a_init), std::log(b_init));
    double a = std::exp(res.x1);
    double b = std::exp(res.x2);

    // Bound sanity check
    a = std::max(0.01, std::min(a, 1e6));
    b = std::max(0.01, std::min(b, 1e6));
    return {a, b};
}

// =====================================================================
// Beta-Binomial MLE fitting
// =====================================================================

namespace {
struct BBMLEData {
    const std::vector<int>* k_obs;
    const std::vector<int>* n_obs;
    const mitoquest::LogFactorial* lf;
};

double bb_nll_objective(double log_a, double log_b, const void* data) {
    const auto* d = static_cast<const BBMLEData*>(data);
    double a = std::exp(log_a);
    double b = std::exp(log_b);
    if (a <= 0.0 || b <= 0.0) return 1e15;

    double nll = 0.0;
    int sz = static_cast<int>(d->k_obs->size());
    for (int i = 0; i < sz; ++i) {
        int k = (*d->k_obs)[i];
        int nn = (*d->n_obs)[i];
        double log_p = d->lf->log_betabinom_pmf(nn, k, a, b);
        if (!std::isfinite(log_p)) log_p = -1e10;
        nll -= log_p;
    }
    return nll;
}
} // anonymous namespace

VariantQC::BetaParams VariantQC::fit_betabinom_mle(
    const std::vector<double>& vafs,
    const std::vector<int>& k_obs,
    const std::vector<int>& n_obs,
    const LogFactorial& lf)
{
    if (vafs.empty()) {
        return {0.5, 0.2}; // default fallback
    }

    // Method of moments initialisation from VAFs
    double sum = 0.0, sum_sq = 0.0;
    int n_valid = 0;
    for (double v : vafs) {
        if (v > 0.0 && v < 1.0) { sum += v; sum_sq += v * v; n_valid++; }
    }
    if (n_valid < 1) return {0.5, 0.2};

    double mu = sum / n_valid;
    double var = sum_sq / n_valid - mu * mu;
    if (var <= 0.0) var = 1e-10;
    double max_var = mu * (1.0 - mu);
    if (var >= max_var) var = max_var * 0.99;
    double common = (mu * (1.0 - mu) / var) - 1.0;
    double a_init = std::max(mu * common, 0.01);
    double b_init = std::max((1.0 - mu) * common, 0.01);

    BBMLEData data{k_obs.empty() ? nullptr : &k_obs,
                   n_obs.empty() ? nullptr : &n_obs,
                   &lf};

    NMResult res = nelder_mead_2d(bb_nll_objective, &data,
                                   std::log(a_init), std::log(b_init),
                                   /*max_iter=*/500, /*tol=*/1e-8);
    double a = std::exp(res.x1);
    double b = std::exp(res.x2);
    a = std::max(0.01, std::min(a, 100.0));
    b = std::max(0.01, std::min(b, 1000.0));
    return {a, b};
}

// =====================================================================
// Bayesian filter
// =====================================================================

double VariantQC::bayesian_filter(
    int A, int D,
    double srf, double hq, int hq_threshold,
    double q_alpha, double q_beta,
    double alpha_h1, double beta_h1,
    double pi)
{
    if (D <= 0 || A < 0 || A > D) return 0.0;

    constexpr double hq_k = 10.0;

    // SRF penalty: srf is min/max strand ratio.
    // srf near 1 = balanced strands (no penalty), near 0 = biased (heavy penalty).
    double penalty_srf = srf;

    // HQ sigmoid penalty
    double penalty_hq = (hq < hq_threshold)
        ? 1.0 / (1.0 + std::exp(hq_k * (hq_threshold - hq)))
        : 1.0;

    double log_penalty = std::log(penalty_hq * penalty_srf + 1e-12);

    // Log-likelihood under H0 (Beta-Binomial, without combinatorial
    // constant which cancels in the ratio):
    //   log_L0 = log B(A + q_alpha, D - A + q_beta) - log B(q_alpha, q_beta)
    double log_L0 = log_beta_fn(A + q_alpha, D - A + q_beta) - log_beta_fn(q_alpha, q_beta);

    // Log-likelihood under H1
    double log_L1 = log_beta_fn(A + alpha_h1, D - A + beta_h1) - log_beta_fn(alpha_h1, beta_h1);
    log_L1 += log_penalty;

    // Log prior odds
    pi = std::max(1e-12, std::min(pi, 1.0 - 1e-12));
    double log_prior_odds = std::log(pi) - std::log1p(-pi);

    // Log-posterior odds -> sigmoid
    double log_posterior_odds = log_prior_odds + log_L1 - log_L0;

    double posterior;
    if (std::isfinite(log_posterior_odds)) {
        posterior = 1.0 / (1.0 + std::exp(-log_posterior_odds));
    } else {
        posterior = (log_posterior_odds > 0) ? 1.0 : 0.0;
    }
    return std::max(0.0, std::min(posterior, 1.0));
}

// =====================================================================
// KL divergence
// =====================================================================

double VariantQC::kl_divergence(const std::vector<double>& vafs,
                                 double q_alpha, double q_beta,
                                 const std::vector<double>& bin_edges)
{
    if (vafs.empty()) return 0.0;
    int n_bins = static_cast<int>(bin_edges.size()) - 1;
    if (n_bins <= 0) return 0.0;

    // Histogram (density=True equivalent)
    std::vector<double> px(n_bins, 0.0);
    double delta = bin_edges[1] - bin_edges[0];
    for (double v : vafs) {
        int idx = static_cast<int>((v - bin_edges[0]) / delta);
        if (idx >= 0 && idx < n_bins) px[idx] += 1.0;
    }
    double total = static_cast<double>(vafs.size());
    if (total == 0.0) return 0.0;
    for (double& p : px) p /= (total * delta);

    // Q(x) from Beta distribution
    std::vector<double> qx(n_bins);
    for (int i = 0; i < n_bins; ++i) {
        double center = (bin_edges[i] + bin_edges[i + 1]) / 2.0;
        qx[i] = beta_pdf(center, q_alpha, q_beta);
        if (qx[i] <= 0.0) qx[i] = 1e-10;
    }

    // D_KL(P || Q) = sum P(x) * log(P(x)/Q(x)) * delta
    double kl = 0.0;
    for (int i = 0; i < n_bins; ++i) {
        if (px[i] > 0.0) {
            kl += px[i] * std::log(px[i] / qx[i] + 1e-16) * delta;
        }
    }
    return kl > 0.0 ? kl : 0.0;
}


// =====================================================================
// Sample info collection from VCF record
// =====================================================================
void VariantQC::_collect_samples_info(
    const ngslib::VCFRecord& rec,
    const ngslib::VCFHeader& hdr,
    const std::vector<std::string>& sample_names,
    std::vector<SampleInfo>& out)
{
    int n_samples = static_cast<int>(sample_names.size());
    out.resize(n_samples);

    // Unpack FORMAT fields
    rec.unpack(BCF_UN_FMT);

    // Extract GT
    std::vector<std::vector<int>> genotypes;
    rec.get_genotypes(hdr, genotypes);

    // Extract DP (one value per sample)
    std::vector<int32_t> dp_vals;
    rec.get_format_int(hdr, "DP", dp_vals);

    // Extract AD (variable per sample, stored sample-major)
    std::vector<int32_t> ad_values;
    int ad_per_sample = rec.get_format_int(hdr, "AD", ad_values);

    // Extract AF (variable per sample, stored sample-major)
    std::vector<float> af_values;
    int af_per_sample = rec.get_format_float(hdr, "AF", af_values);

    // Extract AQ (variable per sample, stored sample-major)
    std::vector<int32_t> aq_values;
    int aq_per_sample = rec.get_format_int(hdr, "AQ", aq_values);

    // Extract FS (variable per sample)
    std::vector<float> fs_values;
    int fs_per_sample = rec.get_format_float(hdr, "FS", fs_values);

    // Extract SOR (variable per sample)
    std::vector<float> sor_values;
    int sor_per_sample = rec.get_format_float(hdr, "SOR", sor_values);

    // Extract SB (string per sample)
    std::vector<std::string> sb_values;
    rec.get_format_string(hdr, "SB", sb_values);

    for (int s = 0; s < n_samples; ++s) {
        SampleInfo& si = out[s];
        si.sample_name = sample_names[s];
        si.pre_filtered = false;

        // GT
        if (s < static_cast<int>(genotypes.size())) {
            si.gt = genotypes[s];
            for (auto& g : si.gt) {
                if (g < 0) g = -1; // mark missing
            }
        }
        si.ploidy = static_cast<int>(si.gt.size());

        // Check for missing GT
        bool gt_missing = false;
        for (int g : si.gt) {
            if (g < 0) { gt_missing = true; break; }
        }
        if (gt_missing) {
            si.gt.clear();
            si.ploidy = 0;
        }

        // DP
        si.dp = (s < static_cast<int>(dp_vals.size()) &&
                 dp_vals[s] != ngslib::VCFRecord::INT_MISSING)
                ? dp_vals[s] : -1;

        // AD (per allele, variable length)
        si.ad.clear();
        if (ad_per_sample > 0) {
            int offset = s * ad_per_sample;
            for (int j = 0; j < ad_per_sample; ++j) {
                int val = (offset + j < static_cast<int>(ad_values.size()))
                          ? ad_values[offset + j] : ngslib::VCFRecord::INT_MISSING;
                si.ad.push_back(val == ngslib::VCFRecord::INT_MISSING ? -1 : val);
            }
        }

        // AF
        si.af.clear();
        if (af_per_sample > 0) {
            int offset = s * af_per_sample;
            for (int j = 0; j < af_per_sample; ++j) {
                float val = (offset + j < static_cast<int>(af_values.size()))
                            ? af_values[offset + j] : ngslib::VCFRecord::FLOAT_MISSING;
                si.af.push_back(std::isnan(val) ? -1.0 : static_cast<double>(val));
            }
        }

        // AQ
        si.aq.clear();
        if (aq_per_sample > 0) {
            int offset = s * aq_per_sample;
            for (int j = 0; j < aq_per_sample; ++j) {
                int val = (offset + j < static_cast<int>(aq_values.size()))
                          ? aq_values[offset + j] : ngslib::VCFRecord::INT_MISSING;
                si.aq.push_back(val == ngslib::VCFRecord::INT_MISSING ? -1 : val);
            }
        }

        // FS
        si.fs.clear();
        if (fs_per_sample > 0) {
            int offset = s * fs_per_sample;
            for (int j = 0; j < fs_per_sample; ++j) {
                float val = (offset + j < static_cast<int>(fs_values.size()))
                            ? fs_values[offset + j] : ngslib::VCFRecord::FLOAT_MISSING;
                si.fs.push_back(std::isnan(val) ? -1.0 : static_cast<double>(val));
            }
        }

        // SOR
        si.sor.clear();
        if (sor_per_sample > 0) {
            int offset = s * sor_per_sample;
            for (int j = 0; j < sor_per_sample; ++j) {
                float val = (offset + j < static_cast<int>(sor_values.size()))
                            ? sor_values[offset + j] : ngslib::VCFRecord::FLOAT_MISSING;
                si.sor.push_back(std::isnan(val) ? -1.0 : static_cast<double>(val));
            }
        }

        // Parse SB string to compute SRF
        si.srf.clear();
        if (s < static_cast<int>(sb_values.size()) && !sb_values[s].empty()) {
            std::vector<std::string> allele_sbs = ngslib::split(sb_values[s], ";");
            for (const auto& sb_str : allele_sbs) {
                if (sb_str.empty() || sb_str == ".") {
                    si.srf.push_back(0.5); // default neutral
                    continue;
                }
                std::vector<std::string> parts = ngslib::split(sb_str, ",");
                if (parts.size() >= 2) {
                    double fwd = std::atof(parts[0].c_str());
                    double rev = std::atof(parts[1].c_str());
                    double mx = std::max(fwd, rev) + 1e-10;
                    double mn = std::min(fwd, rev) + 1e-10;
                    si.srf.push_back(mn / mx);
                } else {
                    si.srf.push_back(0.5);
                }
            }
        }

        // Pre-filter: check for missing essential fields
        if (si.dp < 0 || si.ploidy == 0) {
            si.pre_filtered = true;
            continue;
        }
        if (si.af.empty() || si.ad.empty()) {
            si.pre_filtered = true;
            continue;
        }
        for (double a : si.af) {
            if (a < 0.0) { si.pre_filtered = true; break; }
        }
        for (int a : si.ad) {
            if (a < 0) { si.pre_filtered = true; break; }
        }
        if (si.aq.empty()) {
            si.pre_filtered = true;
        }
    }
}

// =====================================================================
// Background error rate collection
// =====================================================================
std::vector<double> VariantQC::_fetch_background_error_rates(
    const std::vector<SampleInfo>& samples_info)
{
    std::vector<double> rates;
    for (const auto& si : samples_info) {
        if (si.ploidy != 1) continue;           // haploid only
        if (si.dp < _config.dp_threshold) continue;
        if (si.pre_filtered) continue;

        // Exclude samples with low allele quality (HQ) - matches Python behavior
        bool hq_ok = true;
        for (int q : si.aq) {
            if (q < _config.hq_threshold) { hq_ok = false; break; }
        }
        if (!hq_ok) continue;

        // Exclude samples with extreme FS or SOR (NUMT indicators)
        bool skip = false;
        for (double f : si.fs) {
            if (f > 100.0) { skip = true; break; }
        }
        if (skip) continue;
        for (double s : si.sor) {
            if (s > 5.0) { skip = true; break; }
        }
        if (skip) continue;

        // Use reference-homozygous samples (GT = [0])
        if (si.gt.size() == 1 && si.gt[0] == 0) {
            if (!si.af.empty() && si.af[0] > 0.0 && si.af[0] < 1.0) {
                rates.push_back(1.0 - si.af[0]);
            }
        }
        // Also include alt-homoplasmic samples (GT = [1], [2], etc.)
        if (si.gt.size() == 1 && si.gt[0] > 0) {
            if (!si.af.empty() && si.af[0] > 0.0 && si.af[0] < 1.0) {
                rates.push_back(1.0 - si.af[0]);
            }
        }
    }
    return rates;
}

// =====================================================================
// Collect high-confidence mutation VAFs
// =====================================================================
void VariantQC::_collect_mutation_vafs(
    const std::vector<SampleInfo>& samples_info,
    double p_error,
    const LogFactorial& lf,
    std::vector<double>& vafs,
    std::vector<int>& k_obs,
    std::vector<int>& n_obs)
{
    const double ep17_alpha = 0.01;

    for (const auto& si : samples_info) {
        if (si.pre_filtered) continue;
        if (si.dp < _config.dp_threshold) continue;

        // Skip reference-only samples
        bool all_ref = true;
        for (int g : si.gt) {
            if (g > 0) { all_ref = false; break; }
        }
        if (all_ref) continue;

        // Check HQ threshold
        bool hq_ok = true;
        for (int q : si.aq) {
            if (q >= 0 && q < _config.hq_threshold) { hq_ok = false; break; }
        }
        if (!hq_ok) continue;

        // Binomial threshold under H0
        int lob = binomial_ppf(si.dp, p_error, 1.0 - ep17_alpha, lf);

        // Collect alt allele observations exceeding the threshold
        int n_alleles = std::min({
            static_cast<int>(si.gt.size()),
            static_cast<int>(si.ad.size()),
            static_cast<int>(si.af.size())
        });
        for (int j = 0; j < n_alleles; ++j) {
            if (si.gt[j] > 0 && si.ad[j] > lob) {
                vafs.push_back(si.af[j]);
                k_obs.push_back(si.ad[j]);
                n_obs.push_back(si.dp);
            }
        }
    }
}

// =====================================================================
// Iterative Beta fit and Bayesian calling
// =====================================================================
void VariantQC::_iterative_fit_and_call(
    std::vector<ResultRecord>& results,
    const std::vector<double>& vafs,
    const std::vector<int>& k_obs,
    const std::vector<int>& n_obs,
    const LogFactorial& lf)
{
    // Balance high-VAF and low-VAF observations for fitting
    std::vector<size_t> high_idx, low_idx;
    for (size_t i = 0; i < vafs.size(); ++i) {
        if (vafs[i] > 0.0 && vafs[i] < 1.0) {
            if (vafs[i] > 0.95) high_idx.push_back(i);
            else low_idx.push_back(i);
        }
    }

    size_t sampling_num = std::min(high_idx.size(), low_idx.size());
    std::vector<double> vafs_sampled;
    std::vector<int> k_sampled, n_sampled;

    if (sampling_num > 0) {
        size_t by_step = std::max(high_idx.size(), low_idx.size()) / sampling_num;
        if (by_step < 1) by_step = 1;

        auto downsample = [&](const std::vector<size_t>& idx, size_t target) {
            std::vector<size_t> selected;
            if (idx.size() > target) {
                for (size_t i = 0; i < idx.size() && selected.size() < target; i += by_step) {
                    selected.push_back(idx[i]);
                }
            } else {
                selected = idx;
            }
            return selected;
        };

        auto sel_low  = downsample(low_idx, sampling_num);
        auto sel_high = downsample(high_idx, sampling_num);

        for (size_t i : sel_low) {
            vafs_sampled.push_back(vafs[i]);
            k_sampled.push_back(k_obs[i]);
            n_sampled.push_back(n_obs[i]);
        }
        for (size_t i : sel_high) {
            vafs_sampled.push_back(vafs[i]);
            k_sampled.push_back(k_obs[i]);
            n_sampled.push_back(n_obs[i]);
        }
    } else {
        vafs_sampled = vafs;
        k_sampled = k_obs;
        n_sampled = n_obs;
    }

    // Initial fit
    BetaParams h1_params = fit_betabinom_mle(vafs_sampled, k_sampled, n_sampled, lf);
    std::cerr << "Training: " << vafs_sampled.size() << " VAFs. "
              << "Initial Beta: alpha=" << h1_params.alpha
              << " beta=" << h1_params.beta << "\n";

    // EM-style iteration: fit -> call -> refit -> ...
    for (int iter = 0; iter < _config.max_iter; ++iter) {
        size_t n_changed = 0;
        for (auto& r : results) {
            if (r.pre_filtered) continue;

            int total_d = 0;
            for (int ad : r.allele_depths) {
                if (ad >= 0) total_d += ad;
            }

            int n_alleles = std::min({
                static_cast<int>(r.allele_depths.size()),
                static_cast<int>(r.allele_srf.size()),
                static_cast<int>(r.allele_aq.size())
            });

            std::vector<double> pps;
            for (int j = 0; j < n_alleles; ++j) {
                double hq_val = (r.allele_aq[j] >= 0)
                    ? static_cast<double>(r.allele_aq[j])
                    : static_cast<double>(_config.hq_threshold);
                double srf_val = (j < static_cast<int>(r.allele_srf.size()))
                    ? r.allele_srf[j] : 0.5;
                double pp = bayesian_filter(
                    r.allele_depths[j], total_d,
                    srf_val, hq_val, _config.hq_threshold,
                    r.q_alpha, r.q_beta,
                    h1_params.alpha, h1_params.beta,
                    _config.pi);
                pps.push_back(pp);
            }

            double new_posterior = pps.empty() ? 0.0 :
                *std::max_element(pps.begin(), pps.end());
            bool new_is_mut = new_posterior > _config.threshold;

            if (new_is_mut != r.is_mutation) n_changed++;

            r.posterior = new_posterior;
            r.is_mutation = new_is_mut;

            // Update GT: zero out alleles with low posterior
            for (int j = 0; j < static_cast<int>(r.gt.size()) &&
                            j < static_cast<int>(pps.size()); ++j) {
                if (pps[j] <= _config.threshold && r.gt[j] > 0) {
                    r.gt[j] = -1; // mark as filtered
                }
            }
        }

        double change_frac = results.empty() ? 0.0 :
            static_cast<double>(n_changed) / static_cast<double>(results.size());
        std::cerr << "  Iteration " << iter + 1 << ": "
                  << n_changed << " calls changed (" << change_frac << ")\n";

        if (change_frac < _config.convergence_eps) {
            std::cerr << "  Converged after " << iter + 1 << " iterations.\n";
            break;
        }

        // Re-collect VAFs from is_mutation=True records for re-fitting
        // Use original_gt (not modified gt) to match Python behavior: all non-ref alleles are included
        std::vector<double> new_vafs;
        std::vector<int> new_k, new_n;
        for (const auto& r : results) {
            if (!r.is_mutation || r.pre_filtered) continue;
            int total_d = 0;
            for (int ad : r.allele_depths) {
                if (ad >= 0) total_d += ad;
            }
            for (size_t j = 0; j < r.original_gt.size() && j < r.allele_depths.size(); ++j) {
                if (r.original_gt[j] > 0 && r.allele_depths[j] > 0 && total_d > 0) {
                    double vaf = static_cast<double>(r.allele_depths[j]) / total_d;
                    if (vaf > 0.0 && vaf < 1.0) {
                        new_vafs.push_back(vaf);
                        new_k.push_back(r.allele_depths[j]);
                        new_n.push_back(total_d);
                    }
                }
            }
        }

        if (new_vafs.size() < 3) {
            std::cerr << "  Too few mutations for re-fitting. Stopping.\n";
            break;
        }

        h1_params = fit_betabinom_mle(new_vafs, new_k, new_n, lf);
    }

    std::cerr << "Final Beta: alpha=" << h1_params.alpha
              << " beta=" << h1_params.beta << "\n";
}

// =====================================================================
// Main pipeline
// =====================================================================
void VariantQC::run() {
    _process();
}

void VariantQC::_process() {
    const int LF_CACHE_SIZE = 100000;
    LogFactorial lf(LF_CACHE_SIZE);

    // Open input VCF
    ngslib::VCFFile vcf_in(_config.input_vcf, "r");
    const ngslib::VCFHeader& hdr = vcf_in.header();
    std::vector<std::string> sample_names = hdr.sample_names();
    int n_samples = hdr.n_samples();

    std::cerr << "[variant-qc] Input: " << _config.input_vcf
              << " (" << n_samples << " samples)\n";

    std::vector<ResultRecord> all_results;
    std::vector<double> global_vafs;
    std::vector<int> global_k, global_n;
    std::map<std::pair<std::string, uint32_t>, double> pos_kl_div;

    ngslib::VCFRecord rec;
    while (vcf_in.read(rec) >= 0) {
        rec.unpack(BCF_UN_ALL);
        std::string chrom = rec.chrom(hdr);
        uint32_t pos = static_cast<uint32_t>(rec.pos() + 1); // 0->1 based
        std::string ref = rec.ref();
        std::vector<std::string> alts = rec.alt();

        // Skip blacklisted positions
        bool blacklisted = false;
        for (size_t i = 0; i < ref.size(); ++i) {
            if (is_blacklisted(pos + static_cast<uint32_t>(i))) {
                blacklisted = true;
                break;
            }
        }
        if (blacklisted) continue;

        // Skip sites with too many ALT alleles
        if (static_cast<int>(alts.size()) > _config.max_alt_alleles) continue;

        std::vector<std::string> ref_alts = {ref};
        ref_alts.insert(ref_alts.end(), alts.begin(), alts.end());

        // Collect per-sample info
        std::vector<SampleInfo> samples_info;
        _collect_samples_info(rec, hdr, sample_names, samples_info);

        // Background error estimation
        std::vector<double> bg_rates = _fetch_background_error_rates(samples_info);
        if (bg_rates.empty()) {
            std::cerr << "[WARNING] No background samples at "
                      << chrom << ":" << pos << ". Skipping.\n";
            continue;
        }

        double p_error = 0.0;
        for (double r : bg_rates) p_error += r;
        p_error /= bg_rates.size();

        BetaParams bg_params = fit_beta_mle(bg_rates);

        // Bin edges for KL divergence
        int n_bins = _config.bins;
        std::vector<double> bin_edges(n_bins + 1);
        for (int i = 0; i <= n_bins; ++i) {
            bin_edges[i] = static_cast<double>(i) / n_bins;
        }

        // Collect mutation VAFs for global fitting
        _collect_mutation_vafs(samples_info, p_error, lf,
                               global_vafs, global_k, global_n);

        // Build result records for each sample
        std::vector<double> site_vafs;
        for (const auto& si : samples_info) {
            ResultRecord rr;
            rr.sample_name = si.sample_name;
            rr.chrom = chrom;
            rr.pos = pos;
            rr.ref = ref;
            rr.ploidy = si.ploidy;
            rr.q_alpha = bg_params.alpha;
            rr.q_beta = bg_params.beta;
            rr.p_error = p_error;
            rr.posterior = 0.0;
            rr.is_mutation = false;
            rr.pre_filtered = si.pre_filtered;

            rr.alt.clear();
            if (!si.gt.empty()) {
                for (int g : si.gt) {
                    if (g >= 0 && g < static_cast<int>(ref_alts.size())) {
                        rr.alt.push_back(ref_alts[g]);
                    } else {
                        rr.alt.push_back(".");
                    }
                }
            } else {
                rr.alt.push_back(".");
            }
            rr.gt = si.gt;
            rr.original_gt = si.gt;  // Preserve original GT for VAF re-collection during iteration
            rr.allele_depths = si.ad;
            rr.allele_srf = si.srf;
            rr.allele_aq = si.aq;

            // Collect site VAFs for KL divergence
            if (!si.pre_filtered) {
                for (double v : si.af) {
                    if (v > 0.0 && v < 1.0) site_vafs.push_back(v);
                }
            }
            all_results.push_back(std::move(rr));
        }

        // KL divergence for this site
        double kl = kl_divergence(site_vafs, bg_params.alpha, bg_params.beta, bin_edges);
        pos_kl_div[{chrom, pos}] = std::isfinite(kl) ? kl : 1000.0;
    }
    vcf_in.close();

    std::cerr << "[variant-qc] Total records: " << all_results.size()
              << ", Global VAFs collected: " << global_vafs.size() << "\n";

    // Iterative fitting and calling
    _iterative_fit_and_call(all_results, global_vafs, global_k, global_n, lf);

    // Write outputs
    _write_vcf(all_results, pos_kl_div);
    _write_tsv(all_results);

    std::cerr << "[variant-qc] Done. Output VCF: " << _config.output_vcf
              << ", TSV: " << _config.output_tsv << "\n";
}

// =====================================================================
// VCF output
// =====================================================================
void VariantQC::_write_vcf(
    const std::vector<ResultRecord>& results,
    const std::map<std::pair<std::string, uint32_t>, double>& pos_kl_div)
{
    // Build lookup: (chrom, pos, sample) -> (gt, posterior, is_mutation)
    struct LookupEntry {
        std::vector<int> gt;
        double posterior;
        bool is_mutation;
    };
    std::map<std::tuple<std::string, uint32_t, std::string>, LookupEntry> lookup;
    for (const auto& r : results) {
        lookup[{r.chrom, r.pos, r.sample_name}] = {r.gt, r.posterior, r.is_mutation};
    }

    // Site-level QC status
    std::map<std::pair<std::string, uint32_t>, bool> site_pass;
    for (const auto& r : results) {
        auto key = std::make_pair(r.chrom, r.pos);
        if (r.is_mutation) site_pass[key] = true;
        else if (site_pass.find(key) == site_pass.end()) site_pass[key] = false;
    }

    // Re-read input VCF and write output
    ngslib::VCFFile vcf_in(_config.input_vcf, "r");
    ngslib::VCFHeader out_hdr = vcf_in.header().copy_header();

    out_hdr.add_header_line("##FORMAT=<ID=PP,Number=1,Type=Float,Description=\"Posterior probability of true genotype\">");
    out_hdr.add_header_line("##FORMAT=<ID=GOOD_CALL,Number=1,Type=String,Description=\"Good call (True/False)\">");
    out_hdr.add_header_line("##FILTER=<ID=PASS,Description=\"Passed mtDNA QC filter\">");
    out_hdr.add_header_line("##FILTER=<ID=BLACKLISTED_SITE,Description=\"Site is in the blacklisted regions\">");
    out_hdr.add_header_line("##FILTER=<ID=LOW_QUALITY,Description=\"Low quality site based on KL divergence\">");
    out_hdr.add_header_line("##variant_qc_command=" + _cmdline_string);

    ngslib::VCFFile vcf_out(_config.output_vcf, out_hdr, "w");

    static const char* TRUE_STR  = "True";
    static const char* FALSE_STR = "False";

    ngslib::VCFRecord rec;
    while (vcf_in.read(rec) >= 0) {
        rec.unpack(BCF_UN_ALL);
        std::string chrom = rec.chrom(out_hdr);
        uint32_t pos = static_cast<uint32_t>(rec.pos() + 1);

        int n_samp = out_hdr.n_samples();
        std::vector<float> pp_values(n_samp, 0.0f);
        std::vector<const char*> gc_values(n_samp, FALSE_STR);

        for (int s = 0; s < n_samp; ++s) {
            std::string sname = out_hdr.sample_names()[s];
            auto it = lookup.find({chrom, pos, sname});
            if (it != lookup.end()) {
                pp_values[s] = static_cast<float>(it->second.posterior);
                gc_values[s] = it->second.is_mutation ? TRUE_STR : FALSE_STR;
            }
        }

        rec.update_format_float(out_hdr, "PP", pp_values.data(), 1);
        rec.update_format_string(out_hdr, "GOOD_CALL", gc_values.data());

        // Set FILTER
        rec.clear_filters(out_hdr);
        if (is_blacklisted(pos)) {
            rec.add_filter(out_hdr, "BLACKLISTED_SITE");
        } else {
            auto sp = site_pass.find({chrom, pos});
            if (sp != site_pass.end() && sp->second) {
                rec.add_filter(out_hdr, "PASS");
            } else {
                rec.add_filter(out_hdr, "LOW_QUALITY");
            }
        }

        // Set QUAL to KL divergence (rounded)
        auto kl_it = pos_kl_div.find({chrom, pos});
        float qual_val = (kl_it != pos_kl_div.end())
            ? static_cast<float>(std::round(kl_it->second)) : 0.0f;
        rec.set_qual(qual_val);

        vcf_out.write(rec);
    }
    vcf_in.close();
    vcf_out.close();
}

// =====================================================================
// TSV output
// =====================================================================
void VariantQC::_write_tsv(const std::vector<ResultRecord>& results) {
    if (_config.output_tsv.empty()) return;

    std::ofstream ofs(_config.output_tsv);
    if (!ofs.is_open()) {
        throw std::runtime_error("Cannot open output TSV: " + _config.output_tsv);
    }

    ofs << "sample\tchrom\tpos\tref\talt\tploidy\t"
        << "q_alpha\tq_beta\tp_error\tposterior\tis_mutation\n";

    for (const auto& r : results) {
        bool has_alt = false;
        std::string alt_str;
        for (size_t j = 0; j < r.alt.size(); ++j) {
            if (j > 0) alt_str += ",";
            alt_str += r.alt[j];
            if (r.alt[j] != r.ref && r.alt[j] != ".") has_alt = true;
        }
        if (!has_alt) continue;
        if (!r.is_mutation) continue;

        ofs << r.sample_name << "\t"
            << r.chrom << "\t"
            << r.pos << "\t"
            << r.ref << "\t"
            << alt_str << "\t"
            << r.ploidy << "\t"
            << format_double(r.q_alpha, 4) << "\t"
            << format_double(r.q_beta, 4) << "\t"
            << format_double(r.p_error, 8) << "\t"
            << format_double(r.posterior, 6) << "\t"
            << (r.is_mutation ? "True" : "False") << "\n";
    }
    ofs.close();
}
