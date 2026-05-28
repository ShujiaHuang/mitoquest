/**
 * @file ne_estimate.cpp
 * @brief Implementation of the `mitoquest ne-estimate` subcommand.
 *
 * @author Shujia Huang (hshujia@qq.com)
 * @date 2026-05-28
 */
#include "ne_estimate.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <thread>

#include "external/thread_pool.h"
#include "io/utils.h"           // ngslib::is_readable

// ---------------------------------------------------------------------
// LogFactorial implementation lives in src/log_factorial.cpp; we only
// expose `mitoquest::LogFactorial` via the type alias declared in
// ne_estimate.h.
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Numerical helpers
// ---------------------------------------------------------------------

double NeEstimator::log_sum_exp_pair(double a, double b) {
    if (std::isinf(a) && a < 0) return b;
    if (std::isinf(b) && b < 0) return a;
    if (a > b) {
        return a + std::log1p(std::exp(b - a));
    }
    return b + std::log1p(std::exp(a - b));
}

// ---------------------------------------------------------------------
// Core likelihoods
// ---------------------------------------------------------------------

// Per-pair log-likelihood:
//
//   p0      ~ Beta(alpha = m_alt + 1, beta = m_ref + 1)   (uniform prior)
//   k       ~ Binomial(Ne, p0)
//   <=>  k  ~ BetaBinomial(Ne, alpha, beta)        (analytic marginalisation)
//   c_alt  ~ Binomial(c_dp, k / Ne)
//
//   logL(Ne) = log Sigma_{k=0..Ne} [ logBB(k|Ne,a,b) + logBin(c_alt|c_dp, k/Ne) ]
//
// The inner sum has Ne + 1 terms; we accumulate it in log-space with
// log_sum_exp_pair to avoid underflow when individual terms are
// extremely small (which is the common case at high depth).
double NeEstimator::compute_ll_single(const PairData& pd, int ne,
                                      const LogFactorial& lf) {
    if (ne < 1) return -std::numeric_limits<double>::infinity();

    // Maternal posterior on the true VAF: Beta(alpha, beta) given a
    // Beta(1,1) (=Uniform) prior and m_alt successes / m_ref failures.
    const double alpha = static_cast<double>(pd.m_ad_alt) + 1.0;
    const double beta  = static_cast<double>(pd.m_dp - pd.m_ad_alt) + 1.0;

    double total = -std::numeric_limits<double>::infinity();
    for (int k = 0; k <= ne; ++k) {
        const double log_bb  = lf.log_betabinom_pmf(ne, k, alpha, beta);
        const double p1      = static_cast<double>(k) / static_cast<double>(ne);
        const double log_seq = lf.log_binomial_pmf(pd.c_dp, pd.c_ad_alt, p1);
        const double log_jt  = log_bb + log_seq;
        total = log_sum_exp_pair(total, log_jt);
    }
    return total;
}

double NeEstimator::compute_global_ll(int ne,
                                      const std::vector<PairData>& data,
                                      const LogFactorial& lf) {
    double total = 0.0;
    for (const PairData& pd : data) {
        total += compute_ll_single(pd, ne, lf);
    }
    return total;
}

double NeEstimator::compute_global_ll_parallel(int ne,
                                               const std::vector<PairData>& data,
                                               const LogFactorial& lf,
                                               int threads) {
    if (threads <= 1 || data.size() < 64) {
        return compute_global_ll(ne, data, lf);
    }

    ThreadPool pool(static_cast<size_t>(threads));
    const size_t n = data.size();
    const size_t chunk = (n + threads - 1) / threads;

    std::vector<std::future<double>> futures;
    futures.reserve(static_cast<size_t>(threads));
    for (size_t start = 0; start < n; start += chunk) {
        const size_t end = std::min(start + chunk, n);
        futures.emplace_back(pool.submit([&, start, end, ne]() {
            double s = 0.0;
            for (size_t i = start; i < end; ++i) {
                s += compute_ll_single(data[i], ne, lf);
            }
            return s;
        }));
    }
    double total = 0.0;
    for (auto& f : futures) total += f.get();
    return total;
}

// Brute-force scan across [min_ne, max_ne].
//
// We tried an integer golden-section pre-bracket, but the discrete
// log-likelihood is NOT unimodal: when true_ne is small, fitted Ne
// values that are integer multiples of the truth (Ne = 2*true_ne,
// 3*true_ne, ...) align with the observed VAF grid and produce
// secondary local maxima.  Golden-section then converges to one of
// those side maxima and misses the true peak.  A full integer scan
// is O(max_ne) global-LL evaluations and is plenty fast for the
// typical range max_ne <= a few hundred.
int NeEstimator::find_optimal_ne(const std::vector<PairData>& data,
                                 const LogFactorial& lf,
                                 int min_ne, int max_ne, int threads) {
    if (min_ne < 1)          min_ne = 1;
    if (max_ne < min_ne)     max_ne = min_ne;

    int    best_ne = min_ne;
    double best_ll = compute_global_ll_parallel(min_ne, data, lf, threads);
    for (int ne = min_ne + 1; ne <= max_ne; ++ne) {
        const double ll = compute_global_ll_parallel(ne, data, lf, threads);
        if (ll > best_ll) {
            best_ll = ll;
            best_ne = ne;
        }
    }
    return best_ne;
}

// 95% profile-likelihood confidence interval threshold:
//   under regularity, -2 * (logL(Ne) - logL_max) ~~ chi2_1
// so the 95% CI is { Ne : logL(Ne) >= logL_max - chi2_{1, 0.95}/2 }
//             with  chi2_{1, 0.95} / 2 = 3.841 / 2 ~~ 1.92.
static constexpr double kProfileLLThresholdDelta = 1.92;

NeEstimator::Result NeEstimator::estimate(const std::vector<PairData>& data,
                                          int min_ne, int max_ne, int threads) {
    Result r;
    r.n_pairs = data.size();
    if (data.empty()) {
        throw std::runtime_error("[ne-estimate] No transmission pairs to fit.");
    }
    if (min_ne < 1)            min_ne = 1;
    if (max_ne < min_ne)       max_ne = min_ne;

    const int cache_size = required_cache_size(data, max_ne);
    LogFactorial lf(cache_size);

    r.ne          = find_optimal_ne(data, lf, min_ne, max_ne, threads);
    r.max_log_lik = compute_global_ll_parallel(r.ne, data, lf, threads);

    const double thr = r.max_log_lik - kProfileLLThresholdDelta;

    // Walk leftward for ci_low; flag when we run off the search boundary.
    r.ci_low = r.ne;
    r.ci_low_clipped = true;
    for (int ne = r.ne - 1; ne >= min_ne; --ne) {
        const double ll = compute_global_ll_parallel(ne, data, lf, threads);
        if (ll < thr) {
            r.ci_low = ne + 1;
            r.ci_low_clipped = false;
            break;
        }
        r.ci_low = ne;
        if (ne == min_ne) {
            r.ci_low_clipped = true;
            break;
        }
    }

    // Walk rightward for ci_high.
    r.ci_high = r.ne;
    r.ci_high_clipped = true;
    for (int ne = r.ne + 1; ne <= max_ne; ++ne) {
        const double ll = compute_global_ll_parallel(ne, data, lf, threads);
        if (ll < thr) {
            r.ci_high = ne - 1;
            r.ci_high_clipped = false;
            break;
        }
        r.ci_high = ne;
        if (ne == max_ne) {
            r.ci_high_clipped = true;
            break;
        }
    }

    return r;
}

int NeEstimator::required_cache_size(const std::vector<PairData>& data,
                                     int max_ne) {
    int max_n = (max_ne > 0) ? max_ne : 1;
    for (const PairData& pd : data) {
        if (pd.c_dp > max_n) max_n = pd.c_dp;
        if (pd.m_dp > max_n) max_n = pd.m_dp;
    }
    return max_n;
}

// ---------------------------------------------------------------------
// TSV input loader
// ---------------------------------------------------------------------

namespace {

// Locate column index by case-sensitive name; throws when missing.
int col_index(const std::vector<std::string>& cols, const std::string& name) {
    for (size_t i = 0; i < cols.size(); ++i) {
        if (cols[i] == name) return static_cast<int>(i);
    }
    throw std::runtime_error("[ne-estimate] Required column missing in TSV: " + name);
}

std::vector<std::string> split_tab(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == '\t') { out.push_back(cur); cur.clear(); }
        else cur.push_back(c);
    }
    out.push_back(cur);
    return out;
}

}  // namespace

std::vector<NeEstimator::PairData>
NeEstimator::load_pairs(const std::string& tsv_path,
                        double min_vaf, double max_vaf) {
    std::ifstream in_file(tsv_path);
    if (!in_file.is_open()) {
        throw std::runtime_error("[ne-estimate] Failed to open input TSV: " + tsv_path);
    }

    std::string line;
    // Skip leading provenance comment lines beginning with '#'.
    std::vector<std::string> header_cols;
    while (std::getline(in_file, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '#') continue;        // comment line written by `trans-prep`
        header_cols = split_tab(line);
        break;
    }
    if (header_cols.empty()) {
        throw std::runtime_error("[ne-estimate] TSV header missing.");
    }

    const int idx_m_dp     = col_index(header_cols, "MOTHER_DP");
    const int idx_m_ad_alt = col_index(header_cols, "MOTHER_AD_ALT");
    const int idx_m_vaf    = col_index(header_cols, "MOTHER_VAF");
    const int idx_c_dp     = col_index(header_cols, "CHILD_DP");
    const int idx_c_ad_alt = col_index(header_cols, "CHILD_AD_ALT");
    const int idx_qc       = col_index(header_cols, "QC");

    std::vector<PairData> data;
    while (std::getline(in_file, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        std::vector<std::string> tk = split_tab(line);
        if (static_cast<int>(tk.size()) <= idx_qc) continue;
        if (tk[idx_qc] != "PASS") continue;

        double m_vaf;
        try {
            m_vaf = std::stod(tk[idx_m_vaf]);
        } catch (const std::exception&) {
            continue;
        }
        if (m_vaf < min_vaf || m_vaf > max_vaf) continue;

        try {
            PairData pd;
            pd.m_dp     = std::stoi(tk[idx_m_dp]);
            pd.m_ad_alt = std::stoi(tk[idx_m_ad_alt]);
            pd.c_dp     = std::stoi(tk[idx_c_dp]);
            pd.c_ad_alt = std::stoi(tk[idx_c_ad_alt]);
            if (pd.m_dp <= 0 || pd.c_dp <= 0) continue;
            if (pd.m_ad_alt < 0 || pd.m_ad_alt > pd.m_dp) continue;
            if (pd.c_ad_alt < 0 || pd.c_ad_alt > pd.c_dp) continue;
            data.push_back(pd);
        } catch (const std::exception&) {
            continue;
        }
    }
    return data;
}

// ---------------------------------------------------------------------
// Wonnapinij / Kimura cross-check (optional)
// ---------------------------------------------------------------------
//
// For a single transmission generation, the Kimura distribution
// parameter b satisfies
//
//     1 - b = Var(p_child - p_mother) / (p_mother * (1 - p_mother))
//
// (Wonnapinij 2008, Eq. 4).  Naive plug-in of the read-frequency point
// estimates inflates the numerator by the binomial sampling noise.
// The Wonnapinij 2010 sampling-error correction subtracts that noise:
//
//     d_i  =  (p_c_i - p_m_i)^2
//     s_i  =  p_m_i (1 - p_m_i) / c_dp_i  +  p_c_i (1 - p_c_i) / m_dp_i
//     V    =  Sigma_i (d_i - s_i)  /  Sigma_i p_m_i (1 - p_m_i)
//     b    =  1 - V
//     Ne_kimura  =  1 / (1 - b)               (single generation, g = 1)
//
// We clip b to (eps, 1 - eps) to avoid divide-by-zero on degenerate
// cohorts.  The value reported here is *only* a sanity cross-check
// against the deCODE 2024 Cell paper -- the Beta-Binomial MLE remains
// the primary estimator.
NeEstimator::KimuraCheck
NeEstimator::compute_kimura_check(const std::vector<PairData>& data) {
    KimuraCheck out;
    out.computed = true;
    out.n_informative = 0;

    double num = 0.0;   // Sigma_i (d_i - s_i)
    double den = 0.0;   // Sigma_i p_m_i (1 - p_m_i)

    for (const PairData& pd : data) {
        if (pd.m_dp <= 0 || pd.c_dp <= 0) continue;
        const double pm = static_cast<double>(pd.m_ad_alt) / static_cast<double>(pd.m_dp);
        const double pc = static_cast<double>(pd.c_ad_alt) / static_cast<double>(pd.c_dp);
        const double w  = pm * (1.0 - pm);
        if (w <= 0.0) continue;             // mother homoplasmic -> no info

        const double d = (pc - pm) * (pc - pm);
        const double s = pm * (1.0 - pm) / static_cast<double>(pd.c_dp)
                       + pc * (1.0 - pc) / static_cast<double>(pd.m_dp);

        num += (d - s);
        den += w;
        ++out.n_informative;
    }

    if (out.n_informative == 0 || den <= 0.0) {
        out.b         = 1.0;
        out.ne_kimura = std::numeric_limits<double>::infinity();
        out.note      = "insufficient informative pairs for Wonnapinij b";
        return out;
    }

    double V = num / den;
    if (V < 0.0) {
        // Sampling correction overshot the observed variance: this means
        // the cohort variance is statistically indistinguishable from
        // pure sampling noise (i.e. very small drift). Report b ~~ 1.
        V = 0.0;
        out.note = "variance after sampling correction <= 0; b clipped to ~1 (Ne -> infinity)";
    }
    double b = 1.0 - V;
    constexpr double kEps = 1e-9;
    if (b < kEps)        { b = kEps;        out.note = "b clipped to lower bound"; }
    else if (b > 1.0 - kEps) { b = 1.0 - kEps; if (out.note.empty()) out.note = "b clipped to upper bound"; }

    out.b         = b;
    out.ne_kimura = 1.0 / (1.0 - b);
    return out;
}

// ---------------------------------------------------------------------
// Constructors / CLI
// ---------------------------------------------------------------------

NeEstimator::NeEstimator(int argc, char* argv[]) {
    _parse_args(argc, argv);
}

NeEstimator::NeEstimator(Config config) : _config(std::move(config)) {
    if (_config.input_tsv.empty()) {
        throw std::invalid_argument("[ne-estimate] input_tsv is required.");
    }
    if (_config.min_ne < 1)       _config.min_ne = 1;
    if (_config.max_ne < _config.min_ne) _config.max_ne = _config.min_ne;
    if (_config.threads < 1)      _config.threads = 1;
    // kimura_check defaults to false; callers may flip it explicitly.
}

void NeEstimator::usage() {
    std::cerr << "Usage: mitoquest ne-estimate [options] -i <pairs.tsv>\n\n"
                 "Description:\n"
                 "  Estimate the mitochondrial DNA bottleneck size (Ne) from mother-child\n"
                 "  transmission pairs using a Beta-Binomial maximum-likelihood framework.\n"
                 "\n"
                 "  Maternal true VAF p0 ~ Beta(m_alt + 1, m_ref + 1)  (uniform prior).\n"
                 "  Bottleneck count    k ~ BetaBinom(Ne, alpha, beta).\n"
                 "  Child read counts c_alt ~ Binomial(c_dp, k/Ne).\n"
                 "\n"
                 "  Reports the maximum-likelihood Ne and 95%% profile-likelihood CI\n"
                 "  (LL_max - 1.92).\n"
                 "\nRequired options:\n"
                 "  -i, --input     FILE   Input transmission pairs TSV produced by\n"
                 "                         `mitoquest trans-prep`.\n"
                 "\nOptional options:\n"
                 "  -o, --output    FILE   JSON output file (default: stdout).\n"
                 "      --min-vaf   FLOAT  Lower maternal VAF gate, inclusive [0.10].\n"
                 "      --max-vaf   FLOAT  Upper maternal VAF gate, inclusive [0.90].\n"
                 "      --min-ne    INT    Smallest Ne value to consider [1].\n"
                 "      --max-ne    INT    Largest Ne value to consider  [200].\n"
                 "  -t, --threads     INT    Worker threads for the inner sum [1].\n"
                 "      --cross-check NAME   Optional secondary estimator alongside the\n"
                 "                           Beta-Binomial MLE.  Supported value: `kimura`,\n"
                 "                           which computes the Wonnapinij b and the implied\n"
                 "                           Ne (single-generation approximation).\n"
                 "  -h, --help               Print this help message.\n\n"
                 "Notes:\n"
                 "  * Sites with maternal VAF near 0 or 1 carry virtually no information\n"
                 "    about Ne, hence the default 0.10 - 0.90 window.\n"
                 "  * Confidence-interval bounds are flagged with `_clipped = true` in the\n"
                 "    JSON output when they hit the search boundary.\n"
                 "  * The Kimura cross-check is approximate; the Beta-Binomial MLE is the\n"
                 "    exact discrete-Wright-Fisher likelihood for a single transmission and\n"
                 "    is the reported primary estimate.\n\n"
              << "Version: " << MITOQUEST_VERSION << "\n"
              << std::endl;
}

void NeEstimator::_parse_args(int argc, char* argv[]) {
    _config.input_tsv.clear();
    _config.output_file.clear();
    _config.min_vaf      = 0.10;
    _config.max_vaf      = 0.90;
    _config.min_ne       = 1;
    _config.max_ne       = 200;
    _config.threads      = 1;
    _config.kimura_check = false;

    _cmdline_string = "#mitoquest_ne_estimate_command=";
    for (int i = 0; i < argc; ++i) {
        _cmdline_string += (i > 0) ? " " + std::string(argv[i]) : std::string(argv[i]);
    }

    static const struct option long_options[] = {
        {"input",        required_argument, 0, 'i'},
        {"output",       required_argument, 0, 'o'},
        {"min-vaf",      required_argument, 0,  1 },
        {"max-vaf",      required_argument, 0,  2 },
        {"min-ne",       required_argument, 0,  3 },
        {"max-ne",       required_argument, 0,  4 },
        {"cross-check",  required_argument, 0,  5 },
        {"threads",      required_argument, 0, 't'},
        {"help",         no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    int c;
    optind = 1;
    while ((c = getopt_long(argc, argv, "i:o:t:h", long_options, &option_index)) != -1) {
        switch (c) {
            case 'i': _config.input_tsv   = optarg;            break;
            case 'o': _config.output_file = optarg;            break;
            case  1 : _config.min_vaf     = std::stod(optarg); break;
            case  2 : _config.max_vaf     = std::stod(optarg); break;
            case  3 : _config.min_ne      = std::stoi(optarg); break;
            case  4 : _config.max_ne      = std::stoi(optarg); break;
            case  5 : {
                const std::string v(optarg);
                if (v == "kimura") {
                    _config.kimura_check = true;
                } else if (v == "none" || v == "off") {
                    _config.kimura_check = false;
                } else {
                    throw std::runtime_error(
                        "[ne-estimate] Unknown --cross-check value: " + v
                        + " (supported: kimura, none).");
                }
                break;
            }
            case 't': _config.threads     = std::stoi(optarg); break;
            case 'h': usage(); std::exit(EXIT_SUCCESS);
            case '?':
            default:  usage(); std::exit(EXIT_FAILURE);
        }
    }

    if (_config.input_tsv.empty()) {
        std::cerr << "[ne-estimate] Missing required option: -i/--input.\n";
        usage();
        std::exit(EXIT_FAILURE);
    }
    if (!ngslib::is_readable(_config.input_tsv.c_str())) {
        throw std::runtime_error("[ne-estimate] Input not readable: " + _config.input_tsv);
    }
    if (_config.min_vaf < 0.0 || _config.max_vaf > 1.0 ||
        _config.min_vaf > _config.max_vaf) {
        throw std::runtime_error("[ne-estimate] Invalid VAF window.");
    }
    if (_config.min_ne < 1) _config.min_ne = 1;
    if (_config.max_ne < _config.min_ne) {
        throw std::runtime_error("[ne-estimate] --max-ne must be >= --min-ne.");
    }
    if (_config.threads < 1) _config.threads = 1;
}

void NeEstimator::_write_json(const Result& r, std::ostream& out) const {
    out << "{\n"
        << "  \"Ne\":              " << r.ne          << ",\n"
        << "  \"CI_95_Low\":       " << r.ci_low      << ",\n"
        << "  \"CI_95_High\":      " << r.ci_high     << ",\n"
        << "  \"CI_Low_Clipped\":  " << (r.ci_low_clipped  ? "true" : "false") << ",\n"
        << "  \"CI_High_Clipped\": " << (r.ci_high_clipped ? "true" : "false") << ",\n"
        << "  \"Pairs_Used\":      " << r.n_pairs     << ",\n"
        << "  \"Max_LogLik\":      " << std::setprecision(8) << r.max_log_lik << ",\n"
        << "  \"Min_VAF\":         " << _config.min_vaf << ",\n"
        << "  \"Max_VAF\":         " << _config.max_vaf << ",\n"
        << "  \"Search_Min_Ne\":   " << _config.min_ne  << ",\n"
        << "  \"Search_Max_Ne\":   " << _config.max_ne;

    if (r.kimura.computed) {
        out << ",\n"
            << "  \"Kimura_Cross_Check\": {\n"
            << "    \"b\":             " << std::setprecision(8) << r.kimura.b         << ",\n"
            << "    \"Ne_Kimura\":     " << std::setprecision(8) << r.kimura.ne_kimura << ",\n"
            << "    \"N_Informative\": " << r.kimura.n_informative << ",\n"
            << "    \"Note\":          \"" << r.kimura.note << "\",\n"
            << "    \"Method\":        \"Wonnapinij 2008/2010 with sampling-error correction; "
            << "single-generation Ne = 1 / (1 - b)\"\n"
            << "  }";
    }

    out << "\n}\n";
}

NeEstimator::Result NeEstimator::run() {
    std::vector<PairData> data = load_pairs(_config.input_tsv,
                                            _config.min_vaf,
                                            _config.max_vaf);
    if (data.empty()) {
        throw std::runtime_error("[ne-estimate] No PASS pairs survived the maternal-VAF "
                                 "filter; nothing to fit.");
    }

    std::cerr << "[ne-estimate] Fitting Ne on " << data.size()
              << " pairs (maternal VAF in [" << _config.min_vaf
              << ", " << _config.max_vaf << "]).\n";

    Result r = estimate(data, _config.min_ne, _config.max_ne, _config.threads);

    if (_config.kimura_check) {
        r.kimura = compute_kimura_check(data);
    }

    std::ofstream out_file;
    std::ostream* out = &std::cout;
    if (!_config.output_file.empty()) {
        out_file.open(_config.output_file);
        if (!out_file.is_open()) {
            throw std::runtime_error("[ne-estimate] Failed to open output: "
                                     + _config.output_file);
        }
        out = &out_file;
    }

    if (!_cmdline_string.empty()) (*out) << _cmdline_string << "\n";
    _write_json(r, *out);

    std::cerr << "[ne-estimate] Optimal Ne = " << r.ne
              << " (95% CI: " << r.ci_low << " - " << r.ci_high
              << (r.ci_low_clipped  ? " [low clipped]"  : "")
              << (r.ci_high_clipped ? " [high clipped]" : "")
              << "), max logL = " << r.max_log_lik
              << " on " << r.n_pairs << " pairs.\n";

    if (r.kimura.computed) {
        std::cerr << "[ne-estimate] Kimura cross-check: b = " << r.kimura.b
                  << ", Ne_kimura = " << r.kimura.ne_kimura
                  << " on " << r.kimura.n_informative << " informative pairs";
        if (!r.kimura.note.empty()) std::cerr << "  [" << r.kimura.note << "]";
        std::cerr << "\n";
    }

    return r;
}
