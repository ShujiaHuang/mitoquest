/**
 * @file ne_estimate.cpp
 * @brief Implementation of the `mitoquest ne-estimate` subcommand.
 *
 * @author Shujia Huang (hshujia@qq.com)
 * @date 2026-05-28
 */
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
#include <random>
#include <sstream>
#include <stdexcept>
#include <thread>

#include "ne_estimate.h"
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

// Continuous (Beta-diffusion) per-pair log-likelihood.
//
//   p_child | p_mother  ~  Beta(p_m * (Ne-1), (1-p_m) * (Ne-1))
//   c_alt   | p_child   ~  Binomial(c_dp, p_child)
//
// Marginalising p_child analytically:
//   c_alt  |  p_m  ~  BetaBinomial(c_dp, p_m*(Ne-1), (1-p_m)*(Ne-1))
//
// At Ne=1 the Kimura diffusion degenerates (complete drift to fixation)
// so we fall back to the discrete model which handles k in {0, 1}.
double NeEstimator::compute_ll_single_continuous(const PairData& pd, int ne,
                                                  const LogFactorial& lf) {
    if (ne < 1) return -std::numeric_limits<double>::infinity();
    if (ne == 1) return compute_ll_single(pd, 1, lf);

    const double pm = static_cast<double>(pd.m_ad_alt)
                    / static_cast<double>(pd.m_dp);
    // Mother homoplasmic: no information about drift.
    if (pm <= 0.0 || pm >= 1.0) return 0.0;

    const double ne1 = static_cast<double>(ne - 1);
    return lf.log_betabinom_pmf(pd.c_dp, pd.c_ad_alt,
                                pm * ne1, (1.0 - pm) * ne1);
}

// Real-valued Ne overload for continuous model.
double NeEstimator::compute_ll_single_continuous(const PairData& pd, double ne,
                                                  const LogFactorial& lf) {
    if (ne <= 1.0) return compute_ll_single(pd, 1, lf);

    const double pm = static_cast<double>(pd.m_ad_alt)
                    / static_cast<double>(pd.m_dp);
    if (pm <= 0.0 || pm >= 1.0) return 0.0;

    const double ne1 = ne - 1.0;
    return lf.log_betabinom_pmf(pd.c_dp, pd.c_ad_alt,
                                pm * ne1, (1.0 - pm) * ne1);
}

double NeEstimator::compute_global_ll(int ne,
                                      const std::vector<PairData>& data,
                                      const LogFactorial& lf,
                                      bool continuous) {
    double total = 0.0;
    for (const PairData& pd : data) {
        total += continuous ? compute_ll_single_continuous(pd, ne, lf)
                           : compute_ll_single(pd, ne, lf);
    }
    return total;
}

double NeEstimator::compute_global_ll_parallel(int ne,
                                               const std::vector<PairData>& data,
                                               const LogFactorial& lf,
                                               int threads,
                                               bool continuous) {
    if (threads <= 1 || data.size() < 64) {
        return compute_global_ll(ne, data, lf, continuous);
    }

    ThreadPool pool(static_cast<size_t>(threads));
    const size_t n = data.size();
    const size_t chunk = (n + threads - 1) / threads;

    std::vector<std::future<double>> futures;
    futures.reserve(static_cast<size_t>(threads));
    for (size_t start = 0; start < n; start += chunk) {
        const size_t end = std::min(start + chunk, n);
        futures.emplace_back(pool.submit([&, start, end, ne, continuous]() {
            double s = 0.0;
            for (size_t i = start; i < end; ++i) {
                s += continuous ? compute_ll_single_continuous(data[i], ne, lf)
                               : compute_ll_single(data[i], ne, lf);
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
                                 int min_ne, int max_ne, int threads,
                                 bool continuous) {
    if (min_ne < 1)          min_ne = 1;
    if (max_ne < min_ne)     max_ne = min_ne;

    int    best_ne = min_ne;
    double best_ll = compute_global_ll_parallel(min_ne, data, lf, threads, continuous);
    for (int ne = min_ne + 1; ne <= max_ne; ++ne) {
        const double ll = compute_global_ll_parallel(ne, data, lf, threads, continuous);
        if (ll > best_ll) {
            best_ll = ll;
            best_ne = ne;
        }
    }
    return best_ne;
}

// -----------------------------------------------------------------
// Real-valued continuous model optimization
// -----------------------------------------------------------------

// 95% profile-likelihood confidence interval threshold:
//   under regularity, -2 * (logL(Ne) - logL_max) ~~ chi2_1
// so the 95% CI is { Ne : logL(Ne) >= logL_max - chi2_{1, 0.95}/2 }
//             with  chi2_{1, 0.95} / 2 = 3.841 / 2 ~~ 1.92.
static constexpr double kProfileLLThresholdDelta = 1.92;

double NeEstimator::compute_global_ll_continuous(
    double ne, const std::vector<PairData>& data,
    const LogFactorial& lf, int threads) {
    if (threads <= 1 || data.size() < 64) {
        double total = 0.0;
        for (const PairData& pd : data) {
            total += compute_ll_single_continuous(pd, ne, lf);
        }
        return total;
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
                s += compute_ll_single_continuous(data[i], ne, lf);
            }
            return s;
        }));
    }
    double total = 0.0;
    for (auto& f : futures) total += f.get();
    return total;
}

// Golden-section search to find the maximum of a unimodal function
// on [a, b] to within tolerance `tol`.  Returns the Ne that maximizes
// compute_global_ll_continuous.
static double golden_section_max(
    double a, double b, double tol,
    const std::vector<NeEstimator::PairData>& data,
    const NeEstimator::LogFactorial& lf, int threads) {
    constexpr double phi = 0.6180339887498949;  // (sqrt(5)-1)/2
    double x1 = b - phi * (b - a);
    double x2 = a + phi * (b - a);
    double f1 = NeEstimator::compute_global_ll_continuous(x1, data, lf, threads);
    double f2 = NeEstimator::compute_global_ll_continuous(x2, data, lf, threads);
    while ((b - a) > tol) {
        if (f1 < f2) {
            a  = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + phi * (b - a);
            f2 = NeEstimator::compute_global_ll_continuous(x2, data, lf, threads);
        } else {
            b  = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - phi * (b - a);
            f1 = NeEstimator::compute_global_ll_continuous(x1, data, lf, threads);
        }
    }
    return (a + b) * 0.5;
}

double NeEstimator::find_optimal_ne_continuous(
    const std::vector<PairData>& data, const LogFactorial& lf,
    int min_ne, int max_ne, int threads) {
    if (min_ne < 1) min_ne = 1;
    if (max_ne < min_ne) max_ne = min_ne;

    // Phase 1: coarse integer scan to bracket the peak.
    int    best_int = min_ne;
    double best_ll  = compute_global_ll_continuous(
        static_cast<double>(min_ne), data, lf, threads);
    for (int ne = min_ne + 1; ne <= max_ne; ++ne) {
        const double ll = compute_global_ll_continuous(
            static_cast<double>(ne), data, lf, threads);
        if (ll > best_ll) {
            best_ll  = ll;
            best_int = ne;
        }
    }

    // Phase 2: golden-section refinement in [best_int - 1, best_int + 1].
    const double lo = std::max(static_cast<double>(min_ne),
                               static_cast<double>(best_int) - 1.0);
    const double hi = std::min(static_cast<double>(max_ne),
                               static_cast<double>(best_int) + 1.0);
    return golden_section_max(lo, hi, 0.01, data, lf, threads);
}

NeEstimator::Result NeEstimator::estimate_continuous(
    const std::vector<PairData>& data, int min_ne, int max_ne, int threads) {
    Result r;
    r.n_pairs = data.size();
    if (data.empty()) {
        throw std::runtime_error("[ne-estimate] No transmission pairs to fit.");
    }
    if (min_ne < 1)        min_ne = 1;
    if (max_ne < min_ne)   max_ne = min_ne;

    const int cache_size = required_cache_size(data, max_ne);
    LogFactorial lf(cache_size);

    r.ne          = find_optimal_ne_continuous(data, lf, min_ne, max_ne, threads);
    r.max_log_lik = compute_global_ll_continuous(r.ne, data, lf, threads);

    const double thr = r.max_log_lik - kProfileLLThresholdDelta;
    constexpr double step = 0.01;

    // Walk leftward for ci_low.
    r.ci_low = r.ne;
    r.ci_low_clipped = true;
    for (double ne = r.ne - step; ne >= static_cast<double>(min_ne); ne -= step) {
        const double ll = compute_global_ll_continuous(ne, data, lf, threads);
        if (ll < thr) {
            r.ci_low = ne + step;
            r.ci_low_clipped = false;
            break;
        }
        r.ci_low = ne;
    }
    if (r.ci_low_clipped) r.ci_low = static_cast<double>(min_ne);

    // Walk rightward for ci_high.
    r.ci_high = r.ne;
    r.ci_high_clipped = true;
    for (double ne = r.ne + step; ne <= static_cast<double>(max_ne); ne += step) {
        const double ll = compute_global_ll_continuous(ne, data, lf, threads);
        if (ll < thr) {
            r.ci_high = ne - step;
            r.ci_high_clipped = false;
            break;
        }
        r.ci_high = ne;
    }
    if (r.ci_high_clipped) r.ci_high = static_cast<double>(max_ne);

    return r;
}

NeEstimator::Result NeEstimator::estimate(const std::vector<PairData>& data,
                                          int min_ne, int max_ne, int threads,
                                          bool continuous) {
    // For the continuous model, delegate to the real-valued optimizer.
    if (continuous) {
        return estimate_continuous(data, min_ne, max_ne, threads);
    }

    Result r;
    r.n_pairs = data.size();
    if (data.empty()) {
        throw std::runtime_error("[ne-estimate] No transmission pairs to fit.");
    }
    if (min_ne < 1)            min_ne = 1;
    if (max_ne < min_ne)       max_ne = min_ne;

    const int cache_size = required_cache_size(data, max_ne);
    LogFactorial lf(cache_size);

    const int best_ne = find_optimal_ne(data, lf, min_ne, max_ne, threads, false);
    r.ne          = static_cast<double>(best_ne);
    r.max_log_lik = compute_global_ll_parallel(best_ne, data, lf, threads, false);

    const double thr = r.max_log_lik - kProfileLLThresholdDelta;

    // Walk leftward for ci_low; flag when we run off the search boundary.
    r.ci_low = r.ne;
    r.ci_low_clipped = true;
    for (int ne = best_ne - 1; ne >= min_ne; --ne) {
        const double ll = compute_global_ll_parallel(ne, data, lf, threads, false);
        if (ll < thr) {
            r.ci_low = static_cast<double>(ne + 1);
            r.ci_low_clipped = false;
            break;
        }
        r.ci_low = static_cast<double>(ne);
        if (ne == min_ne) {
            r.ci_low_clipped = true;
            break;
        }
    }

    // Walk rightward for ci_high.
    r.ci_high = r.ne;
    r.ci_high_clipped = true;
    for (int ne = best_ne + 1; ne <= max_ne; ++ne) {
        const double ll = compute_global_ll_parallel(ne, data, lf, threads, false);
        if (ll < thr) {
            r.ci_high = static_cast<double>(ne - 1);
            r.ci_high_clipped = false;
            break;
        }
        r.ci_high = static_cast<double>(ne);
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
        header_cols = ngslib::split(line, "\t");
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
        std::vector<std::string> tk = ngslib::split(line, "\t");
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
//     s_i  =  p_m_i (1 - p_m_i) / m_dp_i  +  p_c_i (1 - p_c_i) / c_dp_i
//                       ^ mother variance is        ^ child variance is
//                         divided by the mother       divided by the child
//                         depth (where p_m was        depth (where p_c was
//                         actually sampled)           actually sampled)
//     V    =  Sigma_i (d_i - s_i)  /  Sigma_i p_m_i (1 - p_m_i)
//     b    =  1 - V
//     Ne_kimura  =  1 / (1 - b)               (single generation, g = 1)
//
// We clip b to (eps, 1 - eps) to avoid divide-by-zero on degenerate
// cohorts.  When `n_bootstrap > 0` we also do a pair-level non-parametric
// bootstrap and return the 2.5/97.5 percentile CI for b and Ne_kimura.
// The value reported here is *only* a sanity cross-check against the
// deCODE 2024 Cell paper -- the Beta-Binomial MLE remains the primary
// estimator.  The two can diverge on real data because:
//
//   * The MLE needs the discrete grid k/Ne to land near observed VAFs;
//     hundreds of concordant heteroplasmic pairs pull Ne *upward*
//     (grid-resolution effect).
//   * The Wonnapinij b is variance-only, so a small number of high-drift
//     outliers (errors / NUMTs / mixed populations) can pull Ne_kimura
//     *downward*.
//
// A single-pair contribution helper, factored out so the bootstrap loop
// can reuse it cheaply.
namespace {

struct PairContribution {
    size_t idx;  // 0-based index back into the original PairData vector
    double pm;   // mother point-estimate VAF
    double pc;   // child  point-estimate VAF
    double w;    // p_m (1 - p_m) -- denominator term
    double r;    // (p_c - p_m)^2 - sampling-noise correction
    bool   informative;  // p_m strictly in (0, 1) so that w > 0
};

// Pre-compute (w_i, r_i, informative_i) for every pair so that the
// bootstrap inner loop is just a vector index + accumulate.
std::vector<PairContribution>
prepare_pair_contributions(const std::vector<NeEstimator::PairData>& data) {
    std::vector<PairContribution> out;
    out.reserve(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        const auto& pd = data[i];
        PairContribution c;
        c.idx = i;
        c.informative = false;
        c.pm = c.pc = c.w = c.r = 0.0;
        if (pd.m_dp <= 0 || pd.c_dp <= 0) { out.push_back(c); continue; }

        c.pm = static_cast<double>(pd.m_ad_alt) / static_cast<double>(pd.m_dp);
        c.pc = static_cast<double>(pd.c_ad_alt) / static_cast<double>(pd.c_dp);
        c.w  = c.pm * (1.0 - c.pm);
        if (c.w <= 0.0) { out.push_back(c); continue; }   // mother homoplasmic

        const double d = (c.pc - c.pm) * (c.pc - c.pm);
        const double s = c.pm * (1.0 - c.pm) / static_cast<double>(pd.m_dp)
                       + c.pc * (1.0 - c.pc) / static_cast<double>(pd.c_dp);
        c.r = d - s;
        c.informative = true;
        out.push_back(c);
    }
    return out;
}

// Compute b (clipped to (eps, 1-eps)) from an aggregated (num, den) pair.
// Returns NaN when den <= 0 (e.g., all-homoplasmic resample).
double b_from_aggregates(double num, double den) {
    constexpr double kEps = 1e-9;
    if (!(den > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    double V = num / den;
    if (V < 0.0)              V = 0.0;       // sampling noise dominates
    double b = 1.0 - V;
    if (b < kEps)             b = kEps;
    else if (b > 1.0 - kEps)  b = 1.0 - kEps;
    return b;
}

// Quantile by linear interpolation on a sorted vector.  Caller must sort
// `xs` ascending and pass a probability in [0, 1].  Returns NaN for
// empty input.
double quantile(const std::vector<double>& xs, double q) {
    if (xs.empty()) return std::numeric_limits<double>::quiet_NaN();
    if (q <= 0.0)   return xs.front();
    if (q >= 1.0)   return xs.back();
    const double pos = q * static_cast<double>(xs.size() - 1);
    const size_t lo  = static_cast<size_t>(std::floor(pos));
    const size_t hi  = static_cast<size_t>(std::ceil(pos));
    if (lo == hi)   return xs[lo];
    const double w   = pos - static_cast<double>(lo);
    return xs[lo] * (1.0 - w) + xs[hi] * w;
}

}  // namespace

NeEstimator::KimuraCheck
NeEstimator::compute_kimura_check(const std::vector<PairData>& data,
                                  int n_bootstrap, uint64_t seed,
                                  double trim_frac, int top_drift_k) {
    KimuraCheck out;
    out.computed = true;
    out.n_informative = 0;

    const std::vector<PairContribution> contribs = prepare_pair_contributions(data);

    double num = 0.0;   // Sigma_i (d_i - s_i)
    double den = 0.0;   // Sigma_i p_m_i (1 - p_m_i)
    for (const auto& c : contribs) {
        if (!c.informative) continue;
        num += c.r;
        den += c.w;
        ++out.n_informative;
    }

    if (out.n_informative == 0 || den <= 0.0) {
        out.b         = 1.0;
        out.ne_kimura = std::numeric_limits<double>::infinity();
        out.note      = "insufficient informative pairs for Wonnapinij b";
        return out;
    }

    constexpr double kEps = 1e-9;
    double V_raw = num / den;
    if (V_raw < 0.0) {
        // Sampling correction overshot the observed variance: this means
        // the cohort variance is statistically indistinguishable from
        // pure sampling noise (i.e. very small drift). Report b ~~ 1.
        out.note = "variance after sampling correction <= 0; b clipped to ~1 (Ne -> infinity)";
    }
    double b = b_from_aggregates(num, den);
    if (b == kEps && out.note.empty())              out.note = "b clipped to lower bound";
    else if (b == 1.0 - kEps && out.note.empty())   out.note = "b clipped to upper bound";

    out.b         = b;
    out.ne_kimura = 1.0 / (1.0 - b);

    // ---- Optional non-parametric bootstrap CI on (b, Ne_kimura) ------
    if (n_bootstrap > 0) {
        out.ci_computed    = true;
        out.n_bootstrap    = n_bootstrap;
        out.bootstrap_seed = seed;

        const size_t n = contribs.size();
        std::vector<double> b_samples;
        b_samples.reserve(static_cast<size_t>(n_bootstrap));

        std::mt19937_64 rng(seed);
        std::uniform_int_distribution<size_t> pick(0, (n > 0) ? n - 1 : 0);

        for (int rep = 0; rep < n_bootstrap; ++rep) {
            double bs_num = 0.0;
            double bs_den = 0.0;
            for (size_t i = 0; i < n; ++i) {
                const PairContribution& c = contribs[pick(rng)];
                if (!c.informative) continue;
                bs_num += c.r;
                bs_den += c.w;
            }
            const double bb = b_from_aggregates(bs_num, bs_den);
            if (std::isnan(bb)) continue;        // all-homoplasmic resample
            b_samples.push_back(bb);
        }

        if (b_samples.size() < 2) {
            // Degenerate bootstrap (essentially all resamples were
            // uninformative).  Fall back to the point estimate.
            out.b_ci_low          = b;
            out.b_ci_high         = b;
            out.ne_kimura_ci_low  = out.ne_kimura;
            out.ne_kimura_ci_high = out.ne_kimura;
            if (!out.note.empty()) out.note += "; ";
            out.note += "bootstrap CI degenerate (too few informative resamples)";
        } else {
            std::sort(b_samples.begin(), b_samples.end());
            const double b_lo = quantile(b_samples, 0.025);
            const double b_hi = quantile(b_samples, 0.975);
            out.b_ci_low          = b_lo;
            out.b_ci_high         = b_hi;
            // Note: b is monotonically increasing in Ne, so the CI for
            // Ne_kimura is just 1/(1-b_lo), 1/(1-b_hi).
            out.ne_kimura_ci_low  = 1.0 / (1.0 - b_lo);
            out.ne_kimura_ci_high = 1.0 / (1.0 - b_hi);
        }
    }

    // ---- Robust trimmed Kimura (drops top `trim_frac` of high-drift pairs)
    // The standard Wonnapinij b is variance-of-moments and is *not* robust
    // to outliers (NUMTs / sequencing errors / mixed populations): a small
    // fraction of high-drift pairs can collapse Ne_kimura by an order of
    // magnitude even when the bulk of pairs are well-behaved.  Sort the
    // informative pairs by their per-pair Wonnapinij contribution
    //     F_i = (d_i - s_i) / w_i
    // (high F_i ~~ high apparent drift), drop the top `trim_frac` of them,
    // and recompute b on the survivors.  When the gap between trimmed and
    // untrimmed is large, the data contain heavy-tailed outliers that the
    // unweighted Wonnapinij is over-fitting.
    if (trim_frac > 0.0 && trim_frac < 1.0 && out.n_informative > 0) {
        std::vector<PairContribution> info;
        info.reserve(out.n_informative);
        for (const auto& c : contribs) {
            if (c.informative) info.push_back(c);
        }
        // Sort by F_i descending (highest-drift first).
        std::sort(info.begin(), info.end(),
                  [](const PairContribution& a, const PairContribution& b) {
                      return (a.r / a.w) > (b.r / b.w);
                  });
        const size_t n_drop = static_cast<size_t>(
            std::floor(trim_frac * static_cast<double>(info.size())));
        const size_t n_keep = (n_drop >= info.size()) ? 0 : info.size() - n_drop;

        double tnum = 0.0, tden = 0.0;
        for (size_t i = n_drop; i < info.size(); ++i) {
            tnum += info[i].r;
            tden += info[i].w;
        }
        out.trimmed_computed = true;
        out.trim_frac        = trim_frac;
        out.n_after_trim     = n_keep;
        if (n_keep > 0 && tden > 0.0) {
            const double b_t = b_from_aggregates(tnum, tden);
            out.b_trimmed         = b_t;
            out.ne_kimura_trimmed = 1.0 / (1.0 - b_t);
        } else {
            out.b_trimmed         = std::numeric_limits<double>::quiet_NaN();
            out.ne_kimura_trimmed = std::numeric_limits<double>::quiet_NaN();
            if (!out.note.empty()) out.note += "; ";
            out.note += "trimmed Kimura: no informative pairs after trimming";
        }
    }

    // ---- Per-pair drift outlier diagnostic (top-K by F_i, descending) -----
    if (top_drift_k > 0 && out.n_informative > 0) {
        std::vector<PairContribution> info;
        info.reserve(out.n_informative);
        for (const auto& c : contribs) {
            if (c.informative) info.push_back(c);
        }
        std::sort(info.begin(), info.end(),
                  [](const PairContribution& a, const PairContribution& b) {
                      return (a.r / a.w) > (b.r / b.w);
                  });
        const size_t k = std::min(static_cast<size_t>(top_drift_k), info.size());
        out.top_drift_outliers.reserve(k);
        for (size_t i = 0; i < k; ++i) {
            const PairContribution& c  = info[i];
            const PairData&         pd = data[c.idx];
            KimuraCheck::DriftOutlier o;
            o.pair_index = c.idx;
            o.m_dp       = pd.m_dp;
            o.m_ad_alt   = pd.m_ad_alt;
            o.c_dp       = pd.c_dp;
            o.c_ad_alt   = pd.c_ad_alt;
            o.m_vaf      = c.pm;
            o.c_vaf      = c.pc;
            o.f_i        = c.r / c.w;
            out.top_drift_outliers.push_back(o);
        }
    }

    return out;
}

// ---------------------------------------------------------------------
// Per-bin observed drift summary (deCODE 2024 Cell Figure 5 reproduction)
// ---------------------------------------------------------------------
//
// For each equal-width maternal-VAF bin in [vaf_low, vaf_high] we report
// three independent estimates of the per-bin drift, all computed on the
// same informative pairs that drive the MLE / Wonnapinij b:
//
//   obs_var      =  mean_i (p_c_i - p_m_i)^2                 (raw)
//   obs_var_corr =  mean_i (d_i - s_i)                       (sampling-corrected)
//   obs_F        =  mean_i (d_i - s_i) / [ p_m_i (1 - p_m_i) ]   (= 1 - b per bin)
//
// Under one-generation Wright-Fisher with effective size Ne the
// theoretical predictions at the *bin center* p_bar are:
//
//   E[obs_var]      ~~ p_bar (1 - p_bar) / Ne
//   E[obs_var_corr] ~~ p_bar (1 - p_bar) / Ne
//   E[obs_F]        =  1 / Ne
//
// The plotting script overlays these theoretical parabolas on the
// observed bin means using the fitted Ne and its 95% CI, reproducing
// the deCODE Figure 5 layout.
std::vector<NeEstimator::BinSimulationRow>
NeEstimator::compute_bin_simulation(const std::vector<PairData>& data,
                                    double vaf_low, double vaf_high, int n_bins) {
    if (n_bins < 1) n_bins = 1;
    if (vaf_high <= vaf_low) {
        // Caller-supplied window is degenerate; force a single bin
        // covering the (eps, 1 - eps) interior and emit an empty row.
        vaf_low  = 0.0;
        vaf_high = 1.0;
    }

    std::vector<BinSimulationRow> rows(static_cast<size_t>(n_bins));
    const double width = (vaf_high - vaf_low) / static_cast<double>(n_bins);
    for (int b = 0; b < n_bins; ++b) {
        rows[b].bin_idx    = b;
        rows[b].bin_low    = vaf_low + width * static_cast<double>(b);
        rows[b].bin_high   = (b == n_bins - 1) ? vaf_high
                                               : vaf_low + width * static_cast<double>(b + 1);
        rows[b].bin_center = 0.5 * (rows[b].bin_low + rows[b].bin_high);
    }

    // Per-bin running sums for mean and (sample) variance via Welford-
    // style accumulation.  We need the SE of mean F_i, so collect
    // sum_F and sum_F_sq for each bin.
    std::vector<double> sum_pm  (n_bins, 0.0);
    std::vector<double> sum_pc  (n_bins, 0.0);
    std::vector<double> sum_d   (n_bins, 0.0);
    std::vector<double> sum_dc  (n_bins, 0.0);   // d - s (corrected)
    std::vector<double> sum_F   (n_bins, 0.0);
    std::vector<double> sum_F_sq(n_bins, 0.0);
    std::vector<size_t> n_in    (n_bins, 0);

    const std::vector<PairContribution> contribs = prepare_pair_contributions(data);
    for (const auto& c : contribs) {
        if (!c.informative) continue;
        if (c.pm < vaf_low || c.pm > vaf_high) continue;
        // Locate bin index; clamp to [0, n_bins - 1] to absorb the upper
        // boundary point.
        int b = static_cast<int>(std::floor((c.pm - vaf_low) / width));
        if (b < 0)        b = 0;
        if (b >= n_bins)  b = n_bins - 1;

        const double dsq = (c.pc - c.pm) * (c.pc - c.pm);
        const double F_i = c.r / c.w;     // (d - s) / w  ; w > 0 by informativeness

        sum_pm  [b] += c.pm;
        sum_pc  [b] += c.pc;
        sum_d   [b] += dsq;
        sum_dc  [b] += c.r;
        sum_F   [b] += F_i;
        sum_F_sq[b] += F_i * F_i;
        ++n_in[b];
    }

    for (int b = 0; b < n_bins; ++b) {
        BinSimulationRow& r = rows[b];
        r.n_pairs = n_in[b];
        if (n_in[b] == 0) continue;
        const double n = static_cast<double>(n_in[b]);
        r.mean_pm      = sum_pm[b] / n;
        r.mean_pc      = sum_pc[b] / n;
        r.obs_var      = sum_d   [b] / n;
        r.obs_var_corr = sum_dc  [b] / n;
        r.obs_F        = sum_F   [b] / n;
        if (n_in[b] >= 2) {
            const double var_F = (sum_F_sq[b] - n * r.obs_F * r.obs_F)
                                  / (n - 1.0);
            r.obs_F_se = (var_F > 0.0) ? std::sqrt(var_F / n) : 0.0;
        } else {
            r.obs_F_se = 0.0;
        }
    }
    return rows;
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
    if (_config.kimura_bootstrap < 0) _config.kimura_bootstrap = 0;
    if (_config.kimura_trim < 0.0 || _config.kimura_trim >= 1.0) _config.kimura_trim = 0.0;
    if (_config.top_drift_k < 0)  _config.top_drift_k = 0;
    if (_config.bin_simulation_n_bins < 1) _config.bin_simulation_n_bins = 1;
    if (_config.model.empty())    _config.model = "continuous";
    // kimura_check defaults to false; callers may flip it explicitly.
}

void NeEstimator::usage() {
    std::cerr << "Usage: mitoquest ne-estimate [options] -i <pairs.tsv>\n\n"
                 "Description:\n"
                 "  Estimate the mitochondrial DNA bottleneck size (Ne) from mother-child\n"
                 "  transmission pairs using maximum-likelihood inference.\n"
                 "\n"
                 "  Default model (continuous / Beta-diffusion):\n"
                 "    p_child | p_mother  ~  Beta(p_m*(Ne-1), (1-p_m)*(Ne-1))\n"
                 "    c_alt   | p_child   ~  Binomial(c_dp, p_child)\n"
                 "  This models mtDNA drift as a continuous Kimura diffusion, which\n"
                 "  naturally accounts for post-bottleneck vegetative segregation.\n"
                 "\n"
                 "  Alternative model (discrete / BetaBinomial-Binomial):\n"
                 "    k ~ BetaBinomial(Ne, m_alt+1, m_ref+1)\n"
                 "    c_alt ~ Binomial(c_dp, k/Ne)\n"
                 "  This restricts child heteroplasmy to the grid {0, 1/Ne, ..., 1}\n"
                 "  and can suffer from upward bias at high sequencing depth.\n"
                 "\n"
                 "  Reports the maximum-likelihood Ne and 95%% profile-likelihood CI\n"
                 "  (LL_max - 1.92).\n"
                 "\nRequired options:\n"
                 "  -i, --input     FILE   Input transmission pairs TSV produced by\n"
                 "                         `mitoquest trans-prep`.\n"
                 "\nOptional options:\n"
                 "  -o, --output    FILE   JSON output file (default: stdout).\n"
                 "      --model     NAME   Likelihood model: `continuous` (default,\n"
                 "                         recommended for mtDNA) or `discrete`.\n"
                 "      --min-vaf   FLOAT  Lower maternal VAF gate, inclusive [0.10].\n"
                 "      --max-vaf   FLOAT  Upper maternal VAF gate, inclusive [0.90].\n"
                 "      --min-ne    INT    Smallest Ne value to consider [1].\n"
                 "      --max-ne    INT    Largest Ne value to consider  [200].\n"
                 "  -t, --threads   INT    Worker threads for the inner sum [1].\n"
                 "      --cross-check NAME   Optional secondary estimator alongside the\n"
                 "                           MLE. Supported value: `kimura`, which\n"
                 "                           computes the Wonnapinij b and the implied\n"
                 "                           Ne (single-generation approximation).\n"
                "      --kimura-bootstrap  INT  Non-parametric bootstrap iterations for the\n"
                 "                              Kimura cross-check 95%% CI.  0 disables [1000].\n"
                 "      --kimura-seed      INT  RNG seed for the Kimura bootstrap [42].\n"
                 "      --kimura-trim    FLOAT  Fraction of highest-drift pairs to drop before\n"
                 "                              recomputing the Kimura b (robust trimmed mean,\n"
                 "                              0.0 disables, recommended 0.10) [0.0].\n"
                 "      --top-drift-k     INT   Emit the top-K highest-drift pairs in the JSON\n"
                 "                              output for outlier inspection (NUMTs / errors).\n"
                 "                              0 disables [0].\n"
                 "      --bin-simulation FILE   Emit a deCODE-style \"Figure 5\" TSV: per\n"
                 "                              maternal-VAF bin observed mean drift vs\n"
                 "                              theoretical p_m(1 - p_m) / Ne under the\n"
                 "                              fitted Ne (and its 95%% CI).  Bin range =\n"
                 "                              [--min-vaf, --max-vaf].\n"
                 "      --bin-simulation-bins INT  Number of equal-width maternal-VAF bins\n"
                 "                                 for --bin-simulation [10].\n"
                 "  -h, --help               Print this help message.\n\n"
                 "Notes:\n"
                 "  * Sites with maternal VAF near 0 or 1 carry virtually no information\n"
                 "    about Ne, hence the default 0.10 - 0.90 window.\n"
                 "  * Confidence-interval bounds are flagged with `_clipped = true` in the\n"
                 "    JSON output when they hit the search boundary.\n"
                 "  * The Kimura cross-check is approximate; the Beta-Binomial MLE is the\n"
                 "    exact discrete-Wright-Fisher likelihood for a single transmission and\n"
                 "    is the reported primary estimate.  On real data the two can diverge\n"
                 "    when many concordant heteroplasmic pairs co-exist with a few high-\n"
                 "    drift outliers (errors / NUMTs / mixed populations).  Use\n"
                 "    --kimura-trim 0.10 (drops the top 10%% of high-drift pairs) and/or\n"
                 "    --top-drift-k 20 (lists the worst pairs) to diagnose this -- see\n"
                 "    release notes for v1.8.2.\n\n"
              << "Version: " << MITOQUEST_VERSION << "\n"
              << std::endl;
}

void NeEstimator::_parse_args(int argc, char* argv[]) {
    _config.input_tsv.clear();
    _config.output_file.clear();
    _config.min_vaf          = 0.10;
    _config.max_vaf          = 0.90;
    _config.min_ne           = 1;
    _config.max_ne           = 200;
    _config.threads          = 1;
    _config.model            = "continuous";
    _config.kimura_check     = false;
    _config.kimura_bootstrap = 1000;
    _config.kimura_seed      = 42;
    _config.kimura_trim      = 0.0;
    _config.top_drift_k      = 0;
    _config.bin_simulation_file.clear();
    _config.bin_simulation_n_bins = 10;

    _cmdline_string = "#mitoquest_ne_estimate_command=";
    for (int i = 0; i < argc; ++i) {
        _cmdline_string += (i > 0) ? " " + std::string(argv[i]) : std::string(argv[i]);
    }

    static const struct option long_options[] = {
        {"input",            required_argument, 0, 'i'},
        {"output",           required_argument, 0, 'o'},
        {"model",            required_argument, 0, 10 },
        {"min-vaf",          required_argument, 0,  1 },
        {"max-vaf",          required_argument, 0,  2 },
        {"min-ne",           required_argument, 0,  3 },
        {"max-ne",           required_argument, 0,  4 },
        {"cross-check",      required_argument, 0,  5 },
        {"kimura-bootstrap", required_argument, 0,  6 },
        {"kimura-seed",      required_argument, 0,  7 },
        {"kimura-trim",      required_argument, 0,  8 },
        {"top-drift-k",      required_argument, 0,  9 },
        {"bin-simulation",       required_argument, 0, 11 },
        {"bin-simulation-bins",  required_argument, 0, 12 },
        {"threads",          required_argument, 0, 't'},
        {"help",             no_argument,       0, 'h'},
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
            case  6 : _config.kimura_bootstrap = std::stoi(optarg); break;
            case  7 : _config.kimura_seed      = static_cast<uint64_t>(
                                                    std::stoull(optarg)); break;
            case  8 : _config.kimura_trim      = std::stod(optarg); break;
            case  9 : _config.top_drift_k      = std::stoi(optarg); break;
            case 10 : {
                const std::string m(optarg);
                if (m == "continuous" || m == "discrete") {
                    _config.model = m;
                } else {
                    throw std::runtime_error(
                        "[ne-estimate] Unknown --model value: " + m
                        + " (supported: continuous, discrete).");
                }
                break;
            }
            case 11 : _config.bin_simulation_file   = optarg;          break;
            case 12 : _config.bin_simulation_n_bins = std::stoi(optarg); break;
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
    if (_config.threads < 1)            _config.threads = 1;
    if (_config.kimura_bootstrap < 0)   _config.kimura_bootstrap = 0;
    if (_config.kimura_trim < 0.0 || _config.kimura_trim >= 1.0) {
        throw std::runtime_error("[ne-estimate] --kimura-trim must be in [0, 1).");
    }
    if (_config.top_drift_k < 0) _config.top_drift_k = 0;
    if (_config.bin_simulation_n_bins < 1) _config.bin_simulation_n_bins = 1;
}

void NeEstimator::_write_json(const Result& r, std::ostream& out) const {
    const bool is_continuous = (_config.model == "continuous");
    // For the continuous model, emit Ne/CI with 2 decimal places;
    // for discrete, emit as integers (no trailing decimals).
    auto emit_ne = [&](double v) {
        if (is_continuous) {
            out << std::fixed << std::setprecision(2) << v;
        } else {
            out << static_cast<int>(std::round(v));
        }
    };
    // Reset float format after emit_ne calls (std::fixed is sticky).
    auto reset_fmt = [&]() {
        out << std::defaultfloat;
    };
    out << "{\n"
        << "  \"Ne\":              "; emit_ne(r.ne);      out << ",\n"
        << "  \"CI_95_Low\":       "; emit_ne(r.ci_low);  out << ",\n"
        << "  \"CI_95_High\":      "; emit_ne(r.ci_high); reset_fmt(); out << ",\n"
        << "  \"CI_Low_Clipped\":  " << (r.ci_low_clipped  ? "true" : "false") << ",\n"
        << "  \"CI_High_Clipped\": " << (r.ci_high_clipped ? "true" : "false") << ",\n"
        << "  \"Pairs_Used\":      " << r.n_pairs << ",\n"
        << "  \"Max_LogLik\":      " << std::setprecision(8) << r.max_log_lik << ",\n"
        << "  \"Model\":           \"" << _config.model << "\",\n"
        << "  \"Min_VAF\":         " << _config.min_vaf << ",\n"
        << "  \"Max_VAF\":         " << _config.max_vaf << ",\n"
        << "  \"Search_Min_Ne\":   " << _config.min_ne  << ",\n"
        << "  \"Search_Max_Ne\":   " << _config.max_ne;

    if (r.kimura.computed) {
        out << ",\n"
            << "  \"Kimura_Cross_Check\": {\n"
            << "    \"b\":             " << std::setprecision(8) << r.kimura.b         << ",\n"
            << "    \"Ne_Kimura\":     " << std::setprecision(8) << r.kimura.ne_kimura << ",\n";

        if (r.kimura.ci_computed) {
            // Helper lambda: emit a finite double as a number, an infinite
            // value as the string "Infinity" (non-strict JSON, but easy
            // for downstream consumers to detect).
            auto emit_num = [&](double v) {
                if (std::isfinite(v)) {
                    out << std::setprecision(8) << v;
                } else {
                    out << "\"" << (v < 0 ? "-Infinity" : "Infinity") << "\"";
                }
            };
            out << "    \"b_CI_95_Low\":         "; emit_num(r.kimura.b_ci_low);          out << ",\n"
                << "    \"b_CI_95_High\":        "; emit_num(r.kimura.b_ci_high);         out << ",\n"
                << "    \"Ne_Kimura_CI_95_Low\": "; emit_num(r.kimura.ne_kimura_ci_low);  out << ",\n"
                << "    \"Ne_Kimura_CI_95_High\":"; emit_num(r.kimura.ne_kimura_ci_high); out << ",\n"
                << "    \"N_Bootstrap\":         " << r.kimura.n_bootstrap     << ",\n"
                << "    \"Bootstrap_Seed\":      " << r.kimura.bootstrap_seed  << ",\n";
        }

        out << "    \"N_Informative\": " << r.kimura.n_informative << ",\n";

        // ---- Trimmed Kimura diagnostic ----
        if (r.kimura.trimmed_computed) {
            auto emit_num = [&](double v) {
                if (std::isfinite(v)) {
                    out << std::setprecision(8) << v;
                } else if (std::isnan(v)) {
                    out << "\"NaN\"";
                } else {
                    out << "\"" << (v < 0 ? "-Infinity" : "Infinity") << "\"";
                }
            };
            out << "    \"Trimmed_Kimura\": {\n"
                << "      \"Trim_Frac\":        " << std::setprecision(4) << r.kimura.trim_frac   << ",\n"
                << "      \"N_After_Trim\":     " << r.kimura.n_after_trim  << ",\n"
                << "      \"b_Trimmed\":        "; emit_num(r.kimura.b_trimmed);                  out << ",\n"
                << "      \"Ne_Kimura_Trimmed\": "; emit_num(r.kimura.ne_kimura_trimmed);          out << "\n"
                << "    },\n";
        }

        // ---- Top-K drift outlier diagnostic ----
        if (!r.kimura.top_drift_outliers.empty()) {
            out << "    \"Top_Drift_Outliers\": [\n";
            for (size_t oi = 0; oi < r.kimura.top_drift_outliers.size(); ++oi) {
                const auto& d = r.kimura.top_drift_outliers[oi];
                out << "      { \"Pair_Index\": " << d.pair_index
                    << ", \"M_DP\": " << d.m_dp
                    << ", \"M_AD_ALT\": " << d.m_ad_alt
                    << ", \"C_DP\": " << d.c_dp
                    << ", \"C_AD_ALT\": " << d.c_ad_alt
                    << ", \"M_VAF\": " << std::setprecision(6) << d.m_vaf
                    << ", \"C_VAF\": " << std::setprecision(6) << d.c_vaf
                    << ", \"F_i\": "   << std::setprecision(8) << d.f_i
                    << " }";
                if (oi + 1 < r.kimura.top_drift_outliers.size()) out << ",";
                out << "\n";
            }
            out << "    ],\n";
        }

        out << "    \"Note\":          \"" << r.kimura.note << "\",\n"
            << "    \"Method\":        \"Wonnapinij 2008/2010 with sampling-error correction; "
            << "single-generation Ne = 1 / (1 - b); 95% CI by non-parametric pair-level bootstrap\"\n"
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
              << ", " << _config.max_vaf << "], model=" << _config.model << ").\n";

    Result r = estimate(data, _config.min_ne, _config.max_ne, _config.threads,
                        _config.model == "continuous");

    if (_config.kimura_check) {
        r.kimura = compute_kimura_check(data,
                                        _config.kimura_bootstrap,
                                        _config.kimura_seed,
                                        _config.kimura_trim,
                                        _config.top_drift_k);
    }

    // ---------------------------------------------------------------
    // Optional: deCODE-style "Figure 5" per-bin observed-vs-simulated
    // drift TSV.  Computed on the same pair set used by the MLE; the
    // theoretical curves p_m(1 - p_m) / Ne are derived from the fitted
    // Ne and (when available) its 95% profile-likelihood CI bounds.
    // ---------------------------------------------------------------
    if (!_config.bin_simulation_file.empty()) {
        const std::vector<BinSimulationRow> bins = compute_bin_simulation(
            data, _config.min_vaf, _config.max_vaf, _config.bin_simulation_n_bins);

        std::ofstream sim_out(_config.bin_simulation_file);
        if (!sim_out.is_open()) {
            throw std::runtime_error(
                "[ne-estimate] Failed to open --bin-simulation output: "
                + _config.bin_simulation_file);
        }

        // Commented-header metadata so downstream Python loaders can
        // recover the fitted Ne / CI without re-parsing the JSON.
        if (!_cmdline_string.empty()) sim_out << _cmdline_string << "\n";
        sim_out << "#mitoquest_version="   << MITOQUEST_VERSION   << "\n"
                << "#model="               << _config.model       << "\n"
                << "#min_vaf="             << _config.min_vaf     << "\n"
                << "#max_vaf="             << _config.max_vaf     << "\n"
                << "#n_bins="              << _config.bin_simulation_n_bins << "\n"
                << "#n_pairs_used="        << r.n_pairs           << "\n"
                << "#fitted_ne="           << std::setprecision(8) << r.ne      << "\n"
                << "#fitted_ne_ci_low="    << std::setprecision(8) << r.ci_low  << "\n"
                << "#fitted_ne_ci_high="   << std::setprecision(8) << r.ci_high << "\n"
                << "#max_log_lik="         << std::setprecision(8) << r.max_log_lik << "\n";
        if (r.kimura.computed) {
            sim_out << "#kimura_b="        << std::setprecision(8) << r.kimura.b         << "\n"
                    << "#kimura_ne="       << std::setprecision(8) << r.kimura.ne_kimura << "\n";
            if (r.kimura.ci_computed) {
                sim_out << "#kimura_ne_ci_low="  << std::setprecision(8) << r.kimura.ne_kimura_ci_low  << "\n"
                        << "#kimura_ne_ci_high=" << std::setprecision(8) << r.kimura.ne_kimura_ci_high << "\n";
            }
        }

        // TSV header + rows.  All quantities are computed at the bin
        // *center* p_bar; the plotting script can re-derive the per-pair
        // theoretical curve from the data column.
        sim_out << "bin_idx\tbin_low\tbin_high\tbin_center\tn_pairs"
                   "\tmean_pm\tmean_pc\tobs_var\tobs_var_corr"
                   "\tobs_F\tobs_F_se"
                   "\texpected_var_at_ne\texpected_var_at_ne_ci_low"
                   "\texpected_var_at_ne_ci_high"
                   "\texpected_F_at_ne\texpected_F_at_ne_ci_low"
                   "\texpected_F_at_ne_ci_high\n";
        sim_out << std::setprecision(8);
        const double ne_pt   = r.ne;
        // Note: smaller Ne -> larger F = 1/Ne and larger expected
        // variance.  Hence the "expected_*_ci_low" column uses
        // ci_high (looser bound on drift) and "expected_*_ci_high"
        // uses ci_low (tighter bound on drift).
        const double ne_lo_f = (r.ci_high > 0) ? r.ci_high : ne_pt;
        const double ne_hi_f = (r.ci_low  > 0) ? r.ci_low  : ne_pt;
        for (const auto& b : bins) {
            const double pbar    = b.bin_center;
            const double pbar1mp = pbar * (1.0 - pbar);
            const double exp_var_pt = (ne_pt   > 0) ? pbar1mp / ne_pt   : 0.0;
            const double exp_var_lo = (ne_lo_f > 0) ? pbar1mp / ne_lo_f : 0.0;
            const double exp_var_hi = (ne_hi_f > 0) ? pbar1mp / ne_hi_f : 0.0;
            const double exp_F_pt   = (ne_pt   > 0) ? 1.0 / ne_pt       : 0.0;
            const double exp_F_lo   = (ne_lo_f > 0) ? 1.0 / ne_lo_f     : 0.0;
            const double exp_F_hi   = (ne_hi_f > 0) ? 1.0 / ne_hi_f     : 0.0;
            sim_out << b.bin_idx       << "\t"
                    << b.bin_low       << "\t"
                    << b.bin_high      << "\t"
                    << b.bin_center    << "\t"
                    << b.n_pairs       << "\t"
                    << b.mean_pm       << "\t"
                    << b.mean_pc       << "\t"
                    << b.obs_var       << "\t"
                    << b.obs_var_corr  << "\t"
                    << b.obs_F         << "\t"
                    << b.obs_F_se      << "\t"
                    << exp_var_pt      << "\t"
                    << exp_var_lo      << "\t"
                    << exp_var_hi      << "\t"
                    << exp_F_pt        << "\t"
                    << exp_F_lo        << "\t"
                    << exp_F_hi        << "\n";
        }
        std::cerr << "[ne-estimate] Wrote bin-simulation TSV ("
                  << bins.size() << " bins) to "
                  << _config.bin_simulation_file << "\n";
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
        if (r.kimura.ci_computed) {
            std::cerr << "  [95% CI: Ne_kimura "
                      << r.kimura.ne_kimura_ci_low << " - "
                      << r.kimura.ne_kimura_ci_high
                      << " via " << r.kimura.n_bootstrap << " bootstraps]";
        }
        if (r.kimura.trimmed_computed && std::isfinite(r.kimura.ne_kimura_trimmed)) {
            std::cerr << "\n[ne-estimate] Trimmed Kimura (trim "
                      << r.kimura.trim_frac * 100.0 << "%%): Ne_kimura_trimmed = "
                      << r.kimura.ne_kimura_trimmed
                      << " on " << r.kimura.n_after_trim << " pairs";
        }
        if (!r.kimura.note.empty()) std::cerr << "  [" << r.kimura.note << "]";
        std::cerr << "\n";

        // Warn when MLE and Kimura disagree by more than 3x.
        if (r.kimura.ne_kimura > 0 && std::isfinite(r.kimura.ne_kimura)) {
            const double ratio = r.ne / r.kimura.ne_kimura;
            if (ratio > 3.0 || ratio < 1.0/3.0) {
                std::cerr << "\n*** WARNING *** Ne_MLE (" << r.ne
                          << ") and Ne_Kimura (" << r.kimura.ne_kimura
                          << ") disagree by >3x.\n"
                          << "    This usually indicates the data contain high-drift outlier pairs\n"
                          << "    (NUMTs / sequencing errors / mixed populations) that collapse the\n"
                          << "    variance-of-moments Kimura estimator but do not affect the MLE.\n"
                          << "    Recommendation: re-run with --kimura-trim 0.10 --top-drift-k 20\n"
                          << "    to inspect the outlier pairs and compare the trimmed Ne_Kimura\n"
                          << "    against the MLE.\n\n";
            }
        }
    }

    return r;
}
