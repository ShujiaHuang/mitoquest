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
#include <map>
#include <mutex>
#include <numeric>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "ne_estimate.h"
#include "external/thread_pool.h"
#include "io/utils.h"           // ngslib::is_readable
#include "version.h"

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
double NeEstimator::compute_ll_single(const PairData& pd, int ne, const LogFactorial& lf) {
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
// with endpoint masses 1-p_m and p_m. This is distinct from the legacy
// discrete model, which integrates a separate maternal Beta posterior.
static double compute_ll_single_continuous_degenerate(const NeEstimator::PairData& pd, double p_m) {
    if (pd.c_ad_alt == 0) return std::log1p(-p_m);
    if (pd.c_ad_alt == pd.c_dp) return std::log(p_m);
    
    return -std::numeric_limits<double>::infinity();
}

double NeEstimator::compute_ll_single_continuous(const PairData& pd, int ne,
                                                 const LogFactorial& lf) {
    if (ne < 1) return -std::numeric_limits<double>::infinity();
    if (pd.m_dp <= 0 || pd.c_dp <= 0 || pd.m_ad_alt < 0 || pd.m_ad_alt > pd.m_dp
        || pd.c_ad_alt < 0 || pd.c_ad_alt > pd.c_dp) 
    {
        return -std::numeric_limits<double>::infinity();
    }

    const double pm = static_cast<double>(pd.m_ad_alt) / static_cast<double>(pd.m_dp);
    if (ne == 1) return compute_ll_single_continuous_degenerate(pd, pm);
    // Mother homoplasmic: no information about drift.
    if (pm <= 0.0 || pm >= 1.0) return 0.0;

    const double ne1 = static_cast<double>(ne - 1);
    return lf.log_betabinom_pmf(pd.c_dp, pd.c_ad_alt,
                                pm * ne1, (1.0 - pm) * ne1);
}

// Real-valued Ne overload for continuous model.
double NeEstimator::compute_ll_single_continuous(const PairData& pd, double ne,
                                                 const LogFactorial& lf) {
    if (!std::isfinite(ne) || ne < 1.0 || pd.m_dp <= 0 || pd.c_dp <= 0
        || pd.m_ad_alt < 0 || pd.m_ad_alt > pd.m_dp
        || pd.c_ad_alt < 0 || pd.c_ad_alt > pd.c_dp) {
        return -std::numeric_limits<double>::infinity();
    }

    const double pm = static_cast<double>(pd.m_ad_alt) / static_cast<double>(pd.m_dp);
    if (ne == 1.0) return compute_ll_single_continuous_degenerate(pd, pm);
    if (pm <= 0.0 || pm >= 1.0) return 0.0;

    const double ne1 = ne - 1.0;
    return lf.log_betabinom_pmf(pd.c_dp, pd.c_ad_alt,
                                pm * ne1, (1.0 - pm) * ne1);
}

constexpr double kTrioQuadratureTolerance = 1e-10;
constexpr int kTrioInitialQuadratureOrder = 16;
constexpr int kTrioMinimumAdaptiveOrder = 64;
constexpr int kTrioRequiredStableRefinements = 2;
constexpr int kTrioMaxEigenIterations = 100;

double log_beta(double alpha, double beta) {
    return std::lgamma(alpha) + std::lgamma(beta) - std::lgamma(alpha + beta);
}

double log_sum_exp_values(const std::vector<double>& values) {
    double max_log = -std::numeric_limits<double>::infinity();
    for (double value : values) {
        if (value > max_log) max_log = value;
    }
    if (!std::isfinite(max_log)) return max_log;

    double sum = 0.0;
    for (double value : values) {
        if (std::isfinite(value)) sum += std::exp(value - max_log);
    }
    return (sum > 0.0) ? max_log + std::log(sum)
                       : -std::numeric_limits<double>::infinity();
}

struct BetaGaussJacobiRule {
    std::vector<double> nodes;
    std::vector<double> log_weights;
};

// Defined below, after the tridiagonal eigensolver.
BetaGaussJacobiRule make_beta_gauss_jacobi_rule(double alpha, double beta,
                                                 int n_nodes);

// Gauss-Jacobi rule cache.
//
// The quadrature rule for the maternal posterior Beta(posterior_alpha,
// posterior_beta) depends only on those two shape parameters and the node
// count -- not on the child read counts.  Rows that share the same
// (grandmother, mother) read signature therefore reuse one rule per
// refinement order.
//
// Two lifetimes are in use, because the two paths have different keys:
//
//   * trio rows are keyed on (alpha_G + k_M, beta_G + d_M - k_M, order) with
//     alpha_G = p_hat_G (Ne-1) and beta_G = (1-p_hat_G) (Ne-1), so the key
//     grid MOVES with Ne.  Retaining those rules across Ne values has almost
//     no hit rate while growing without bound during the CI sweep, so the
//     trio cache stays scoped to a single global-likelihood evaluation
//     (constructed unbounded, as before).
//   * maternal-marginalised rows are keyed on (k_M+1, d_M-k_M+1, order),
//     which carries NO Ne at all.  One rule therefore serves every Ne in the
//     scan, and hoisting this cache to scan lifetime is what makes the
//     marginalised likelihood affordable as the default: rule construction is
//     the Golub-Welsch tridiagonal eigendecomposition and dominates the
//     per-row cost (measured 97.5% of it at d_M=d_C=2000, where an
//     order-128 rule alone costs ~329 us against ~11 us for all 240
//     Beta-Binomial node evaluations).
//
// A scan-lifetime cache is bounded by `node_budget`: once that many nodes are
// stored, further misses are built locally and not retained, so a cohort with
// tens of thousands of distinct maternal signatures degrades to the old
// rebuild-per-Ne behaviour instead of exhausting memory.  Eviction is
// deliberately absent -- it would change nothing numerically (same key always
// yields the same rule) but would cost hit rate.
struct TrioQuadratureCache {
    struct Key {
        double alpha;
        double beta;
        int order;
        bool operator==(const Key& other) const {
            return order == other.order && alpha == other.alpha && beta == other.beta;
        }
    };
    struct KeyHash {
        size_t operator()(const Key& k) const {
            size_t h = std::hash<double>{}(k.alpha);
            h ^= std::hash<double>{}(k.beta) + 0x9e3779b9 + (h << 6) + (h >> 2);
            h ^= static_cast<size_t>(k.order) + 0x9e3779b9 + (h << 6) + (h >> 2);
            return h;
        }
    };

    std::mutex mu;
    std::unordered_map<Key, BetaGaussJacobiRule, KeyHash> rules;
    size_t node_budget  = 0;   // 0 = unlimited (per-evaluation caches)
    size_t nodes_stored = 0;

    explicit TrioQuadratureCache(size_t node_budget_ = 0) : node_budget(node_budget_) {}

    // Returns the rule for (alpha, beta, order), constructing it at most
    // once.  Construction happens outside the lock so parallel workers do
    // not serialize on rule construction; duplicate constructions are
    // possible only when threads race on the identical key, which is rare.
    //
    // Returns nullptr when the node budget is exhausted; the caller must then
    // build the rule itself.  A non-null pointer stays valid for the lifetime
    // of the cache: std::unordered_map keeps element addresses stable across
    // rehash, and insertion is the only mutation.
    const BetaGaussJacobiRule* try_get(double alpha, double beta, int order) {
        const Key key{alpha, beta, order};
        {
            std::lock_guard<std::mutex> lock(mu);
            const auto it = rules.find(key);
            if (it != rules.end()) return &it->second;
            if (node_budget != 0
                && nodes_stored + static_cast<size_t>(order) > node_budget) {
                return nullptr;
            }
        }
        BetaGaussJacobiRule rule = make_beta_gauss_jacobi_rule(alpha, beta, order);
        std::lock_guard<std::mutex> lock(mu);
        const auto it = rules.find(key);
        if (it != rules.end()) return &it->second;
        if (node_budget != 0 && nodes_stored + static_cast<size_t>(order) > node_budget) {
            return nullptr;
        }
        nodes_stored += static_cast<size_t>(order);

        return &rules.emplace(key, std::move(rule)).first->second;
    }
};

// Node budget for a scan-lifetime maternal cache.  4M nodes x 16 bytes
// (node + log-weight) is ~64 MB, enough to hold every rule of a ~16k-row
// cohort (up to 4 refinement orders per distinct maternal signature) while
// keeping pathological inputs bounded.
constexpr size_t kMaternalScanCacheNodeBudget = 4u << 20;

// Diagonalize a real symmetric tridiagonal matrix with implicit QL steps.
// We only need the first component of every normalized eigenvector: its
// square is the normalized Gauss quadrature weight (Golub-Welsch theorem).
void tridiagonal_eigenpairs_first_components(
    std::vector<double>& diagonal,
    const std::vector<double>& off_diagonal,
    std::vector<double>& first_components) {
    const int n = static_cast<int>(diagonal.size());
    if (n < 1 || static_cast<int>(off_diagonal.size()) != n - 1) {
        throw std::runtime_error("[ne-estimate] Invalid Jacobi quadrature matrix.");
    }

    std::vector<double> e(static_cast<size_t>(n), 0.0);
    for (int i = 0; i < n - 1; ++i) e[static_cast<size_t>(i)] = off_diagonal[static_cast<size_t>(i)];

    first_components.assign(static_cast<size_t>(n), 0.0);
    first_components[0] = 1.0;
    const double epsilon = std::numeric_limits<double>::epsilon();

    for (int left = 0; left < n; ++left) {
        int iterations = 0;
        while (true) {
            int right = left;
            for (; right < n - 1; ++right) {
                const double scale = std::abs(diagonal[static_cast<size_t>(right)])
                                   + std::abs(diagonal[static_cast<size_t>(right + 1)]);
                if (std::abs(e[static_cast<size_t>(right)]) <= epsilon * scale) break;
            }
            if (right == left) break;
            if (++iterations > kTrioMaxEigenIterations) {
                throw std::runtime_error("[ne-estimate] Gauss-Jacobi eigensolver did not converge.");
            }

            double shift = (diagonal[static_cast<size_t>(left + 1)]
                         - diagonal[static_cast<size_t>(left)]) / (2.0 * e[static_cast<size_t>(left)]);
            const double hypotenuse = std::hypot(shift, 1.0);
            shift = diagonal[static_cast<size_t>(right)] - diagonal[static_cast<size_t>(left)]
                  + e[static_cast<size_t>(left)] / (shift + std::copysign(hypotenuse, shift));

            double sine = 1.0;
            double cosine = 1.0;
            double correction = 0.0;
            double rotation_norm = 0.0;
            int index = right - 1;
            for (; index >= left; --index) {
                const double f = sine * e[static_cast<size_t>(index)];
                const double b = cosine * e[static_cast<size_t>(index)];
                rotation_norm = std::hypot(f, shift);
                e[static_cast<size_t>(index + 1)] = rotation_norm;
                if (rotation_norm == 0.0) {
                    diagonal[static_cast<size_t>(index + 1)] -= correction;
                    e[static_cast<size_t>(right)] = 0.0;
                    break;
                }

                sine = f / rotation_norm;
                cosine = shift / rotation_norm;
                shift = diagonal[static_cast<size_t>(index + 1)] - correction;
                rotation_norm = (diagonal[static_cast<size_t>(index)] - shift) * sine + 2.0 * cosine * b;
                correction = sine * rotation_norm;
                diagonal[static_cast<size_t>(index + 1)] = shift + correction;
                shift = cosine * rotation_norm - b;

                const double next_component = first_components[static_cast<size_t>(index + 1)];
                first_components[static_cast<size_t>(index + 1)] =
                    sine * first_components[static_cast<size_t>(index)] + cosine * next_component;
                first_components[static_cast<size_t>(index)] =
                    cosine * first_components[static_cast<size_t>(index)] - sine * next_component;
            }
            if (rotation_norm == 0.0 && index >= left) continue;

            diagonal[static_cast<size_t>(left)] -= correction;
            e[static_cast<size_t>(left)] = shift;
            e[static_cast<size_t>(right)] = 0.0;
        }
    }
}

BetaGaussJacobiRule make_beta_gauss_jacobi_rule(double alpha, double beta, int n_nodes) {
    if (!(alpha > 0.0) || !(beta > 0.0) || n_nodes < 1) {
        throw std::runtime_error("[ne-estimate] Invalid Beta parameters for Gauss-Jacobi quadrature.");
    }
    if (n_nodes == 1) {
        return {{alpha / (alpha + beta)}, {0.0}};
    }

    const double total = alpha + beta;
    std::vector<double> diagonal(static_cast<size_t>(n_nodes), 0.0);
    std::vector<double> off_diagonal(static_cast<size_t>(n_nodes - 1), 0.0);
    diagonal[0] = alpha / total;

    for (int degree = 1; degree < n_nodes; ++degree) {
        const double n = static_cast<double>(degree);
        const double center = 2.0 * n + total - 2.0;
        diagonal[static_cast<size_t>(degree)] = 0.5 * (1.0
            + (alpha - beta) * (total - 2.0) / (center * (center + 2.0)));

        double off_diagonal_squared;
        if (degree == 1) {
            off_diagonal_squared = alpha * beta / (total * total * (total + 1.0));
        } else {
            off_diagonal_squared = n * (n + alpha - 1.0) * (n + beta - 1.0)
                * (n + total - 2.0) / (center * center * (center - 1.0) * (center + 1.0));
        }
        if (!(off_diagonal_squared > 0.0) || !std::isfinite(off_diagonal_squared)) {
            throw std::runtime_error("[ne-estimate] Invalid Gauss-Jacobi recurrence coefficient.");
        }
        off_diagonal[static_cast<size_t>(degree - 1)] = std::sqrt(off_diagonal_squared);
    }

    std::vector<double> first_components;
    tridiagonal_eigenpairs_first_components(diagonal, off_diagonal, first_components);

    std::vector<size_t> order(static_cast<size_t>(n_nodes));
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t left, size_t right) {
        return diagonal[left] < diagonal[right];
    });

    BetaGaussJacobiRule rule;
    rule.nodes.reserve(static_cast<size_t>(n_nodes));
    rule.log_weights.reserve(static_cast<size_t>(n_nodes));
    std::vector<long double> raw_weights;
    raw_weights.reserve(static_cast<size_t>(n_nodes));
    long double weight_sum = 0.0L;
    for (size_t index : order) {
        double node = diagonal[index];
        if (!std::isfinite(node)) {
            throw std::runtime_error("[ne-estimate] Non-finite Gauss-Jacobi node.");
        }
        if (node <= 0.0) node = std::nextafter(0.0, 1.0);
        if (node >= 1.0) node = std::nextafter(1.0, 0.0);
        rule.nodes.push_back(node);

        const long double component = static_cast<long double>(first_components[index]);
        const long double weight = component * component;
        raw_weights.push_back(weight);
        weight_sum += weight;
    }
    if (!(weight_sum > 0.0L) || !std::isfinite(static_cast<double>(weight_sum))) {
        throw std::runtime_error("[ne-estimate] Invalid Gauss-Jacobi quadrature weights.");
    }
    for (long double weight : raw_weights) {
        rule.log_weights.push_back((weight > 0.0L)
            ? static_cast<double>(std::log(weight) - std::log(weight_sum))
            : -std::numeric_limits<double>::infinity());
    }
    return rule;
}

double compute_ll_trio_gauss_jacobi_at_order(const NeEstimator::PairData& pd,
                                              double ne,
                                              const NeEstimator::LogFactorial& lf,
                                              int n_nodes,
                                              TrioQuadratureCache* cache) {
    const double s = ne - 1.0;
    const double pg = static_cast<double>(pd.g_ad_alt) / static_cast<double>(pd.g_dp);
    const double alpha_g = pg * s;
    const double beta_g = (1.0 - pg) * s;
    const double posterior_alpha = alpha_g + static_cast<double>(pd.m_ad_alt);
    const double posterior_beta = beta_g + static_cast<double>(pd.m_dp - pd.m_ad_alt);
    // A miss (no cache handed in, or the cache's node budget is spent) falls
    // back to building the rule here.  Both routes produce the same rule, so
    // the returned log-likelihood is bit-identical either way.
    BetaGaussJacobiRule local_rule;
    const BetaGaussJacobiRule* rule_ptr = (cache != nullptr)
        ? cache->try_get(posterior_alpha, posterior_beta, n_nodes) : nullptr;

    if (rule_ptr == nullptr) {
        local_rule = make_beta_gauss_jacobi_rule(posterior_alpha, posterior_beta, n_nodes);
        rule_ptr = &local_rule;
    }
    const BetaGaussJacobiRule& rule = *rule_ptr;

    std::vector<double> log_terms;
    log_terms.reserve(rule.nodes.size());
    for (size_t i = 0; i < rule.nodes.size(); ++i) {
        const double p_m = rule.nodes[i];
        const double child_ll = lf.log_betabinom_pmf(
            pd.c_dp, pd.c_ad_alt, s * p_m, s * (1.0 - p_m));
        log_terms.push_back(rule.log_weights[i] + child_ll);
    }

    return lf.log_comb(pd.m_dp, pd.m_ad_alt)
         + log_beta(posterior_alpha, posterior_beta)
         - log_beta(alpha_g, beta_g)
         + log_sum_exp_values(log_terms);
}

// Adaptive Gauss-Jacobi order driver, shared by the trio path and the maternal-
// marginalisation path.  `ll_at_order(n)` evaluates the row log-likelihood with
// an n-node rule; the loop starts at kTrioInitialQuadratureOrder, doubles, and
// returns once kTrioRequiredStableRefinements consecutive refinements agree
// within kTrioQuadratureTolerance at an order of at least
// kTrioMinimumAdaptiveOrder, or when the polynomial-exact cap is reached.
//
// Reaching `exact_order` is NOT a correctness guarantee.  The child factor is a
// degree-d_C polynomial in p_M, so a rule with n >= (d_C+1)/2 nodes is exact in
// real arithmetic, but the Golub-Welsch nodes/weights themselves lose double-
// precision fidelity at very high order, so the capped value can be *worse*
// than the converged low-order one.  Measured on the row (d_M,k_M,d_C,k_C) =
// (2000,1000,2000,40) at Ne=100: orders 16-512 all give log E[g] = -102.9067,
// agreeing with an independent 400k-point trapezoid reference, while order 1001
// gives -99.2613 -- off by +3.65 log units.  The cap is therefore a last-resort
// bound, not a target.  Empirically it is never reached at deep coverage (0 of
// 4761 evaluations at d_C=2000 and 0 of 4719 at d_C=500, over Ne in
// {2,30,200}), because the sequence stabilises by order 128; at shallow
// coverage the cap *is* the exact order and is small (16 nodes at d_C=30, 51 at
// d_C=100), where it agrees with a 256-node rule to ~2e-13.
template <typename FnLlAtOrder>
double adaptive_gauss_jacobi_ll(int exact_order, FnLlAtOrder ll_at_order) {
    int order = std::min(kTrioInitialQuadratureOrder, exact_order);
    double previous = ll_at_order(order);
    if (order == exact_order) return previous;

    int stable_refinements = 0;
    while (order < exact_order) {
        const int next_order = std::min(exact_order, order * 2);
        const double current = ll_at_order(next_order);
        if (next_order >= std::min(kTrioMinimumAdaptiveOrder, exact_order)
            && std::abs(current - previous) <= kTrioQuadratureTolerance) {
            ++stable_refinements;
        } else {
            stable_refinements = 0;
        }
        if (next_order == exact_order
            || stable_refinements >= kTrioRequiredStableRefinements) {
            return current;
        }
        previous = current;
        order = next_order;
    }
    return previous;
}

double compute_ll_trio_gauss_jacobi(const NeEstimator::PairData& pd,
                                    double ne,
                                    const NeEstimator::LogFactorial& lf,
                                    TrioQuadratureCache* cache) {
    const int exact_order = (pd.c_dp + 2) / 2;
    return adaptive_gauss_jacobi_ll(exact_order, [&](int n_nodes) {
        return compute_ll_trio_gauss_jacobi_at_order(pd, ne, lf, n_nodes, cache);
    });
}

double compute_ll_trio_degenerate(const NeEstimator::PairData& pd, double pg) {
    const bool all_reference = pd.m_ad_alt == 0 && pd.c_ad_alt == 0;
    const bool all_alternate = pd.m_ad_alt == pd.m_dp && pd.c_ad_alt == pd.c_dp;
    if (all_reference && pg < 1.0) return std::log1p(-pg);
    if (all_alternate && pg > 0.0) return std::log(pg);
    return -std::numeric_limits<double>::infinity();
}

// Trio marginal log-likelihood with an optional per-evaluation quadrature
// rule cache (see TrioQuadratureCache).  The public entry point below passes
// a null cache; the continuous-model optimizer uses a real one.
double compute_ll_trio_continuous_cached(const NeEstimator::PairData& pd,
                                         double ne,
                                         const NeEstimator::LogFactorial& lf,
                                         TrioQuadratureCache* cache) {
    if (pd.has_g == 0) return NeEstimator::compute_ll_single_continuous(pd, ne, lf);
    if (!std::isfinite(ne) || ne < 1.0 || pd.g_dp <= 0
        || pd.g_ad_alt < 0 || pd.g_ad_alt > pd.g_dp) {
        return -std::numeric_limits<double>::infinity();
    }
    const double pg = static_cast<double>(pd.g_ad_alt)
                    / static_cast<double>(pd.g_dp);
    if (ne == 1.0 || pg <= 0.0 || pg >= 1.0) {
        return compute_ll_trio_degenerate(pd, pg);
    }
    return compute_ll_trio_gauss_jacobi(pd, ne, lf, cache);
}

// -----------------------------------------------------------------
// Two-generation M-C pair with the maternal heteroplasmy integrated out
// (the default; --no-maternal-marginalization selects the plug-in).  See the
// declaration comment in ne_estimate.h for the formula, the normalisation
// identities and the Ne == 1 closed form.
// -----------------------------------------------------------------

// Ne == 1 degenerate limit: s = Ne - 1 = 0 collapses the child Beta transition
// onto the two endpoint masses {0, 1} with weights (1 - p_M, p_M), so E[g]
// reduces to E[1 - p_M] or E[p_M] under Beta(k_M+1, d_M-k_M+1), whose mean is
// (k_M+1)/(d_M+2).  The joint normalisation log C(d_M,k_M) +
// log B(k_M+1, d_M-k_M+1) is exactly -log(d_M+1).  e_pm is strictly inside
// (0,1) for every admissible row, so both logs are finite.
double compute_ll_single_marginalized_degenerate(const NeEstimator::PairData& pd) {
    const double e_pm = (static_cast<double>(pd.m_ad_alt) + 1.0)
                      / (static_cast<double>(pd.m_dp) + 2.0);

    const double log_norm = -std::log(static_cast<double>(pd.m_dp) + 1.0);
    if (pd.c_ad_alt == 0)        return log_norm + std::log1p(-e_pm);
    if (pd.c_ad_alt == pd.c_dp)  return log_norm + std::log(e_pm);

    return -std::numeric_limits<double>::infinity();
}

// Fixed-order evaluation.  Structurally identical to
// compute_ll_trio_gauss_jacobi_at_order with the grandmother-informed prior
// replaced by Beta(1,1): the rule weight becomes the maternal posterior
// Beta(k_M+1, d_M-k_M+1) itself, so the `- log_beta(alpha_g, beta_g)` prior
// term drops out (log B(1,1) = 0) and only the posterior normalisation remains.
// The rule parameters do not depend on Ne, so the caller passes a cache whose
// lifetime spans the whole Ne scan (see compute_global_ll_continuous), giving
// a strictly better hit rate than the trio path, whose alpha_g / beta_g move
// with Ne and which therefore keeps a per-evaluation cache.
double compute_ll_single_marginalized_at_order(const NeEstimator::PairData& pd,
                                               double ne,
                                               const NeEstimator::LogFactorial& lf,
                                               int n_nodes,
                                               TrioQuadratureCache* cache) {

    const double s     = ne - 1.0;
    const double alpha = static_cast<double>(pd.m_ad_alt) + 1.0;
    const double beta  = static_cast<double>(pd.m_dp - pd.m_ad_alt) + 1.0;
    // Same local fallback as the trio path; bit-identical on a cache miss.
    BetaGaussJacobiRule local_rule;
    const BetaGaussJacobiRule* rule_ptr = (cache != nullptr)
        ? cache->try_get(alpha, beta, n_nodes) : nullptr;

    if (rule_ptr == nullptr) {
        local_rule = make_beta_gauss_jacobi_rule(alpha, beta, n_nodes);
        rule_ptr = &local_rule;
    }
    const BetaGaussJacobiRule& rule = *rule_ptr;

    std::vector<double> log_terms;
    log_terms.reserve(rule.nodes.size());

    for (size_t i = 0; i < rule.nodes.size(); ++i) {
        const double p_m = rule.nodes[i];
        const double child_ll = lf.log_betabinom_pmf(pd.c_dp, pd.c_ad_alt, s * p_m, s * (1.0 - p_m));
        log_terms.push_back(rule.log_weights[i] + child_ll);
    }

    return lf.log_comb(pd.m_dp, pd.m_ad_alt) + log_beta(alpha, beta) + log_sum_exp_values(log_terms);
}

double compute_ll_single_marginalized_cached(const NeEstimator::PairData& pd,
                                             double ne,
                                             const NeEstimator::LogFactorial& lf,
                                             TrioQuadratureCache* cache) {

    if (!std::isfinite(ne) || ne < 1.0 || pd.m_dp <= 0 || pd.c_dp <= 0
        || pd.m_ad_alt < 0 || pd.m_ad_alt > pd.m_dp
        || pd.c_ad_alt < 0 || pd.c_ad_alt > pd.c_dp) {
        return -std::numeric_limits<double>::infinity();
    }
    if (ne == 1.0) return compute_ll_single_marginalized_degenerate(pd);

    // Deliberately NO homoplasmic shortcut here.  The plug-in path returns the
    // constant 0.0 when p_hat_M is 0 or 1 because a point estimate at the
    // boundary carries no drift information; the whole point of marginalising
    // is that k_M = 0 out of d_M reads still leaves a proper Beta(1, d_M+1)
    // posterior over p_M, against which the child's k_C IS informative.  Such
    // rows are excluded by the default VAF window anyway and only enter under
    // --min-vaf 0, where recovering them is the desired behaviour.
    //
    // No new -infinity risk is introduced: make_beta_gauss_jacobi_rule clamps
    // every node strictly inside (0,1), so s * p_m > 0 for ne > 1 and the
    // Beta-Binomial log-pmf stays finite at all nodes.
    const int exact_order = (pd.c_dp + 2) / 2;
    return adaptive_gauss_jacobi_ll(exact_order, [&](int n_nodes) {
        return compute_ll_single_marginalized_at_order(pd, ne, lf, n_nodes, cache);
    });
}

// Single place that decides which per-row likelihood a continuous-model
// evaluation uses.  Trio rows always marginalise p_M over the grandmother-
// informed posterior, so opts is only consulted for two-generation rows.
//
// The two caches have deliberately different lifetimes (see
// TrioQuadratureCache): `trio_cache` is scoped to one global-likelihood
// evaluation because the trio key moves with Ne, while `maternal_cache` spans
// the whole Ne scan because the maternal key does not.
double compute_row_ll_continuous(const NeEstimator::PairData& pd, double ne,
                                 const NeEstimator::LogFactorial& lf,
                                 TrioQuadratureCache* trio_cache,
                                 TrioQuadratureCache* maternal_cache,
                                 const NeEstimator::ModelOptions& opts) {
    if (pd.has_g == 1) {
        return compute_ll_trio_continuous_cached(pd, ne, lf, trio_cache);
    }
    if (opts.marginalize_maternal) {
        return compute_ll_single_marginalized_cached(pd, ne, lf, maternal_cache);
    }

    return NeEstimator::compute_ll_single_continuous(pd, ne, lf);
}

// -----------------------------------------------------------------
// Three-generation G-M-C trio marginal log-likelihood.
// -----------------------------------------------------------------
//
// With s = Ne - 1, the exact likelihood is the Beta-weighted integral
//
//   C(d_M,k_M) * B(alpha_G+k_M, beta_G+d_M-k_M) / B(alpha_G,beta_G)
//   * E_{p_M ~ Beta(alpha_G+k_M, beta_G+d_M-k_M)}[
//       BetaBin(k_C | d_C, s*p_M, s*(1-p_M))].
//
// The child Beta-Binomial term is a degree-d_C polynomial in p_M.  The
// production path uses adaptive Gauss-Jacobi quadrature for this Beta weight,
// capped at the polynomial-exact order ceil((d_C + 1) / 2).
double NeEstimator::compute_ll_trio_continuous(const PairData& pd, double ne,
                                               const LogFactorial& lf) {
    return compute_ll_trio_continuous_cached(pd, ne, lf, nullptr);
}

// -----------------------------------------------------------------
// Two-generation M-C pair with the maternal heteroplasmy integrated out.
// -----------------------------------------------------------------
//
// Public entry point; the optimizer uses compute_ll_single_marginalized_cached
// with a real quadrature cache.  Validated against the shipped trio path: at
// Ne = 3 with p_hat_G = 0.5 the trio prior Beta(p_hat_G s, (1-p_hat_G) s) is
// exactly Beta(1,1), so compute_ll_trio_continuous() on a row carrying
// g_dp = 2, g_ad_alt = 1 returns precisely this function's value.
double NeEstimator::compute_ll_single_marginalized(const PairData& pd, double ne,
                                                   const LogFactorial& lf) {
    return compute_ll_single_marginalized_cached(pd, ne, lf, nullptr);
}

// -----------------------------------------------------------------
// Gauss-Legendre quadrature for the trio marginal likelihood.
//
// This is retained as a diagnostic implementation of the exact integrand.
// Unlike the production Gauss-Jacobi path, it does not absorb the endpoint
// Beta weight and is not reliable when the posterior Beta density is singular.
//
// The integral over p_M in [0,1] is:
//   int Beta(p_M | alpha_G, beta_G)
//       * Bin(k_M | d_M, p_M)
//       * BetaBin(k_C | d_C, p_M*(Ne-1), (1-p_M)*(Ne-1))
//       dp_M
//
// We map t in [-1,1] to p_M = (1+t)/2 using standard GL nodes/weights.
// -----------------------------------------------------------------
double NeEstimator::compute_ll_trio_quadrature(const PairData& pd, double ne,
                                               const LogFactorial& lf,
                                               int n_nodes) {

    if (pd.has_g == 0) return compute_ll_single_continuous(pd, ne, lf);
    if (!std::isfinite(ne) || ne < 1.0 || n_nodes < 1 || pd.g_dp <= 0
        || pd.g_ad_alt < 0 || pd.g_ad_alt > pd.g_dp) {
        return -std::numeric_limits<double>::infinity();
    }

    const double pg = static_cast<double>(pd.g_ad_alt) / static_cast<double>(pd.g_dp);
    if (ne == 1.0 || pg <= 0.0 || pg >= 1.0) {
        return compute_ll_trio_degenerate(pd, pg);
    }

    const double ne1     = ne - 1.0;
    const double alpha_G = pg * ne1;
    const double beta_G  = (1.0-pg) * ne1;

    // Compute GL nodes and weights on [-1, 1] (Golub-Welsch / Newton).
    std::vector<double> nodes(static_cast<size_t>(n_nodes));
    std::vector<double> weights(static_cast<size_t>(n_nodes));
    const int m = (n_nodes + 1) / 2;
    for (int i = 0; i < m; ++i) {
        double z = std::cos(M_PI * (static_cast<double>(i) + 0.75)
                            / (static_cast<double>(n_nodes) + 0.5));
        for (int iter = 0; iter < 100; ++iter) {
            // Evaluate P_n(z) via three-term recurrence.
            double pjm1 = 1.0, pj = z;
            for (int j = 1; j < n_nodes; ++j) {
                const double pjp1 =
                    ((2.0*j + 1.0)*z*pj - static_cast<double>(j)*pjm1)
                    / static_cast<double>(j + 1);
                pjm1 = pj; pj = pjp1;
            }
            // pj = P_n(z), pjm1 = P_{n-1}(z)
            const double deriv = static_cast<double>(n_nodes)
                                 * (z * pj - pjm1) / (z * z - 1.0);
            const double dz = pj / deriv;
            z -= dz;
            if (std::abs(dz) < 1e-15) break;
        }
        // Recompute P_{n-1}(z) at the converged root for the weight.
        double pjm1 = 1.0, pj = z;
        for (int j = 1; j < n_nodes; ++j) {
            const double pjp1 = ((2.0*j + 1.0)*z*pj - static_cast<double>(j)*pjm1) / static_cast<double>(j + 1);
            pjm1 = pj; pj = pjp1;
        }
        const double deriv = static_cast<double>(n_nodes) * (z * pj - pjm1) / (z * z - 1.0);
        nodes[static_cast<size_t>(i)]               = -z;
        nodes[static_cast<size_t>(n_nodes - 1 - i)] =  z;

        const double w = 2.0 / ((1.0 - z * z) * deriv * deriv);
        weights[static_cast<size_t>(i)]               = w;
        weights[static_cast<size_t>(n_nodes - 1 - i)] = w;
    }

    // Accumulate the integrand at each node.
    // The integrand is the unnormalised product:
    //   p_M^(alpha_G-1) * (1-p_M)^(beta_G-1)
    //   * p_M^k_M * (1-p_M)^(d_M-k_M)
    //   * BetaBin(k_C | d_C, p_M*(Ne-1), (1-p_M)*(Ne-1))
    // divided by B(alpha_G, beta_G) to make it a proper Beta prior.
    // We work in log-space and factor out a scale for numerical stability.
    const int k_M = pd.m_ad_alt, d_M = pd.m_dp;
    const int k_C = pd.c_ad_alt, d_C = pd.c_dp;

    // Precompute the log-normalisation of the Beta prior.
    const double log_beta_norm = std::lgamma(alpha_G) + std::lgamma(beta_G) 
                               - std::lgamma(alpha_G + beta_G);

    // Evaluate log-integrand at each node and track max for log-sum-exp.
    std::vector<double> log_vals(static_cast<size_t>(n_nodes));
    double max_log = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < n_nodes; ++i) {
        const double t  = nodes[static_cast<size_t>(i)];
        const double pM = 0.5 * (1.0 + t);
        if (pM <= 0.0 || pM >= 1.0) {
            log_vals[static_cast<size_t>(i)] = -std::numeric_limits<double>::infinity();
            continue;
        }
        const double log_prior_pdf = (alpha_G - 1.0) * std::log(pM) 
                                   + (beta_G  - 1.0) * std::log(1.0 - pM)
                                   - log_beta_norm;
        const double log_binom_m = static_cast<double>(k_M) * std::log(pM)
                                 + static_cast<double>(d_M - k_M) * std::log(1.0 - pM);
        const double log_betabin_c = lf.log_betabinom_pmf(d_C, k_C, pM * ne1, (1.0 - pM) * ne1);
        // Jacobian dp_M/dt = 0.5; we add log(0.5) = -log(2) at the end.
        log_vals[static_cast<size_t>(i)] = log_prior_pdf + log_binom_m + log_betabin_c;
        
        if (log_vals[static_cast<size_t>(i)] > max_log) {
            max_log = log_vals[static_cast<size_t>(i)];
        }
    }

    // Weighted sum via log-sum-exp for stability.
    double sum_exp = 0.0;
    for (int i = 0; i < n_nodes; ++i) {
        if (!std::isfinite(log_vals[static_cast<size_t>(i)])) continue;
        sum_exp += weights[static_cast<size_t>(i)] * std::exp(log_vals[static_cast<size_t>(i)] - max_log);
    }
    if (sum_exp <= 0.0) return -std::numeric_limits<double>::infinity();

    // log C(d_M, k_M) (binomial coeff for the mother's read sampling)
    const double log_binom_m_coeff = lf.log_comb(d_M, k_M);

    return log_binom_m_coeff + max_log + std::log(sum_exp)
           - std::log(2.0);  // Jacobian dp_M/dt = 0.5
}

double NeEstimator::compute_global_ll(int ne,
                                      const std::vector<PairData>& data,
                                      const LogFactorial& lf,
                                      bool continuous) {
    double total = 0.0;
    for (const PairData& pd : data) {
        if (continuous) {
            total += (pd.has_g == 1)
                ? compute_ll_trio_continuous(pd, static_cast<double>(ne), lf)
                : compute_ll_single_continuous(pd, ne, lf);
        } else {
            total += compute_ll_single(pd, ne, lf);
        }
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
                if (continuous) {
                    s += (data[i].has_g == 1)
                        ? compute_ll_trio_continuous(data[i], static_cast<double>(ne), lf)
                        : compute_ll_single_continuous(data[i], ne, lf);
                } else {
                    s += compute_ll_single(data[i], ne, lf);
                }
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
    if (min_ne < 1)      min_ne = 1;
    if (max_ne < min_ne) max_ne = min_ne;

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

namespace {

// Body of NeEstimator::compute_global_ll_continuous(), split out so the
// scan-level callers (estimate_continuous and the --ne-profile grid walk) can
// hand in a maternal cache that outlives a single evaluation.  It stays in
// this translation unit because TrioQuadratureCache is internal and must not
// leak into the public header.
//
// `maternal_scan_cache` may be null, in which case a per-evaluation cache is
// used -- the behaviour the public entry point had before, and numerically
// identical either way since a rule is a pure function of its key.
double global_ll_continuous_impl(
    double ne, const std::vector<NeEstimator::PairData>& data,
    const NeEstimator::LogFactorial& lf, int threads,
    const NeEstimator::ModelOptions& opts,
    TrioQuadratureCache* maternal_scan_cache) 
{
    // Trio rules are rebuilt at every evaluation: their key moves with Ne.
    TrioQuadratureCache trio_cache;
    TrioQuadratureCache maternal_local(kMaternalScanCacheNodeBudget);
    TrioQuadratureCache* maternal_cache = (maternal_scan_cache != nullptr) ? maternal_scan_cache : &maternal_local;

    if (threads <= 1 || data.size() < 64) {
        double total = 0.0;
        for (const NeEstimator::PairData& pd : data) {
            total += compute_row_ll_continuous(pd, ne, lf, &trio_cache,
                                               maternal_cache, opts);
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
        futures.emplace_back(pool.submit([&, start, end, ne, opts]() {
            double s = 0.0;
            for (size_t i = start; i < end; ++i) {
                s += compute_row_ll_continuous(data[i], ne, lf, &trio_cache,
                                               maternal_cache, opts);
            }
            return s;
        }));
    }
    double total = 0.0;
    for (auto& f : futures) total += f.get();
    return total;
}
}  // namespace

double NeEstimator::compute_global_ll_continuous(
    double ne, const std::vector<PairData>& data,
    const LogFactorial& lf, int threads, ModelOptions opts) 
{
    return global_ll_continuous_impl(ne, data, lf, threads, opts, nullptr);
}

// Golden-section search to find the maximum of a unimodal function
// on [a, b] to within tolerance `tol`.  Returns the Ne that maximizes
// compute_global_ll_continuous.
//
// `maternal_cache` is the scan-lifetime Gauss-Jacobi rule cache owned by
// estimate_continuous(); every evaluation in this refinement reuses it.
static double golden_section_max(
    double a, double b, double tol,
    const std::vector<NeEstimator::PairData>& data,
    const NeEstimator::LogFactorial& lf, int threads,
    NeEstimator::ModelOptions opts,
    TrioQuadratureCache* maternal_cache) 
{
    constexpr double phi = 0.6180339887498949;  // 黄金分割点：(sqrt(5)-1)/2
    double x1 = b - phi * (b - a);
    double x2 = a + phi * (b - a);
    double f1 = global_ll_continuous_impl(x1, data, lf, threads, opts, maternal_cache);
    double f2 = global_ll_continuous_impl(x2, data, lf, threads, opts, maternal_cache);
    while ((b - a) > tol) {
        if (f1 < f2) {
            a  = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + phi * (b - a);
            f2 = global_ll_continuous_impl(x2, data, lf, threads, opts, maternal_cache);
        } else {
            b  = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - phi * (b - a);
            f1 = global_ll_continuous_impl(x1, data, lf, threads, opts, maternal_cache);
        }
    }
    return (a + b) * 0.5;
}

// Body of NeEstimator::find_optimal_ne_continuous(); split out for the same
// reason as global_ll_continuous_impl, so the caller's scan-lifetime cache can
// be threaded through both the coarse integer sweep and the golden-section
// refinement.
static double find_optimal_ne_continuous_impl(
    const std::vector<NeEstimator::PairData>& data,
    const NeEstimator::LogFactorial& lf,
    int min_ne, int max_ne, int threads,
    NeEstimator::ModelOptions opts,
    TrioQuadratureCache* maternal_cache)
{
    if (min_ne < 1) min_ne = 1;
    if (max_ne < min_ne) max_ne = min_ne;

    // Phase 1: coarse integer scan to bracket the peak.
    int    best_int = min_ne;
    double best_ll  = global_ll_continuous_impl(
        static_cast<double>(min_ne), data, lf, threads, opts, maternal_cache);
    for (int ne = min_ne + 1; ne <= max_ne; ++ne) {
        const double ll = global_ll_continuous_impl(
            static_cast<double>(ne), data, lf, threads, opts, maternal_cache);
        if (ll > best_ll) {
            best_ll  = ll;
            best_int = ne;
        }
    }

    // Golden-section search samples only interior points.  When the coarse
    // scan identifies a search boundary as the maximum, preserve that exact
    // constrained optimum rather than returning a nearby interior value.
    if (best_int == min_ne || best_int == max_ne) {
        return static_cast<double>(best_int);
    }

    // Phase 2: golden-section refinement in [best_int - 1, best_int + 1].
    const double lo = std::max(static_cast<double>(min_ne), static_cast<double>(best_int) - 1.0);
    const double hi = std::min(static_cast<double>(max_ne), static_cast<double>(best_int) + 1.0);
    return golden_section_max(lo, hi, 0.01, data, lf, threads, opts, maternal_cache);
}

double NeEstimator::find_optimal_ne_continuous(
    const std::vector<PairData>& data, 
    const LogFactorial& lf,
    int min_ne, int max_ne, int threads, ModelOptions opts) 
{
    return find_optimal_ne_continuous_impl(data, lf, min_ne, max_ne, 
                                           threads, opts, nullptr);
}

NeEstimator::Result NeEstimator::estimate_continuous(
    const std::vector<PairData>& data, 
    int min_ne, int max_ne, int threads, ModelOptions opts) 
{
    Result r;
    r.n_pairs = data.size();
    r.maternal_marginalization = opts.marginalize_maternal;
    if (data.empty()) {
        throw std::runtime_error("[ne-estimate] No transmission pairs to fit.");
    }
    if (min_ne < 1)        min_ne = 1;
    if (max_ne < min_ne)   max_ne = min_ne;

    const int cache_size = required_cache_size(data, max_ne);
    LogFactorial lf(cache_size);

    // Scan-lifetime Gauss-Jacobi rule cache.  The maternal-marginalised rule
    // key (k_M+1, d_M-k_M+1, order) carries no Ne, so the single cache built
    // here serves the coarse integer sweep, the golden-section refinement, the
    // maximum-LL re-evaluation and both profile-likelihood CI walks -- a few
    // hundred evaluations that would otherwise each rebuild every row's rule
    // from scratch.  Rule construction is the dominant per-row cost (97.5% at
    // d_M=d_C=2000), so this is what keeps the marginalised likelihood, which
    // is the default, affordable.  Bounded by kMaternalScanCacheNodeBudget;
    // trio rows do not use it because their key moves with Ne.
    TrioQuadratureCache maternal_cache(kMaternalScanCacheNodeBudget);

    r.ne = find_optimal_ne_continuous_impl(data, lf, min_ne, max_ne, threads, opts, &maternal_cache);
    r.max_log_lik = global_ll_continuous_impl(r.ne, data, lf, threads, opts, &maternal_cache);

    if (!std::isfinite(r.max_log_lik)) {
        throw std::runtime_error(
            "[ne-estimate] No finite continuous-model likelihood exists in the requested Ne range.");
    }

    const double thr = r.max_log_lik - kProfileLLThresholdDelta;
    constexpr double fine_step   = 0.01;
    constexpr double coarse_step = 0.10;

    auto ll_at = [&](double ne) {
        return global_ll_continuous_impl(ne, data, lf, threads, opts, &maternal_cache);
    };

    // Coarse-to-fine profile-likelihood CI walk.
    //
    // Phase 1 walks away from the MLE on the coarse grid (0.1) to bracket
    // the 1.92-drop crossing.  Phase 2 re-walks the bracket on the fine
    // grid (0.01) and stops at the first below-threshold point, exactly
    // like the legacy full fine walk.  The likelihood surface is smooth
    // around the MLE, so the coarse grid cannot skip the crossing in
    // practice; this yields the same boundary as the legacy walk (up to
    // float grid alignment) with ~10x fewer evaluations.
    auto walk_ci_side = [&](double start_ne, double direction, double bound, bool* clipped_out) {
        // Phase 1: coarse walk.  `cursor` stops at the first
        // below-threshold coarse point or exits the search range.
        double last_above = start_ne;
        double cursor     = start_ne + direction * coarse_step;
        bool   crossed    = false;
        while ((direction > 0.0) ? (cursor <= bound) : (cursor >= bound)) {
            if (ll_at(cursor) < thr) { crossed = true; break; }
            last_above = cursor;
            cursor += direction * coarse_step;
        }
        // The fine-walk endpoint: the coarse crossing point itself (already
        // known below threshold, the fine walk terminates on it) or the
        // search bound (which the legacy walk evaluates, so we do too).
        const double end = crossed ? cursor : bound;

        // Phase 2: fine walk over (last_above, end].
        double fine = last_above + direction * fine_step;
        while ((direction > 0.0) ? (fine <= end) : (fine >= end)) {
            if (ll_at(fine) < thr) {
                *clipped_out = false;
                return fine - direction * fine_step;
            }
            if (fine == end) break;   // endpoint handled below
            fine += direction * fine_step;
        }
        
        if (crossed) {
            // Every fine-grid point above the bracket held; the coarse
            // crossing point is the first below-threshold point.
            *clipped_out = false;
            return cursor - direction * fine_step;
        }
        // No crossing on either grid: CI runs to the search bound.
        *clipped_out = true;
        return bound;
    };

    r.ci_low  = walk_ci_side(r.ne, -1.0, static_cast<double>(min_ne), &r.ci_low_clipped);
    r.ci_high = walk_ci_side(r.ne, +1.0, static_cast<double>(max_ne), &r.ci_high_clipped);

    return r;
}

NeEstimator::Result NeEstimator::estimate(const std::vector<PairData>& data,
                                          int min_ne, int max_ne, int threads,
                                          bool continuous, ModelOptions opts) {
    // For the continuous model, delegate to the real-valued optimizer.
    if (continuous) {
        return estimate_continuous(data, min_ne, max_ne, threads, opts);
    }
    // The discrete model already integrates the maternal frequency over
    // Beta(m_alt+1, m_dp-m_alt+1), so ModelOptions has no meaning here and is
    // a silent no-op.  Reporting that is the CLI's job (_parse_args), because
    // marginalize_maternal is now true by default: warning here would fire on
    // every `--model discrete` run and describe the default configuration as
    // if the user had misconfigured it.

    const auto trio_row = std::find_if(data.begin(), data.end(),
        [](const PairData& pd) { return pd.has_g == 1; });
    
    if (trio_row != data.end()) {
        throw std::runtime_error(
            "[ne-estimate] --model discrete does not support HAS_G=1 trio rows; "
            "use --model continuous or remove trio rows.");
    }

    Result r;
    r.n_pairs = data.size();
    if (data.empty()) {
        throw std::runtime_error("[ne-estimate] No transmission pairs to fit.");
    }
    if (min_ne < 1)      min_ne = 1;
    if (max_ne < min_ne) max_ne = min_ne;

    const int cache_size = required_cache_size(data, max_ne);
    LogFactorial lf(cache_size);

    const int best_ne = find_optimal_ne(data, lf, min_ne, max_ne, threads, false);
    r.ne          = static_cast<double>(best_ne);
    r.max_log_lik = compute_global_ll_parallel(best_ne, data, lf, threads, false);
    if (!std::isfinite(r.max_log_lik)) {
        throw std::runtime_error(
            "[ne-estimate] No finite discrete-model likelihood exists in the requested Ne range.");
    }

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

// Locate column index by case-sensitive name; throws when missing.
static int col_index(const std::vector<std::string>& cols, const std::string& name) {
    for (size_t i = 0; i < cols.size(); ++i) {
        if (cols[i] == name) return static_cast<int>(i);
    }
    throw std::runtime_error("[ne-estimate] Required column missing in TSV: " + name);
}

std::vector<NeEstimator::PairData>
NeEstimator::load_pairs(const std::string& tsv_path, double min_vaf, double max_vaf,
                       NeEstimator::LoadStats* stats, int min_depth) {
    if (!std::isfinite(min_vaf) || !std::isfinite(max_vaf) ||
        min_vaf < 0.0 || max_vaf > 1.0 || min_vaf > max_vaf) {
        throw std::invalid_argument("[ne-estimate] Invalid finite maternal-VAF window.");
    }
    if (min_depth < 0) {
        throw std::invalid_argument("[ne-estimate] --min-depth must be >= 0 (0 disables the filter).");
    }
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

    // Optional trio columns (backward-compatible with legacy 16-col TSVs).
    // When HAS_G == 1 the row is a G-M-C trio; the grandmother's DP and
    // AD_ALT are read from GRANDMOTHER_DP / GRANDMOTHER_AD_ALT.
    auto opt_col = [&](const std::string& name) -> int {
        for (size_t i = 0; i < header_cols.size(); ++i)
            if (header_cols[i] == name) return static_cast<int>(i);
        return -1;
    };
    const int idx_has_g    = opt_col("HAS_G");
    const int idx_g_dp     = opt_col("GRANDMOTHER_DP");
    const int idx_g_ad_alt = opt_col("GRANDMOTHER_AD_ALT");
    const bool has_any_trio_column = idx_has_g >= 0 || idx_g_dp >= 0 || idx_g_ad_alt >= 0;
    const bool has_trio_columns = idx_has_g >= 0 && idx_g_dp >= 0 && idx_g_ad_alt >= 0;
    if (has_any_trio_column && !has_trio_columns) {
        throw std::runtime_error(
            "[ne-estimate] TSV must contain HAS_G, GRANDMOTHER_DP, and "
            "GRANDMOTHER_AD_ALT together.");
    }

    // Family identifier columns (backward-compatible; empty when absent).
    const int idx_fam_id    = opt_col("FAM_ID");
    const int idx_mother_id = opt_col("MOTHER_ID");
    const int idx_child_id  = opt_col("CHILD_ID");
    const int idx_chrom     = opt_col("CHROM");
    const int idx_pos       = opt_col("POS");
    const int idx_alt       = opt_col("ALT");
    const bool has_locus_columns = idx_chrom >= 0 && idx_pos >= 0 && idx_alt >= 0;

    const int max_required_index = std::max({
        idx_m_dp, idx_m_ad_alt, idx_m_vaf, idx_c_dp, idx_c_ad_alt, idx_qc,
        idx_has_g, idx_g_dp, idx_g_ad_alt,
        idx_fam_id, idx_mother_id, idx_child_id,
        has_locus_columns ? idx_chrom : -1,
        has_locus_columns ? idx_pos : -1,
        has_locus_columns ? idx_alt : -1
    });

    std::vector<PairData> data;
    while (std::getline(in_file, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        std::vector<std::string> tk = ngslib::split(line, "\t");
        if (static_cast<int>(tk.size()) <= max_required_index) continue;
        if (tk[idx_qc] != "PASS") continue;

        try {
            PairData pd;
            pd.m_dp     = std::stoi(tk[idx_m_dp]);
            pd.m_ad_alt = std::stoi(tk[idx_m_ad_alt]);
            pd.c_dp     = std::stoi(tk[idx_c_dp]);
            pd.c_ad_alt = std::stoi(tk[idx_c_ad_alt]);
            if (pd.m_dp <= 0 || pd.c_dp <= 0) continue;
            if (pd.m_ad_alt < 0 || pd.m_ad_alt > pd.m_dp) continue;
            if (pd.c_ad_alt < 0 || pd.c_ad_alt > pd.c_dp) continue;

            // Optional minimum-coverage gate, applied to BOTH sides.  A shallow
            // mother makes p_hat_M unreliable (the plug-in bias scales as
            // 1/d_M) and a shallow child makes k_C nearly uninformative about
            // the bottleneck.  This is a QC lever that complements, rather than
            // replaces, the default maternal marginalisation: it removes the
            // shallow rows instead of modelling them, so it also shrinks the
            // sample and moves the estimand with the retained depth
            // composition.  Under the default its remaining job is to gate out
            // the d_M ~ 30 regime where the Beta(1,1) prior mismatch is still
            // ~+14%.
            if (min_depth > 0 && (pd.m_dp < min_depth || pd.c_dp < min_depth)) {
                if (stats) ++stats->min_depth_skipped;
                continue;
            }

            // MOTHER_VAF remains a required provenance column in the TSV
            // schema, but AD/DP is the canonical value used everywhere in
            // the estimator. Filtering on it prevents stale or non-finite
            // text VAF fields from changing the fitted data set.
            const double m_vaf = static_cast<double>(pd.m_ad_alt)
                               / static_cast<double>(pd.m_dp);
            if (m_vaf < min_vaf || m_vaf > max_vaf) continue;

            // A declared G-M-C row must provide valid grandmother counts;
            // otherwise it is not a valid one-generation M-C observation.
            if (has_trio_columns) {
                const int hg = std::stoi(tk[idx_has_g]);
                if (hg != 0 && hg != 1) continue;
                if (hg == 1) {
                    const int gdp  = std::stoi(tk[idx_g_dp]);
                    const int galt = std::stoi(tk[idx_g_ad_alt]);
                    if (gdp <= 0 || galt < 0 || galt > gdp) continue;
                    pd.g_dp     = gdp;
                    pd.g_ad_alt = galt;
                    pd.has_g    = 1;

                    // A homoplasmic grandmother is an absorbing founder
                    // state (see compute_ll_trio_degenerate): compatible
                    // all-homoplasmic descendants contribute the
                    // Ne-independent constant log 1 = 0, and any
                    // segregating descendants are impossible (-inf) and
                    // would abort the whole fit.  Drop both at load time
                    // and count them for QC transparency.
                    const bool g_hom_ref = (galt == 0);
                    const bool g_hom_alt = (galt == gdp);
                    if (g_hom_ref || g_hom_alt) {
                        const bool desc_all_ref =
                            pd.m_ad_alt == 0 && pd.c_ad_alt == 0;
                        const bool desc_all_alt =
                            pd.m_ad_alt == pd.m_dp && pd.c_ad_alt == pd.c_dp;
                        if ((g_hom_ref && desc_all_ref) ||
                            (g_hom_alt && desc_all_alt)) {
                            if (stats) ++stats->trio_founder_hom_skipped;
                        } else {
                            if (stats) ++stats->trio_founder_mismatch_skipped;
                        }
                        continue;
                    }
                }
            }

            // Family identifiers (optional; used by per-family mode).
            if (idx_fam_id >= 0)    pd.fam_id    = tk[idx_fam_id];
            if (idx_mother_id >= 0) pd.mother_id = tk[idx_mother_id];
            if (idx_child_id >= 0)  pd.child_id  = tk[idx_child_id];
            if (has_locus_columns) {
                pd.locus_id = tk[idx_chrom] + "\x1f" + tk[idx_pos] + "\x1f" + tk[idx_alt];
            }

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
//     s_i  =  p_m_i (1 - p_m_i) / (m_dp_i - 1)  +  p_c_i (1 - p_c_i) / (c_dp_i - 1)
//                       ^ mother sampling           ^ child sampling
//                         variance, unbiasedly         variance, unbiasedly
//                         estimated: since             estimated: since
//                         E[p_hat(1-p_hat)] =          E[p_hat(1-p_hat)] =
//                         p(1-p)(d-1)/d, dividing      p(1-p)(d-1)/d, dividing
//                         by (d-1) is unbiased for     by (d-1) is unbiased for
//                         p(1-p)/d = Var(p_hat).       p(1-p)/d = Var(p_hat).
//                         Depth-1 rows are excluded    (variance not identifiable).
//     V    =  Sigma_i (d_i - s_i)  /  Sigma_i p_m_i (1 - p_m_i)
//     b    =  1 - V
//     Ne_kimura  =  1 / (1 - b)               (single generation, g = 1)
//
// We clip b to (eps, 1 - eps) to avoid divide-by-zero on degenerate
// cohorts.  When `n_bootstrap > 0` we also do a pair-level non-parametric
// bootstrap and return the 2.5/97.5 percentile CI for b and Ne_kimura.
// The value reported here is *only* a sanity cross-check against the
// primary continuous MMLE. The two can diverge on real data because:
//
//   * The MMLE uses the full child count likelihood, while the Kimura
//     diagnostic uses only sampling-corrected squared frequency shifts.
//   * The Wonnapinij b is variance-only, so a small number of high-drift
//     outliers (errors / NUMTs / mixed populations) can pull Ne_kimura
//     *downward*.
//
// A single-pair contribution helper, factored out so the bootstrap loop
// can reuse it cheaply.

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
        // p_hat(1-p_hat)/(d-1) is unbiased for p(1-p)/d under binomial
        // sampling. At depth one the plug-in variance is not identifiable,
        // so the row is excluded from Kimura aggregation.
        if (pd.m_dp <= 1 || pd.c_dp <= 1) { out.push_back(c); continue; }

        c.pm = static_cast<double>(pd.m_ad_alt) / static_cast<double>(pd.m_dp);
        c.pc = static_cast<double>(pd.c_ad_alt) / static_cast<double>(pd.c_dp);
        c.w  = c.pm * (1.0 - c.pm);
        if (c.w <= 0.0) { out.push_back(c); continue; }   // mother homoplasmic

        const double d = (c.pc - c.pm) * (c.pc - c.pm);
        const double s = c.pm * (1.0 - c.pm) / static_cast<double>(pd.m_dp - 1)
                       + c.pc * (1.0 - c.pc) / static_cast<double>(pd.c_dp - 1);
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
    if (V < 0.0) V = 0.0;       // sampling noise dominates
    
    double b = 1.0 - V;
    if (b < kEps) {
        b = kEps;
    } else if (b > 1.0 - kEps) {
        b = 1.0 - kEps;
    }
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
    
    if (lo == hi) return xs[lo];
    const double w = pos - static_cast<double>(lo);
    return xs[lo] * (1.0 - w) + xs[hi] * w;
}

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
    if (b == kEps && out.note.empty()) { 
        out.note = "b clipped to lower bound"; 
    } else if (b == 1.0 - kEps && out.note.empty()) { 
        out.note = "b clipped to upper bound"; 
    }

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
        const size_t n_drop = static_cast<size_t>(std::floor(trim_frac * static_cast<double>(info.size())));
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
// Per-bin observed drift summary with analytical Kimura predictions
// ---------------------------------------------------------------------
//
// For each equal-width maternal-VAF bin in [vaf_low, vaf_high] we report
// three independent estimates of the per-bin drift, all computed on the
// same informative pairs that drive the MMLE / Wonnapinij b:
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
// observed bin means using the fitted Ne and its 95% CI.
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
// Ne-profile scan (compare MMLE vs Kimura preferences across Ne)
// ---------------------------------------------------------------------
//
// For each candidate Ne in [min_ne, max_ne] (step `step`) compute two
// independent goodness-of-fit metrics on the *same* informative pair set:
//
//   mmle_log_lik(Ne) =  global marginal log-likelihood under the configured
//                       model (continuous Beta-diffusion or discrete Beta-
//                       Binomial).  Maximised at the fitted Ne_MMLE.
//   kimura_ssr(Ne)   =  Sigma_i ( (d_i - s_i) - p_m_i (1 - p_m_i) / Ne )^2
//                       Per-pair least-squares fit of the one-generation
//                       Wright-Fisher prediction.  Minimised at the
//                       analytic best Ne_kimura_ssr = Sigma w^2 / Sigma rw.
//
// The Kimura SSR is a closed-form least-squares analogue of the
// Kimura distribution-fit approach: it minimises the per-pair quadratic
// residual between observed drift r_i and the Wright-Fisher prediction
// w_i / Ne.  We report both profiles so the user can see whether the
// two estimators agree on the location of the best Ne.
std::vector<NeEstimator::NeProfileRow>
NeEstimator::compute_ne_profile(const std::vector<PairData>& data,
                                const LogFactorial& lf,
                                double min_ne, double max_ne, double step,
                                int threads, bool continuous,
                                ModelOptions opts) {
    if (!std::isfinite(min_ne) || !std::isfinite(max_ne) || !std::isfinite(step)) {
        throw std::invalid_argument("[ne-estimate] Ne-profile bounds and step must be finite.");
    }
    if (step <= 0.0)         step = 0.1;
    if (min_ne < 1.0)        min_ne = 1.0;
    if (max_ne < min_ne)     max_ne = min_ne;

    // Pre-compute (w_i, r_i) for every informative pair so the inner
    // Kimura-SSR loop is O(n) regardless of how dense the Ne grid is.
    const std::vector<PairContribution> contribs = prepare_pair_contributions(data);

    std::vector<NeProfileRow> profile;
    size_t continuous_steps = 0;
    if (continuous) {
        const double step_count = std::floor((max_ne - min_ne) / step);
        if (step_count > static_cast<double>(std::numeric_limits<size_t>::max() - 2)) {
            throw std::invalid_argument("[ne-estimate] Ne-profile grid is too large.");
        }
        continuous_steps = static_cast<size_t>(step_count);
        profile.reserve(continuous_steps + 2);
    } else {
        const int first_ne = std::max(1, static_cast<int>(std::ceil(min_ne)));
        const int last_ne = static_cast<int>(std::floor(max_ne));
        if (last_ne >= first_ne) {
            profile.reserve(static_cast<size_t>(last_ne - first_ne + 1));
        }
    }

    double max_ll  = -std::numeric_limits<double>::infinity();
    double min_ssr =  std::numeric_limits<double>::infinity();

    // The profile is itself a scan over Ne (up to (max_ne - min_ne) / step
    // evaluations, thousands at the default step of 0.1), so it gets the same
    // scan-lifetime maternal rule cache the optimizer uses.  Captured by
    // reference into append_profile_row below.
    TrioQuadratureCache profile_maternal_cache(kMaternalScanCacheNodeBudget);

    auto append_profile_row = [&](double ne, int discrete_ne) {
        NeProfileRow row;
        row.ne_candidate = ne;

        if (continuous) {
            row.mmle_log_lik = global_ll_continuous_impl(ne, data, lf, threads, opts, &profile_maternal_cache);
        } else {
            row.mmle_log_lik = compute_global_ll_parallel(discrete_ne, data, lf, threads, false);
        }
        if (row.mmle_log_lik > max_ll) max_ll = row.mmle_log_lik;

        // Kimura SSR: Sigma_i (r_i - w_i / Ne)^2 over informative pairs.
        const double inv_ne = 1.0 / ne;
        double ssr = 0.0;
        for (const auto& c : contribs) {
            if (!c.informative) continue;
            const double resid = c.r - c.w * inv_ne;
            ssr += resid * resid;
        }
        row.kimura_ssr = ssr;
        if (ssr < min_ssr) min_ssr = ssr;

        profile.push_back(row);
    };

    if (continuous) {
        const double tolerance = 8.0 * std::numeric_limits<double>::epsilon()
                               * std::max({1.0, std::abs(min_ne), std::abs(max_ne)});
        for (size_t index = 0; index <= continuous_steps; ++index) {
            double ne = min_ne + step * static_cast<double>(index);
            if (ne > max_ne) {
                if (ne - max_ne > tolerance) break;
                ne = max_ne;
            }
            append_profile_row(ne, 0);
        }
        if (profile.empty() || max_ne - profile.back().ne_candidate > tolerance) {
            append_profile_row(max_ne, 0);
        }
    } else {
        const int first_ne = std::max(1, static_cast<int>(std::ceil(min_ne)));
        const int last_ne = static_cast<int>(std::floor(max_ne));
        for (int ne = first_ne; ne <= last_ne; ++ne) {
            append_profile_row(static_cast<double>(ne), ne);
        }
    }

    // Post-process: normalise both metrics so they are on comparable
    // "distance from best fit" scales for plotting.
    for (auto& row : profile) {
        row.mmle_delta_2ll  = -2.0 * (row.mmle_log_lik - max_ll);
        row.kimura_norm_ssr = (min_ssr > 0.0) ? row.kimura_ssr / min_ssr : 1.0;
    }
    return profile;
}

double NeEstimator::kimura_ssr_best_ne(const std::vector<PairData>& data) {
    // Closed-form least-squares estimator for Ne under the prediction
    //     E[d_i - s_i]  =  p_m_i (1 - p_m_i) / Ne.
    // Setting d/dNe Sigma (r_i - w_i / Ne)^2  =  0 and solving:
    //     Ne_best  =  Sigma w_i^2  /  Sigma r_i w_i.
    const std::vector<PairContribution> contribs = prepare_pair_contributions(data);
    double sum_w_sq    = 0.0;
    double sum_r_w     = 0.0;
    for (const auto& c : contribs) {
        if (!c.informative) continue;
        sum_w_sq += c.w * c.w;
        sum_r_w  += c.r * c.w;
    }
    if (!(sum_r_w > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    return sum_w_sq / sum_r_w;
}

// ---------------------------------------------------------------------
// Per-family estimation
// ---------------------------------------------------------------------

std::vector<NeEstimator::FamilyData>
NeEstimator::group_into_families(const std::vector<PairData>& data) {
    // Key = (fam_id, mother_id) -> index into result vector.
    std::map<std::pair<std::string,std::string>, size_t> key_to_idx;
    std::vector<FamilyData> families;

    for (const auto& pd : data) {
        // Legacy TSVs without FAM_ID: group everything into one family.
        const std::string fid = pd.fam_id.empty() ? "ALL" : pd.fam_id;
        const std::string mid = pd.mother_id.empty() ? "ALL" : pd.mother_id;
        auto key = std::make_pair(fid, mid);

        auto it = key_to_idx.find(key);
        if (it == key_to_idx.end()) {
            FamilyData fd;
            fd.fam_id    = fid;
            fd.mother_id = mid;
            fd.pairs.push_back(pd);
            if (!pd.child_id.empty()) fd.child_ids.push_back(pd.child_id);
            key_to_idx[key] = families.size();
            families.push_back(std::move(fd));
        } else {
            auto& fam = families[it->second];
            fam.pairs.push_back(pd);
            // Track unique child IDs.
            if (!pd.child_id.empty()) {
                bool found = false;
                for (const auto& cid : fam.child_ids) {
                    if (cid == pd.child_id) { found = true; break; }
                }
                if (!found) fam.child_ids.push_back(pd.child_id);
            }
        }
    }
    return families;
}

NeEstimator::FamilyResult
NeEstimator::estimate_family(const FamilyData& fam,
                              int min_ne, int max_ne,
                              int min_family_sites, ModelOptions opts) {
    FamilyResult fr;
    fr.fam_id    = fam.fam_id;
    fr.mother_id = fam.mother_id;
    fr.n_children = fam.child_ids.size();
    fr.n_pairs   = fam.pairs.size();

    // Count informative sites (mother heteroplasmic: 0 < p_M < 1).
    // Note this count is a property of the plug-in view of the data: under the
    // default marginalisation a homoplasmic mother (k_M = 0 or k_M = d_M)
    // still leaves a proper Beta posterior over p_M and does carry information,
    // so min_family_sites stays a deliberately conservative gate either way.
    size_t n_info = 0;
    double sum_m_dp = 0.0, sum_c_dp = 0.0;
    double sum_inv_m_dp = 0.0, sum_inv_c_dp = 0.0;
    size_t n_depth = 0;
    for (const auto& pd : fam.pairs) {
        if (pd.m_dp <= 0) continue;
        const double pm = static_cast<double>(pd.m_ad_alt)
                        / static_cast<double>(pd.m_dp);
        if (pm > 0.0 && pm < 1.0) ++n_info;
        sum_m_dp += pd.m_dp;
        sum_c_dp += pd.c_dp;
        sum_inv_m_dp += 1.0 / static_cast<double>(pd.m_dp);
        if (pd.c_dp > 0) sum_inv_c_dp += 1.0 / static_cast<double>(pd.c_dp);
        ++n_depth;
    }
    fr.n_informative = n_info;
    fr.mean_mother_dp = fam.pairs.empty() ? 0.0 : sum_m_dp / fam.pairs.size();
    fr.mean_child_dp  = fam.pairs.empty() ? 0.0 : sum_c_dp / fam.pairs.size();
    fr.harmonic_mother_dp = (n_depth == 0 || !(sum_inv_m_dp > 0.0))
                          ? 0.0 : static_cast<double>(n_depth) / sum_inv_m_dp;
    fr.harmonic_child_dp  = (n_depth == 0 || !(sum_inv_c_dp > 0.0))
                          ? 0.0 : static_cast<double>(n_depth) / sum_inv_c_dp;

    if (static_cast<int>(n_info) < min_family_sites) {
        fr.skipped = true;
        fr.warning = "too few informative sites (" + std::to_string(n_info)
                   + " < " + std::to_string(min_family_sites) + ")";
        return fr;
    }

    // Reuse the existing continuous MMLE estimator on this family's pairs.
    Result est = estimate_continuous(fam.pairs, min_ne, max_ne, /*threads=*/1, opts);
    fr.ne             = est.ne;
    fr.ci_low         = est.ci_low;
    fr.ci_high        = est.ci_high;
    fr.max_log_lik    = est.max_log_lik;
    fr.ci_low_clipped  = est.ci_low_clipped;
    fr.ci_high_clipped = est.ci_high_clipped;

    // Warnings for edge cases.
    if (n_info < 10) {
        fr.warning = "small sample (" + std::to_string(n_info) + " informative sites)";
    }
    if (fr.ci_low_clipped || fr.ci_high_clipped) {
        if (!fr.warning.empty()) fr.warning += "; ";
        fr.warning += "CI clipped at search boundary";
    }

    return fr;
}

std::vector<NeEstimator::FamilyResult>
NeEstimator::estimate_all_families(const std::vector<FamilyData>& families,
                                    int min_ne, int max_ne,
                                    int min_family_sites,
                                    int threads, ModelOptions opts) {
    std::vector<FamilyResult> results(families.size());

    if (threads <= 1 || families.size() <= 1) {
        // Serial path.
        for (size_t i = 0; i < families.size(); ++i) {
            results[i] = estimate_family(families[i], min_ne, max_ne,
                                          min_family_sites, opts);
        }
        return results;
    }

    // Embarrassingly parallel: distribute families across workers.
    const size_t n = families.size();
    auto worker = [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            results[i] = estimate_family(families[i], min_ne, max_ne,
                                          min_family_sites, opts);
        }
    };

    std::vector<std::future<void>> futures;
    const size_t chunk = (n + threads - 1) / threads;
    for (size_t t = 0; t < static_cast<size_t>(threads); ++t) {
        const size_t lo = t * chunk;
        const size_t hi = std::min(lo + chunk, n);
        if (lo >= n) break;
        futures.push_back(std::async(std::launch::async, worker, lo, hi));
    }
    for (auto& f : futures) f.get();

    return results;
}

NeEstimator::KimuraCheck
NeEstimator::compute_family_kimura_check(const FamilyData& fam,
                                          int n_bootstrap, uint64_t seed,
                                          double trim_frac) {
    // The point estimate and trimmed diagnostic remain row-summed. Bootstrap
    // uncertainty, however, resamples one locus with all child observations
    // together when a locus ID is available.
    KimuraCheck out = compute_kimura_check(fam.pairs, 0, seed, trim_frac, 0);
    if (n_bootstrap <= 0 || out.n_informative == 0) return out;

    const std::vector<PairContribution> contribs = prepare_pair_contributions(fam.pairs);
    std::map<std::string, std::vector<const PairContribution*>> clusters;
    bool used_row_fallback = false;
    for (const auto& contribution : contribs) {
        if (!contribution.informative) continue;
        const PairData& pair = fam.pairs[contribution.idx];
        std::string key = pair.locus_id;
        if (key.empty()) {
            key = "__row__" + std::to_string(contribution.idx);
            used_row_fallback = true;
        }
        clusters[key].push_back(&contribution);
    }
    if (clusters.empty()) return out;

    out.ci_computed = true;
    out.n_bootstrap = n_bootstrap;
    out.bootstrap_seed = seed;
    std::vector<const std::vector<const PairContribution*>*> cluster_list;
    cluster_list.reserve(clusters.size());
    for (const auto& cluster : clusters) cluster_list.push_back(&cluster.second);

    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<size_t> pick(0, cluster_list.size() - 1);
    std::vector<double> b_samples;
    b_samples.reserve(static_cast<size_t>(n_bootstrap));
    for (int replicate = 0; replicate < n_bootstrap; ++replicate) {
        double numerator = 0.0;
        double denominator = 0.0;
        for (size_t draw = 0; draw < cluster_list.size(); ++draw) {
            const auto& cluster = *cluster_list[pick(rng)];
            for (const PairContribution* contribution : cluster) {
                numerator += contribution->r;
                denominator += contribution->w;
            }
        }
        const double b = b_from_aggregates(numerator, denominator);
        if (!std::isnan(b)) b_samples.push_back(b);
    }

    if (b_samples.size() < 2) {
        out.b_ci_low = out.b;
        out.b_ci_high = out.b;
        out.ne_kimura_ci_low = out.ne_kimura;
        out.ne_kimura_ci_high = out.ne_kimura;
        if (!out.note.empty()) out.note += "; ";
        out.note += "family bootstrap CI degenerate (too few informative resamples)";
    } else {
        std::sort(b_samples.begin(), b_samples.end());
        const double b_low = quantile(b_samples, 0.025);
        const double b_high = quantile(b_samples, 0.975);
        out.b_ci_low = b_low;
        out.b_ci_high = b_high;
        out.ne_kimura_ci_low = 1.0 / (1.0 - b_low);
        out.ne_kimura_ci_high = 1.0 / (1.0 - b_high);
    }
    if (used_row_fallback) {
        if (!out.note.empty()) out.note += "; ";
        out.note += "family bootstrap used row clusters where locus identifiers were absent";
    }
    return out;
}

void NeEstimator::_write_family_tsv(
    const std::vector<FamilyResult>& results, std::ostream& out) const
{
    // Header.
    out << "FAM_ID\tMOTHER_ID\tN_CHILDREN\tN_PAIRS\tN_INFORMATIVE"
        << "\tNE_MMLE\tCI_95_LOW\tCI_95_HIGH"
        << "\tCI_LOW_CLIPPED\tCI_HIGH_CLIPPED"
        << "\tMAX_LOG_LIK"
        << "\tMEAN_MOTHER_DP\tMEAN_CHILD_DP";
    // Kimura columns (if any family has kimura computed).
    bool any_kimura = false;
    for (const auto& fr : results) {
        if (fr.kimura.computed) { any_kimura = true; break; }
    }
    if (any_kimura) {
        out << "\tKIMURA_b\tKIMURA_NE\tKIMURA_N_INFORMATIVE";
    }
    bool any_kimura_ci = false;
    for (const auto& fr : results) {
        if (fr.kimura.ci_computed) { any_kimura_ci = true; break; }
    }
    if (any_kimura_ci) {
        out << "\tKIMURA_B_CI_95_LOW\tKIMURA_B_CI_95_HIGH"
            << "\tKIMURA_NE_CI_95_LOW\tKIMURA_NE_CI_95_HIGH"
            << "\tKIMURA_N_BOOTSTRAP\tKIMURA_BOOTSTRAP_SEED";
    }
    out << "\tSKIPPED\tWARNING\n";

    // Rows.
    out << std::setprecision(4);
    for (const auto& fr : results) {
        out << fr.fam_id       << "\t"
            << fr.mother_id    << "\t"
            << fr.n_children   << "\t"
            << fr.n_pairs      << "\t"
            << fr.n_informative << "\t";
        if (fr.skipped) {
            out << "NA\tNA\tNA\tNA\tNA\tNA\t"
                << std::setprecision(8) << fr.mean_mother_dp << "\t"
                << fr.mean_child_dp;
        } else {
            out << std::setprecision(4)
                << fr.ne       << "\t"
                << fr.ci_low   << "\t"
                << fr.ci_high  << "\t"
                << (fr.ci_low_clipped  ? "TRUE" : "FALSE") << "\t"
                << (fr.ci_high_clipped ? "TRUE" : "FALSE") << "\t"
                << std::setprecision(8) << fr.max_log_lik << "\t"
                << std::setprecision(8) << fr.mean_mother_dp << "\t"
                << fr.mean_child_dp;
        }
        if (any_kimura) {
            if (fr.kimura.computed) {
                out << "\t" << std::setprecision(8) << fr.kimura.b
                    << "\t" << fr.kimura.ne_kimura
                    << "\t" << fr.kimura.n_informative;
            } else {
                out << "\tNA\tNA\tNA";
            }
        }
        if (any_kimura_ci) {
            if (fr.kimura.ci_computed) {
                out << "\t" << std::setprecision(8) << fr.kimura.b_ci_low
                    << "\t" << fr.kimura.b_ci_high
                    << "\t" << fr.kimura.ne_kimura_ci_low
                    << "\t" << fr.kimura.ne_kimura_ci_high
                    << "\t" << fr.kimura.n_bootstrap
                    << "\t" << fr.kimura.bootstrap_seed;
            } else {
                out << "\tNA\tNA\tNA\tNA\tNA\tNA";
            }
        }
        out << "\t" << (fr.skipped ? "TRUE" : "FALSE")
            << "\t" << fr.warning << "\n";
    }
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
    if (!std::isfinite(_config.min_vaf) || !std::isfinite(_config.max_vaf) ||
        _config.min_vaf < 0.0 || _config.max_vaf > 1.0 ||
        _config.min_vaf > _config.max_vaf) {
        throw std::invalid_argument("[ne-estimate] Invalid finite VAF window.");
    }
    if (!std::isfinite(_config.kimura_trim) ||
        !std::isfinite(_config.ne_profile_step)) {
        throw std::invalid_argument("[ne-estimate] Floating-point options must be finite.");
    }
    if (_config.min_ne < 1)                                      _config.min_ne = 1;
    if (_config.max_ne < _config.min_ne)                         _config.max_ne = _config.min_ne;
    if (_config.threads < 1)                                     _config.threads = 1;
    if (_config.kimura_bootstrap < 0)                            _config.kimura_bootstrap = 0;
    if (_config.kimura_trim < 0.0 || _config.kimura_trim >= 1.0) _config.kimura_trim = 0.0;
    if (_config.top_drift_k < 0)                                 _config.top_drift_k = 0;
    if (_config.bin_simulation_n_bins < 1)                       _config.bin_simulation_n_bins = 1;
    if (_config.ne_profile_step <= 0.0)                          _config.ne_profile_step = 0.1;
    if (_config.model.empty())                                   _config.model = "continuous";
    if (_config.model != "continuous" && _config.model != "discrete") {
        throw std::invalid_argument("[ne-estimate] model must be continuous or discrete.");
    }
    if (_config.per_family && _config.model != "continuous") {
        throw std::invalid_argument(
            "[ne-estimate] per-family estimation currently requires model=continuous.");
    }
    if (_config.min_family_sites < 1) _config.min_family_sites = 1;
    if (_config.min_depth < 0) {
        throw std::invalid_argument(
            "[ne-estimate] min_depth must be >= 0 (0 disables the filter).");
    }
    // kimura_check defaults to false; callers may flip it explicitly.
}

void NeEstimator::usage() {
    std::cerr << "Usage: mitoquest ne-estimate [options] -i <pairs.tsv>\n\n"
                 "Description:\n"
                 "  Estimate the mitochondrial DNA bottleneck size (Ne) from mother-child\n"
                 "  transmission pairs using the Maximum Marginal Likelihood Estimator (MMLE).\n"
                 "  The child's latent true allele frequency is analytically integrated out\n"
                 "  (Beta-Binomial conjugacy) and pairs are treated as independent (composite\n"
                 "  likelihood), yielding a per-pair marginal log-likelihood that is maximised\n"
                 "  jointly over Ne.\n"
                 "\n"
                 "  The maternal frequency p_m is integrated out against its Beta\n"
                 "  read-sampling posterior BY DEFAULT under the continuous model.\n"
                 "  --no-maternal-marginalization falls back to the legacy plug-in\n"
                 "  point estimate p_m = m_alt / m_dp, which folds the mother's own\n"
                 "  read noise into drift and makes the fitted Ne depth-dependent\n"
                 "  (see Notes).  The legacy --model discrete path always integrates\n"
                 "  a maternal Beta posterior.\n"
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
                 "  Reports the MMLE point estimate of Ne and its 95% profile-likelihood\n"
                 "  confidence interval (LL_max - 1.92).\n"
                 "\nRequired options:\n"
                 "  -i, --input     FILE   Input transmission pairs TSV produced by\n"
                 "                         `mitoquest trans-prep`.\n"
                 "\nOptional options:\n"
                 "  -o, --output    FILE   JSON output file (default: stdout).\n"
                 "      --model     NAME   Marginal-likelihood model: `continuous` (default,\n"
                 "                         recommended for mtDNA) or `discrete`.\n"
                 "      --no-maternal-marginalization\n"
                 "                         Fall back to the legacy plug-in p_m = m_alt /\n"
                 "                         m_dp for the continuous model's two-generation\n"
                 "                         rows.  By default p_m is integrated out against\n"
                 "                         its Beta(m_alt+1, m_ref+1) read-sampling\n"
                 "                         posterior by adaptive Gauss-Jacobi quadrature;\n"
                 "                         the plug-in is faster but biased low by\n"
                 "                         ~d_M/(Ne+d_M) and is not comparable across\n"
                 "                         coverage, so it is available only by explicit\n"
                 "                         request.  Continuous model only: --model discrete\n"
                 "                         always integrates p_m and cannot honour this flag\n"
                 "                         (warned on stderr); trio rows also already\n"
                 "                         integrate p_m given the grandmother, so it is a\n"
                 "                         no-op for them.\n"
                 "      --min-vaf   FLOAT  Lower maternal VAF gate, inclusive [0.10].\n"
                 "      --max-vaf   FLOAT  Upper maternal VAF gate, inclusive [0.90].\n"
                 "      --min-depth INT    Drop pairs whose MOTHER or CHILD depth is below\n"
                 "                         this value (both sides are gated: a shallow child\n"
                 "                         makes k_C uninformative just as a shallow mother\n"
                 "                         makes p_m unreliable).  0 disables [0].\n"
                 "      --min-ne    INT    Smallest Ne value to consider [1].\n"
                 "      --max-ne    INT    Largest Ne value to consider  [100].\n"
                 "  -t, --threads   INT    Worker threads for the inner sum [1].\n"
                 "      --cross-check NAME   Optional secondary estimator alongside the\n"
                 "                           MMLE. Supported value: `kimura`, which\n"
                 "                           computes the Wonnapinij b and the implied\n"
                 "                           Ne (single-generation approximation).\n"
                 "      --kimura-bootstrap  INT  Non-parametric bootstrap iterations for the\n"
                 "                              Kimura cross-check 95% CI.  0 disables [1000].\n"
                 "      --kimura-seed      INT  RNG seed for the Kimura bootstrap [42].\n"
                 "      --kimura-trim    FLOAT  Fraction of highest-drift pairs to drop before\n"
                 "                              recomputing the Kimura b (robust trimmed mean,\n"
                 "                              0.0 disables, recommended 0.10) [0.0].\n"
                 "      --top-drift-k     INT   Emit the top-K highest-drift pairs in the JSON\n"
                 "                              output for outlier inspection (NUMTs / errors).\n"
                 "                              0 disables [0].\n"
                 "      --bin-simulation FILE   Emit a per-bin drift summary TSV: per\n"
                 "                              maternal-VAF bin observed mean drift vs\n"
                 "                              analytical Kimura prediction p_m(1 - p_m) / Ne\n"
                 "                              under the fitted Ne (and its 95% CI).  Bin\n"
                 "                              range = [--min-vaf, --max-vaf].\n"
                 "      --bin-simulation-bins INT  Number of equal-width maternal-VAF bins\n"
                 "                                 for --bin-simulation [10].\n"
                 "      --ne-profile     FILE   Emit a TSV that scores every candidate Ne\n"
                 "                              under both the MMLE marginal log-likelihood\n"
                 "                              and the Kimura per-pair SSR metric.  Useful to\n"
                 "                              visually compare which Ne each estimator\n"
                 "                              prefers (dual-objective Ne scan).\n"
                 "                              Grid range = [--min-ne, --max-ne].\n"
                 "      --ne-profile-step FLOAT Grid step on the Ne axis for --ne-profile\n"
                 "                              [0.1].\n"
                 "      --per-family              Enable per-family Ne estimation.\n"
                 "                                Groups pairs by FAM_ID + MOTHER_ID and\n"
                 "                                estimates Ne independently for each family.\n"
                 "      --min-family-sites INT    Minimum informative sites per family [3].\n"
                 "                                Families with fewer sites are skipped.\n"
                 "      --per-family-output FILE  Write per-family Ne results as TSV (one\n"
                 "                                row per family).  If not set, per-family\n"
                 "                                results are embedded in the JSON output.\n"
                 "  -h, --help               Print this help message.\n\n"
                 "Notes:\n"
                 "  * Sites with maternal VAF near 0 or 1 carry virtually no information\n"
                 "    about Ne, hence the default 0.10 - 0.90 window.\n"
                 "  * Confidence-interval bounds are flagged with `_clipped = true` in the\n"
                 "    JSON output when they hit the search boundary.\n"
                 "  * The plug-in maternal VAF (--no-maternal-marginalization) makes the\n"
                 "    fitted value a pseudo-true Ne_app obeying 1/Ne_app = 1/Ne + 1/d_M at a\n"
                 "    SINGLE maternal depth d_M, so plug-in estimates from cohorts sequenced\n"
                 "    at different coverage are not directly comparable.  That identity is not\n"
                 "    invertible on a mixed-depth cohort: no arithmetic, harmonic or\n"
                 "    depth-weighted summary of the reported Mean_Mother_DP /\n"
                 "    Harmonic_Mother_DP reproduces the plug-in fit, so those fields are\n"
                 "    descriptive only.  The default marginalisation removes the bias instead\n"
                 "    of trying to back-correct the reported Ne.\n"
                 "  * Cross-depth comparability of the default, measured at true Ne=30 with\n"
                 "    d_M swept over 30..5000 (240 pairs, d_C=2000, 20 replicate cohorts per\n"
                 "    depth): the plug-in fit moves over [16.5, 30.3], a spread of 46% of the\n"
                 "    true Ne, while the marginalised fit moves over [29.4, 34.2], 16% -- 2.9x\n"
                 "    tighter.  Over d_M in [60, 5000] alone the marginalised fit spans\n"
                 "    [29.4, 30.7], i.e. 4.5%, so at d_M >= 60 the reported Ne really is\n"
                 "    comparable across coverage.\n"
                 "  * Per-family fits need SITES, not depth.  At true Ne=30 with d_C=2000,\n"
                 "    d_M=200 and 8 replicates of nested truths (every site count shares one\n"
                 "    bottleneck realisation), growing a family from 12 to 48 to 192\n"
                 "    segregating sites takes the marginalised median bias +11.2% -> +2.9% ->\n"
                 "    +1.5% and the share of families pinned at --max-ne 9% -> 0% -> 0%,\n"
                 "    converging on the cohort-level +1.7..2.4%.  Over the same sweep the\n"
                 "    plug-in moves -4.7% -> -10.6% -> -11.3%, i.e. towards its own\n"
                 "    pseudo-true -13.0% at d_M=200: its apparent accuracy on 12 sites is a\n"
                 "    small-sample MLE right-skew cancelling the depth bias, and the\n"
                 "    cancellation breaks as soon as either factor changes.  Across depth the\n"
                 "    marginalised per-family median is the comparable quantity -- +1.5% at\n"
                 "    d_M=200 against +1.1% at d_M=2000 on 192 sites, 0.4 points apart,\n"
                 "    where the plug-in is 11.0 points apart.\n"
                 "  * Maternal marginalisation does not repair a mis-specified bottleneck.\n"
                 "    If the true bottleneck draws an integer number of genomes, the\n"
                 "    continuous Beta-diffusion transition is the wrong model form and both\n"
                 "    the plug-in and the marginalised fit stay biased; use --model discrete\n"
                 "    in that regime.  It also assumes a Beta(1,1) prior on the maternal true\n"
                 "    frequency (matching --model discrete).  When the cohort's true maternal\n"
                 "    frequencies are truncated away from 0 and 1 that prior is over-dispersed\n"
                 "    and the fitted Ne comes out HIGH, by an amount that decays with maternal\n"
                 "    depth: measured +14% at d_M=30, <= +5% for d_M >= 60 and <= 2.5% (about\n"
                 "    one standard error) for d_M >= 500.  It is not exactly unbiased, so\n"
                 "    prefer --min-depth over interpreting very shallow maternal rows at face\n"
                 "    value.\n"
                 "  * The VAF window is not the source of that residual.  Whether a row passes\n"
                 "    k_M/d_M in [--min-vaf, --max-vaf] depends on p_m and the prior but not\n"
                 "    on Ne, so the window contributes an Ne-independent factor to the\n"
                 "    likelihood and cannot shift the MLE.  Confirmed with p_m drawn right at\n"
                 "    the window edge (U[0.10,0.30], only 90% of rows passing at d_M=30):\n"
                 "    +13.0% bias with the default window against +11.6% with --min-vaf 0\n"
                 "    --max-vaf 1, a difference inside the standard error at every depth.\n"
                 "  * The pair log-likelihoods are summed as a composite likelihood: pairs are\n"
                 "    treated as independent even when several come from the same mother.  The\n"
                 "    95% profile-likelihood interval therefore applies the chi2_1 threshold\n"
                 "    1.92 without a sandwich correction and can be too narrow for cohorts\n"
                 "    with many siblings per mother; use --per-family for family-resolved\n"
                 "    inference.\n"
                 "  * The Kimura cross-check is approximate; the default continuous MMLE uses\n"
                 "    a plug-in maternal VAF and a Beta-Binomial child-count working likelihood\n"
                 "    and is the reported primary estimate. On real data the two can diverge\n"
                 "    when many concordant heteroplasmic pairs co-exist with a few high-\n"
                 "    drift outliers (errors / NUMTs / mixed populations).  Use\n"
                 "    --kimura-trim 0.10 (drops the top 10% of high-drift pairs) and/or\n"
                 "    --top-drift-k 20 (lists the worst pairs) to diagnose this.\n\n"
              << "Version: " << MITOQUEST_VERSION << "\n"
              << std::endl;
}

void NeEstimator::_parse_args(int argc, char* argv[]) {
    _config.input_tsv.clear();
    _config.output_file.clear();
    _config.min_vaf          = 0.10;
    _config.max_vaf          = 0.90;
    _config.min_ne           = 1;
    _config.max_ne           = 100;
    _config.threads          = 1;
    _config.model            = "continuous";
    _config.kimura_check     = false;
    _config.kimura_bootstrap = 1000;
    _config.kimura_seed      = 42;
    _config.kimura_trim      = 0.0;
    _config.top_drift_k      = 0;
    _config.bin_simulation_file.clear();
    _config.bin_simulation_n_bins = 10;
    _config.ne_profile_file.clear();
    _config.ne_profile_step       = 0.1;
    _config.per_family            = false;
    _config.min_family_sites      = 3;
    _config.per_family_output_file.clear();
    _config.maternal_marginalization = true;
    _config.min_depth           = 0;

    // Set only by --no-maternal-marginalization, the one maternal switch that
    // is left now that marginalisation is the default.  It exists so that
    // `--model discrete` can be told the truth: that model always integrates
    // p_m and therefore cannot honour an opt-out, so without this the tool
    // would silently return a marginalised fit to a user who explicitly asked
    // for the plug-in.  The other direction needs no tracking -- a continuous
    // run marginalises whether or not anyone asked, and a discrete run gets
    // the equivalent treatment regardless.
    bool plug_in_requested = false;

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
        {"ne-profile",           required_argument, 0, 13 },
        {"ne-profile-step",      required_argument, 0, 14 },
        {"per-family",           no_argument,       0, 15 },
        {"min-family-sites",     required_argument, 0, 16 },
        {"per-family-output",    required_argument, 0, 17 },
        {"no-maternal-marginalization", no_argument, 0, 20 },
        {"min-depth",            required_argument, 0, 19 },
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
            case  7 : _config.kimura_seed      = static_cast<uint64_t>(std::stoull(optarg)); break;
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
            case 11 : _config.bin_simulation_file   = optarg;            break;
            case 12 : _config.bin_simulation_n_bins = std::stoi(optarg); break;
            case 13 : _config.ne_profile_file       = optarg;            break;
            case 14 : _config.ne_profile_step       = std::stod(optarg); break;
            case 15 : _config.per_family            = true;              break;
            case 16 : _config.min_family_sites      = std::stoi(optarg); break;
            case 17 : _config.per_family_output_file = optarg;           break;
            case 19 : _config.min_depth = std::stoi(optarg);             break;
            case 20 : _config.maternal_marginalization = false;
                      plug_in_requested = true;                          break;
            case 't': _config.threads   = std::stoi(optarg);             break;
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
    if (!std::isfinite(_config.min_vaf) || !std::isfinite(_config.max_vaf) ||
        _config.min_vaf < 0.0 || _config.max_vaf > 1.0 ||
        _config.min_vaf > _config.max_vaf) {
        throw std::runtime_error("[ne-estimate] Invalid finite VAF window.");
    }
    if (_config.min_ne < 1) _config.min_ne = 1;
    if (_config.max_ne < _config.min_ne) {
        throw std::runtime_error("[ne-estimate] --max-ne must be >= --min-ne.");
    }
    if (_config.threads < 1)            _config.threads = 1;
    if (_config.kimura_bootstrap < 0)   _config.kimura_bootstrap = 0;
    if (!std::isfinite(_config.kimura_trim) ||
        _config.kimura_trim < 0.0 || _config.kimura_trim >= 1.0) {
        throw std::runtime_error("[ne-estimate] --kimura-trim must be finite and in [0, 1).");
    }
    if (_config.top_drift_k < 0) _config.top_drift_k = 0;
    if (_config.bin_simulation_n_bins < 1) _config.bin_simulation_n_bins = 1;
    if (!std::isfinite(_config.ne_profile_step) || _config.ne_profile_step <= 0.0) {
        throw std::runtime_error("[ne-estimate] --ne-profile-step must be finite and > 0.");
    }
    if (_config.min_family_sites < 1) _config.min_family_sites = 1;
    if (_config.per_family && _config.model != "continuous") {
        throw std::runtime_error(
            "[ne-estimate] --per-family currently requires --model continuous.");
    }
    if (_config.min_depth < 0) {
        throw std::runtime_error(
            "[ne-estimate] --min-depth must be >= 0 (0 disables the filter).");
    }
    // The discrete model already integrates the maternal Beta posterior, so it
    // cannot honour an opt-out.  Report that rather than quietly handing back a
    // marginalised fit to someone who asked for the plug-in by name.  Verified
    // on 5 pairs: discrete yields the same Ne=21 and logL=-41.8442 with and
    // without the flag, i.e. the request really is dropped on the floor.
    if (plug_in_requested && _config.model == "discrete") {
        std::cerr << "[ne-estimate] WARNING: --no-maternal-marginalization has no effect "
                     "with --model discrete (that model always integrates the maternal "
                     "frequency over its Beta posterior); the plug-in was NOT selected.\n";
    }
}

namespace {

std::string json_escape(const std::string& value) {
    static constexpr char kHex[] = "0123456789abcdef";
    std::string escaped;
    escaped.reserve(value.size() + 8);
    for (unsigned char character : value) {
        switch (character) {
            case '"':  escaped += "\\\""; break;
            case '\\': escaped += "\\\\"; break;
            case '\b': escaped += "\\b"; break;
            case '\f': escaped += "\\f"; break;
            case '\n': escaped += "\\n"; break;
            case '\r': escaped += "\\r"; break;
            case '\t': escaped += "\\t"; break;
            default:
                if (character < 0x20) {
                    escaped += "\\u00";
                    escaped += kHex[character >> 4];
                    escaped += kHex[character & 0x0f];
                } else {
                    escaped += static_cast<char>(character);
                }
        }
    }
    return escaped;
}

void emit_json_number(std::ostream& out, double value) {
    if (std::isfinite(value)) {
        out << std::setprecision(8) << value;
    } else if (std::isnan(value)) {
        out << "\"NaN\"";
    } else {
        out << "\"" << (value < 0.0 ? "-Infinity" : "Infinity") << "\"";
    }
}

// Per-family Ne summary statistics.
//
// The median and the clipped count are first-class results here, not
// decoration.  A family whose marginal likelihood has no interior maximum gets
// pinned at --max-ne, and under the default maternal marginalisation that is
// common at shallow maternal depth: measured on 20 families x 12 sites at true
// Ne=30 with d_C=2000, 72% of families pin at the boundary when d_M=30 (42% at
// d_M=60, 26% at d_M=100, 10% at d_M=200, <=4% for d_M>=500), against 0% for
// the plug-in at those depths, because marginalising p_M removes the
// read-noise term that the plug-in misreads as drift and so pushes the
// likelihood towards LARGER Ne.  On such a cohort the mean quotes the search
// boundary back to the user as if it were an estimate -- mean 64.9 against
// median 41.4 for a truth of 30.  The mean only becomes the right summary once
// the clip rate is low.
//
// The pinning is a FEW-SITES failure more than a depth failure.  Holding
// d_M=200 and d_C=2000 fixed and growing the sites per family (8 replicates of
// nested truths, so every site count shares one bottleneck realisation) the
// clip rate runs 9% -> 0% -> 0% and the median bias +11.2% -> +2.9% -> +1.5%
// for 12 -> 48 -> 192 sites, converging on the cohort-level +1.7..2.4%.  The
// lever for a family-level study is therefore more segregating sites per
// family, not more maternal depth once d_M >= 200.
//
// Note the clipped families are still counted in both statistics: excluding
// them would silently drop exactly the families the depth made unestimable,
// which is the opposite of what a comparability diagnostic should do.
struct FamilyNeSummary {
    size_t n_estimated = 0;
    size_t n_skipped   = 0;
    size_t n_clipped   = 0;
    double mean_ne     = 0.0;
    double median_ne   = 0.0;
};

FamilyNeSummary summarize_family_ne(const std::vector<NeEstimator::FamilyResult>& families) {
    FamilyNeSummary summary;
    std::vector<double> ne;
    ne.reserve(families.size());

    for (const auto& fr : families) {
        if (fr.skipped) { ++summary.n_skipped; continue; }
        ++summary.n_estimated;
        if (fr.ci_low_clipped || fr.ci_high_clipped) ++summary.n_clipped;
        ne.push_back(fr.ne);
    }

    if (!ne.empty()) {
        summary.mean_ne = std::accumulate(ne.begin(), ne.end(), 0.0) / static_cast<double>(ne.size());
        std::sort(ne.begin(), ne.end());
        const size_t mid = ne.size() / 2;
        summary.median_ne = (ne.size() % 2 == 1) ? ne[mid] : 0.5 * (ne[mid - 1] + ne[mid]);
    }

    return summary;
}

}  // namespace

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
        << "  \"Trio_Founder_Mismatch_Skipped\": " << r.load_stats.trio_founder_mismatch_skipped << ",\n"
        << "  \"Trio_Founder_Hom_Skipped\":      " << r.load_stats.trio_founder_hom_skipped << ",\n"
        << "  \"Min_Depth_Skipped\":             " << r.load_stats.min_depth_skipped << ",\n"
        << "  \"Estimator\":       \"MMLE (composite marginal likelihood)\",\n"
        << "  \"Max_Marginal_LogLik\": " << std::setprecision(8) << r.max_log_lik << ",\n"
        << "  \"Model\":           \"" << json_escape(_config.model) << "\",\n"
        << "  \"Min_VAF\":         " << _config.min_vaf << ",\n"
        << "  \"Max_VAF\":         " << _config.max_vaf << ",\n"
        << "  \"Search_Min_Ne\":   " << _config.min_ne  << ",\n"
        << "  \"Search_Max_Ne\":   " << _config.max_ne  << ",\n"
        << "  \"Min_Depth\":       " << r.min_depth     << ",\n"
        << "  \"Maternal_Marginalization\": "
        << (r.maternal_marginalization ? "true" : "false") << ",\n"
        // Depth descriptors, descriptive only -- see Result::mean_mother_dp in
        // ne_estimate.h for why no aggregate of these is a plug-in
        // back-correction handle.  Ambient precision is already 8 (set by
        // Max_Marginal_LogLik above) and the format is defaultfloat.
        << "  \"Mean_Mother_DP\":     " << std::setprecision(8) << r.mean_mother_dp     << ",\n"
        << "  \"Mean_Child_DP\":      " << r.mean_child_dp      << ",\n"
        << "  \"Harmonic_Mother_DP\": " << r.harmonic_mother_dp << ",\n"
        << "  \"Harmonic_Child_DP\":  " << r.harmonic_child_dp;

    if (r.kimura.computed) {
        out << ",\n"
            << "  \"Kimura_Cross_Check\": {\n"
            << "    \"b\":             " << std::setprecision(8) << r.kimura.b         << ",\n"
            << "    \"Ne_Kimura\":     ";
        emit_json_number(out, r.kimura.ne_kimura);
        out << ",\n";

        if (r.kimura.ci_computed) {
            out << "    \"b_CI_95_Low\":         "; emit_json_number(out, r.kimura.b_ci_low);          out << ",\n"
                << "    \"b_CI_95_High\":        "; emit_json_number(out, r.kimura.b_ci_high);         out << ",\n"
                << "    \"Ne_Kimura_CI_95_Low\" :"; emit_json_number(out, r.kimura.ne_kimura_ci_low);  out << ",\n"
                << "    \"Ne_Kimura_CI_95_High\":"; emit_json_number(out, r.kimura.ne_kimura_ci_high); out << ",\n"
                << "    \"N_Bootstrap\":         " << r.kimura.n_bootstrap     << ",\n"
                << "    \"Bootstrap_Seed\":      " << r.kimura.bootstrap_seed  << ",\n";
        }

        out << "    \"N_Informative\": " << r.kimura.n_informative << ",\n";

        // ---- Trimmed Kimura diagnostic ----
        if (r.kimura.trimmed_computed) {
            out << "    \"Trimmed_Kimura\": {\n"
                << "      \"Trim_Frac\":        " << std::setprecision(4) << r.kimura.trim_frac   << ",\n"
                << "      \"N_After_Trim\":     " << r.kimura.n_after_trim  << ",\n"
                << "      \"b_Trimmed\":        "; emit_json_number(out, r.kimura.b_trimmed);                  out << ",\n"
                << "      \"Ne_Kimura_Trimmed\": "; emit_json_number(out, r.kimura.ne_kimura_trimmed);          out << "\n"
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

        out << "    \"Note\":          \"" << json_escape(r.kimura.note) << "\",\n"
            << "    \"Method\":        \"Wonnapinij 2008/2010 with sampling-error correction; "
            << "single-generation Ne = 1 / (1 - b); 95% CI by non-parametric pair-level bootstrap\"\n"
            << "  }";
    }

    // Per-family estimates (when --per-family is set).
    if (!r.family_results.empty()) {
        const FamilyNeSummary fam_summary = summarize_family_ne(r.family_results);
        out << ",\n  \"Per_Family_Estimates\": [\n";
        bool first_fam = true;
        for (const auto& fr : r.family_results) {
            if (fr.skipped) continue;
            if (!first_fam) out << ",\n";
            first_fam = false;
            out << "    {\n"
                << "      \"FAM_ID\":            \"" << json_escape(fr.fam_id) << "\",\n"
                << "      \"Mother_ID\":         \"" << json_escape(fr.mother_id) << "\",\n"
                << "      \"N_Children\":        " << fr.n_children << ",\n"
                << "      \"N_Sites\":           " << fr.n_pairs << ",\n"
                << "      \"N_Informative\":     " << fr.n_informative << ",\n"
                << "      \"Ne\":                " << std::fixed << std::setprecision(2) << fr.ne << ",\n"
                << "      \"CI_95_Low\":         " << fr.ci_low << ",\n"
                << "      \"CI_95_High\":        " << fr.ci_high << ",\n"
                << std::defaultfloat
                << "      \"CI_Low_Clipped\":    " << (fr.ci_low_clipped ? "true" : "false") << ",\n"
                << "      \"CI_High_Clipped\":   " << (fr.ci_high_clipped ? "true" : "false") << ",\n"
                << "      \"Max_Marginal_LogLik\": " << std::setprecision(8) << fr.max_log_lik << ",\n"
                // Depth descriptors at 8 significant digits under defaultfloat, matching
                // the top-level Mean_/Harmonic_ DP contract. This was setprecision(2),
                // which rendered e.g. 674.75 as 6.7e+02 and silently dropped precision.
                << "      \"Mean_Mother_DP\":    " << std::setprecision(8) << fr.mean_mother_dp << ",\n"
                << "      \"Mean_Child_DP\":     " << fr.mean_child_dp << ",\n"
                << "      \"Harmonic_Mother_DP\": " << fr.harmonic_mother_dp << ",\n"
                << "      \"Harmonic_Child_DP\":  " << fr.harmonic_child_dp;
            if (fr.kimura.computed) {
                out << ",\n      \"Kimura_Cross_Check\": {\n"
                    << "        \"b\":             " << std::setprecision(8) << fr.kimura.b << ",\n"
                    << "        \"Ne_Kimura\":     ";
                emit_json_number(out, fr.kimura.ne_kimura);
                out << ",\n";
                if (fr.kimura.ci_computed) {
                    out << "        \"b_CI_95_Low\":         ";
                    emit_json_number(out, fr.kimura.b_ci_low);
                    out << ",\n        \"b_CI_95_High\":        ";
                    emit_json_number(out, fr.kimura.b_ci_high);
                    out << ",\n        \"Ne_Kimura_CI_95_Low\": ";
                    emit_json_number(out, fr.kimura.ne_kimura_ci_low);
                    out << ",\n        \"Ne_Kimura_CI_95_High\": ";
                    emit_json_number(out, fr.kimura.ne_kimura_ci_high);
                    out << ",\n        \"N_Bootstrap\":         " << fr.kimura.n_bootstrap
                        << ",\n        \"Bootstrap_Seed\":      " << fr.kimura.bootstrap_seed
                        << ",\n";
                }
                out << "        \"N_Informative\": " << fr.kimura.n_informative << ",\n"
                    << "        \"Sampling_Corrected\": true\n"
                    << "      }";
            }
            if (!fr.warning.empty()) {
                out << ",\n      \"Warning\":         \"" << json_escape(fr.warning) << "\"";
            }
            out << "\n    }";
        }
        out << "\n  ],\n"
            << "  \"Per_Family_Summary\": {\n"
            << "    \"N_Families_Estimated\":  " << fam_summary.n_estimated << ",\n"
            << "    \"N_Families_Skipped\":    " << fam_summary.n_skipped   << ",\n"
            // Families whose CI hit the search boundary: their Ne is a bound
            // rather than an estimate, and they inflate Mean_Ne towards
            // --max-ne.  Median_Ne is the summary to quote while this is
            // non-zero (see summarize_family_ne).
            << "    \"N_Families_CI_Clipped\": " << fam_summary.n_clipped   << ",\n"
            // Same fixed two decimals as Per_Family_Estimates[].Ne above, so a
            // reader can recompute both summaries from the array and reproduce
            // the printed values.  emit_json_number()'s 8 significant digits
            // would carry more precision than the entries themselves have, and
            // the recomputation would then disagree in the third decimal --
            // measured 3.4447831 against 3.446 from the rounded entries.
            << "    \"Mean_Ne\":               " << std::fixed << std::setprecision(2)
            << fam_summary.mean_ne << ",\n"
            << "    \"Median_Ne\":             " << fam_summary.median_ne
            << std::defaultfloat << "\n  }";
    }

    out << "\n}\n";
}

NeEstimator::Result NeEstimator::run() {
    NeEstimator::LoadStats load_stats;
    std::vector<PairData> data = load_pairs(_config.input_tsv,
                                            _config.min_vaf,
                                            _config.max_vaf,
                                            &load_stats,
                                            _config.min_depth);
    if (data.empty()) {
        std::string reason =
            "[ne-estimate] No PASS pairs survived the maternal-VAF filter";
        if (_config.min_depth > 0) {
            reason += " and the --min-depth " + std::to_string(_config.min_depth)
                    + " gate";
        }
        throw std::runtime_error(reason + "; nothing to fit.");
    }

    // Model-form switches threaded down to the per-row likelihood.  The
    // discrete model and the trio rows ignore `opts` because both already
    // integrate the maternal frequency; _parse_args warns about the discrete
    // case when the user asked for the switch explicitly, so it is never a
    // silent no-op for someone who chose it.
    const ModelOptions opts{_config.maternal_marginalization};

    std::cerr << "[ne-estimate] Fitting Ne on " << data.size()
              << " pairs (maternal VAF in [" << _config.min_vaf
              << ", " << _config.max_vaf << "], model=" << _config.model << ").\n";
    if (load_stats.trio_founder_mismatch_skipped > 0 ||
        load_stats.trio_founder_hom_skipped > 0) {
        std::cerr << "[ne-estimate] Dropped "
                  << (load_stats.trio_founder_mismatch_skipped +
                      load_stats.trio_founder_hom_skipped)
                  << " trio rows with a homoplasmic grandmother: "
                  << load_stats.trio_founder_mismatch_skipped
                  << " with segregating descendants (model-incompatible, -inf) and "
                  << load_stats.trio_founder_hom_skipped
                  << " with all-homoplasmic descendants (Ne-independent constant).\n";
    }
    // Remaining non-default levers are logged only when set.
    if (_config.min_depth > 0) {
        std::cerr << "[ne-estimate] --min-depth " << _config.min_depth
                  << " dropped " << load_stats.min_depth_skipped
                  << " pair(s) whose mother or child depth fell below the gate.\n";
    }
    // Rows the switch can actually reach.  Trio rows (has_g == 1) integrate
    // p_M over the grandmother-informed posterior whatever the switch says,
    // and --model discrete integrates p_M as well, so on an all-trio
    // continuous run the opt-out is a genuine no-op.  Announcing "the maternal
    // frequency is the plug-in" there was simply false: measured on
    // tests/data/ne_pipeline/smoke_trio.tsv (300/300 rows HAS_G=1, 289
    // retained pairs) the opt-out arm marginalised every pair and produced a
    // fit bit-identical to the default arm (Ne=23.1765, logL=-2516.298), yet
    // it printed the plug-in banner and echoed Maternal_Marginalization=false.
    // A mixed input pins the scope from the other side: 20 trio + 20
    // two-generation rows do change the fit (logL -344.35218 default against
    // -244.66186 opted out).  Counting the reachable rows once here lets both
    // the banner and the JSON echo below describe what was APPLIED.
    size_t n_two_gen = 0;
    for (const PairData& pd : data) {
        if (pd.has_g == 0) ++n_two_gen;
    }
    const bool plug_in_applied = (_config.model == "continuous")
                                && !_config.maternal_marginalization
                                && n_two_gen > 0;

    // The maternal treatment IS reported either way, because it decides
    // whether the Ne below is comparable across coverage -- that belongs in a
    // default invocation's stderr as much as in an opted-out one.
    if (_config.model == "continuous") {
        if (_config.maternal_marginalization) {
            std::cerr << "[ne-estimate] Maternal heteroplasmy is integrated out against "
                         "its Beta(m_alt+1, m_ref+1) read-sampling posterior (the default; "
                         "--no-maternal-marginalization selects the plug-in); "
                         "two-generation rows only.\n";
        } else if (plug_in_applied) {
            std::cerr << "[ne-estimate] --no-maternal-marginalization: the maternal "
                         "frequency is the plug-in p_m = m_alt / m_dp for the "
                      << n_two_gen
                      << " two-generation pair(s), so the mother's read noise is "
                         "folded into drift.  The fit targets Ne_app with "
                         "1/Ne_app = 1/Ne + 1/d_M and is NOT comparable across cohorts "
                         "or families sequenced at different depth.\n";
        } else {
            std::cerr << "[ne-estimate] --no-maternal-marginalization has NO EFFECT on "
                         "this input: all " << data.size() << " retained pair(s) are "
                         "three-generation rows, which always integrate p_M over the "
                         "grandmother-informed posterior.\n";
        }
    }

    // Warn when Kimura-specific options are set without --cross-check kimura.
    if (!_config.kimura_check) {
        if (_config.kimura_trim > 0.0) {
            std::cerr << "[ne-estimate] WARNING: --kimura-trim "
                      << _config.kimura_trim
                      << " is ignored without --cross-check kimura.\n";
        }
        if (_config.top_drift_k > 0) {
            std::cerr << "[ne-estimate] WARNING: --top-drift-k "
                      << _config.top_drift_k
                      << " is ignored without --cross-check kimura.\n";
        }
        if (_config.kimura_bootstrap != 1000) {
            std::cerr << "[ne-estimate] WARNING: --kimura-bootstrap "
                      << _config.kimura_bootstrap
                      << " is ignored without --cross-check kimura.\n";
        }
    }

    Result r = estimate(data, _config.min_ne, _config.max_ne, _config.threads,
                        _config.model == "continuous", opts);
    r.load_stats = load_stats;
    r.min_depth  = _config.min_depth;
    // Authoritative echo of what was *applied*, not merely requested.  Two
    // things can make the request and the application differ, and both are
    // reflected here: --model discrete integrates p_M itself, so it reports
    // false even though the switch is on by default (_parse_args warns on
    // stderr only when the user asked for it by name); and an opt-out on an
    // all-trio input reaches no row at all, so p_M WAS marginalised for every
    // retained pair and the honest echo is true.  The previous form,
    // `_config.maternal_marginalization && _config.model == "continuous"`,
    // handled the first case and reported `false` for the second -- i.e. it
    // echoed the request, contradicting this comment.
    r.maternal_marginalization = (_config.model == "continuous") && !plug_in_applied;

    // Cohort-level depth descriptors, computed on exactly the rows that
    // entered the likelihood.  load_pairs() guarantees m_dp > 0 and c_dp > 0
    // for every retained row, so the reciprocals below are always finite.
    {
        const double n = static_cast<double>(data.size());
        double sum_m = 0.0, sum_c = 0.0, sum_inv_m = 0.0, sum_inv_c = 0.0;
        for (const PairData& pd : data) {
            sum_m     += static_cast<double>(pd.m_dp);
            sum_c     += static_cast<double>(pd.c_dp);
            sum_inv_m += 1.0 / static_cast<double>(pd.m_dp);
            sum_inv_c += 1.0 / static_cast<double>(pd.c_dp);
        }
        r.mean_mother_dp     = sum_m / n;
        r.mean_child_dp      = sum_c / n;
        r.harmonic_mother_dp = (sum_inv_m > 0.0) ? n / sum_inv_m : 0.0;
        r.harmonic_child_dp  = (sum_inv_c > 0.0) ? n / sum_inv_c : 0.0;
    }

    if (_config.kimura_check) {
        r.kimura = compute_kimura_check(data,
                                        _config.kimura_bootstrap,
                                        _config.kimura_seed,
                                        _config.kimura_trim,
                                        _config.top_drift_k);
    }

    // ---------------------------------------------------------------
    // Optional: per-family Ne estimation.
    // Groups pairs by (FAM_ID, MOTHER_ID) and estimates Ne independently
    // for each family using the continuous MMLE.
    // ---------------------------------------------------------------
    if (_config.per_family) {
        auto families = group_into_families(data);
        std::cerr << "[ne-estimate] Per-family mode: " << families.size()
                  << " families detected from " << data.size() << " pairs.\n";

        r.family_results = estimate_all_families(
            families, _config.min_ne, _config.max_ne,
            _config.min_family_sites, _config.threads, opts);

        // Optional: per-family Kimura cross-check.
        if (_config.kimura_check) {
            for (size_t fi = 0; fi < r.family_results.size(); ++fi) {
                auto& fr = r.family_results[fi];
                if (fr.skipped) continue;
                fr.kimura = compute_family_kimura_check(
                    families[fi], _config.kimura_bootstrap,
                    _config.kimura_seed, _config.kimura_trim);
            }
        }

        // Summary stats.  Reported as mean AND median plus the clipped count,
        // because a family with no interior likelihood maximum is pinned at
        // --max-ne and the mean then quotes the search boundary back as if it
        // were an estimate (see summarize_family_ne for the measured clip
        // rates by maternal depth).
        const FamilyNeSummary fam_summary = summarize_family_ne(r.family_results);
        std::cerr << "[ne-estimate] Per-family results: "
                  << fam_summary.n_estimated << " families estimated, "
                  << fam_summary.n_skipped << " skipped (< "
                  << _config.min_family_sites << " informative sites).\n";
        if (fam_summary.n_estimated > 0) {
            std::cerr << "[ne-estimate] Per-family Ne: mean = "
                      << std::setprecision(4) << fam_summary.mean_ne
                      << ", median = " << fam_summary.median_ne << "\n";
        }
        if (fam_summary.n_clipped > 0) {
            std::cerr << "[ne-estimate] WARNING: " << fam_summary.n_clipped
                      << " of " << fam_summary.n_estimated
                      << " families hit the Ne search boundary ["
                      << _config.min_ne << ", " << _config.max_ne
                      << "]; their Ne is a bound, not an estimate, and they "
                         "pull the mean above the median. Quote the median, or "
                         "widen --max-ne / raise --min-depth. The usual cause "
                         "is shallow maternal depth: with p_m marginalised out "
                         "(the default) a few-site family's likelihood can rise "
                         "monotonically in Ne and have no interior maximum.\n";
        }

        // Per-family TSV output.
        if (!_config.per_family_output_file.empty()) {
            std::ofstream fam_out(_config.per_family_output_file);
            if (!fam_out.is_open()) {
                throw std::runtime_error(
                    "[ne-estimate] Failed to open --per-family-output: "
                    + _config.per_family_output_file);
            }
            if (!_cmdline_string.empty()) fam_out << _cmdline_string << "\n";
            fam_out << "#mitoquest_version="  << MITOQUEST_VERSION << "\n"
                    << "#model="              << _config.model     << "\n"
                    << "#min_vaf="            << _config.min_vaf   << "\n"
                    << "#max_vaf="            << _config.max_vaf   << "\n"
                    << "#min_ne="             << _config.min_ne    << "\n"
                    << "#max_ne="             << _config.max_ne    << "\n"
                    << "#n_families="         << r.family_results.size() << "\n"
                    << "#n_families_estimated=" << fam_summary.n_estimated << "\n"
                    << "#n_families_skipped="  << fam_summary.n_skipped   << "\n"
                    // Boundary-pinned families, and the two summaries that
                    // differ exactly because of them: while the clipped count
                    // is non-zero, #median_family_ne is the comparable one.
                    << "#n_families_ci_clipped=" << fam_summary.n_clipped  << "\n"
                    << "#mean_family_ne="      << std::setprecision(8) << fam_summary.mean_ne   << "\n"
                    << "#median_family_ne="    << std::setprecision(8) << fam_summary.median_ne << "\n"
                    << "#population_ne="       << std::setprecision(8) << r.ne      << "\n"
                    << "#population_ne_ci_low=" << std::setprecision(8) << r.ci_low  << "\n"
                    << "#population_ne_ci_high=" << std::setprecision(8) << r.ci_high << "\n";
            _write_family_tsv(r.family_results, fam_out);
            std::cerr << "[ne-estimate] Wrote per-family TSV to "
                      << _config.per_family_output_file << "\n";
        }
    }

    // ---------------------------------------------------------------
    // Optional: per-bin observed-vs-theoretical drift summary TSV.
    // Computed on the same pair set used by the MMLE; the analytical
    // Kimura predictions p_m(1 - p_m) / Ne are derived from the fitted
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
                << "#max_marginal_log_lik=" << std::setprecision(8) << r.max_log_lik << "\n";
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

    // ---------------------------------------------------------------
    // Optional: Ne-profile TSV (dual-objective Ne scan).
    // For each Ne candidate on a fine grid, score it under both the MMLE
    // marginal log-likelihood and the Kimura per-pair SSR.  This lets the
    // user see which Ne each of the two estimators in the program
    // prefers, and whether they converge to the same optimum.
    // ---------------------------------------------------------------
    if (!_config.ne_profile_file.empty()) {
        const int cache_size = required_cache_size(data, _config.max_ne);
        LogFactorial lf(cache_size);
        const bool continuous = (_config.model == "continuous");
        const std::vector<NeProfileRow> profile = compute_ne_profile(
            data, lf,
            static_cast<double>(_config.min_ne),
            static_cast<double>(_config.max_ne),
            _config.ne_profile_step,
            _config.threads, continuous, opts);

        // Locate the best Ne under each metric (already encoded in the
        // normalised columns, but we surface them in metadata for the
        // plotting script).
        double best_ne_mmle    = profile.empty() ? r.ne : profile.front().ne_candidate;
        double best_ne_kimura = profile.empty() ? std::numeric_limits<double>::quiet_NaN()
                                                 : profile.front().ne_candidate;
        double best_ll        = -std::numeric_limits<double>::infinity();
        double best_ssr       =  std::numeric_limits<double>::infinity();
        for (const auto& row : profile) {
            if (row.mmle_log_lik > best_ll)  { best_ll  = row.mmle_log_lik; best_ne_mmle   = row.ne_candidate; }
            if (row.kimura_ssr   < best_ssr) { best_ssr = row.kimura_ssr;   best_ne_kimura = row.ne_candidate; }
        }
        const double best_ne_kimura_analytic = kimura_ssr_best_ne(data);

        std::ofstream prof_out(_config.ne_profile_file);
        if (!prof_out.is_open()) {
            throw std::runtime_error(
                "[ne-estimate] Failed to open --ne-profile output: "
                + _config.ne_profile_file);
        }

        if (!_cmdline_string.empty()) prof_out << _cmdline_string << "\n";
        prof_out << "#mitoquest_version="      << MITOQUEST_VERSION   << "\n"
                 << "#model="                  << _config.model       << "\n"
                 << "#min_vaf="                << _config.min_vaf     << "\n"
                 << "#max_vaf="                << _config.max_vaf     << "\n"
                 << "#n_pairs_used="           << r.n_pairs           << "\n"
                 << "#profile_min_ne="         << _config.min_ne      << "\n"
                 << "#profile_max_ne="         << _config.max_ne      << "\n"
                 << "#profile_step="           << _config.ne_profile_step << "\n"
                 << "#fitted_ne_mmle="          << std::setprecision(8) << r.ne      << "\n"
                 << "#fitted_ne_mmle_ci_low="   << std::setprecision(8) << r.ci_low  << "\n"
                 << "#fitted_ne_mmle_ci_high="  << std::setprecision(8) << r.ci_high << "\n"
                 << "#max_marginal_log_lik="    << std::setprecision(8) << r.max_log_lik << "\n"
                 << "#best_ne_mmle_on_grid="    << std::setprecision(8) << best_ne_mmle   << "\n"
                 << "#best_ne_kimura_on_grid=" << std::setprecision(8) << best_ne_kimura << "\n"
                 << "#best_ne_kimura_analytic="<< std::setprecision(8) << best_ne_kimura_analytic << "\n";
        if (r.kimura.computed) {
            prof_out << "#kimura_b="            << std::setprecision(8) << r.kimura.b         << "\n"
                     << "#kimura_ne_moments="   << std::setprecision(8) << r.kimura.ne_kimura << "\n";
            if (r.kimura.ci_computed) {
                prof_out << "#kimura_ne_ci_low="  << std::setprecision(8) << r.kimura.ne_kimura_ci_low  << "\n"
                         << "#kimura_ne_ci_high=" << std::setprecision(8) << r.kimura.ne_kimura_ci_high << "\n";
            }
        }

        prof_out << "ne_candidate\tmmle_log_lik\tmmle_delta_2ll"
                    "\tkimura_ssr\tkimura_norm_ssr\n";
        prof_out << std::setprecision(8);
        for (const auto& row : profile) {
            prof_out << row.ne_candidate     << "\t"
                     << row.mmle_log_lik     << "\t"
                     << row.mmle_delta_2ll   << "\t"
                     << row.kimura_ssr       << "\t"
                     << row.kimura_norm_ssr  << "\n";
        }
        std::cerr << "[ne-estimate] Wrote Ne-profile TSV ("
                  << profile.size() << " grid points) to "
                  << _config.ne_profile_file << "\n";
        std::cerr << "[ne-estimate]   Best Ne on grid: MMLE = "
                  << best_ne_mmle
                  << ", Kimura SSR = " << best_ne_kimura;
        if (std::isfinite(best_ne_kimura_analytic)) {
            std::cerr << " (analytic = " << best_ne_kimura_analytic << ")";
        }
        std::cerr << "\n";
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
              << "), max marginal logL = " << r.max_log_lik
              << " on " << r.n_pairs << " pairs.\n";

    // Runtime depth advisory for the DEFAULT path.
    //
    // The opt-out banner above tells the user in terms that the plug-in fit is
    // depth-dependent.  The default banner said nothing about depth at all, so
    // making marginalisation the default quietly removed the only depth
    // sensitivity message a shallow-coverage run would print; --help carries
    // the measured numbers, but help is read once and stderr on every run.
    //
    // Marginalising removes the 1/Ne_app = 1/Ne + 1/d_M bias, it does not make
    // the fit depth-blind: the Beta(1,1) prior on the maternal true frequency
    // is over-dispersed against a cohort whose p_M is truncated away from 0
    // and 1, and the resulting over-estimate decays with maternal depth.
    // Measured at true Ne=30, d_C=2000, 240 pairs per cohort and 20 replicate
    // cohorts per depth (.cache/depth_comparability.py, re-run 2026-09-04
    // against the post-flip binary), the marginalised cohort fit carries
    // +14.1% bias at d_M=30 against |bias| <= 2.4% at every d_M >= 60 --
    // individually +1.9 / +1.7 / +2.1 / -2.1 / +2.4 / +0.9 / +1.6% at
    // d_M = 60 / 100 / 200 / 500 / 1000 / 2000 / 5000.  An earlier revision
    // quoted "+1.6..2.4% for d_M >= 500", which dropped the -2.1% at 500 and
    // the +0.9% at 2000 and so mis-signed the former.
    // It never pins at the search boundary at ANY depth, including d_M=30, so
    // the advisory is about residual bias; boundary pinning belongs to few-row
    // fits and is reported by the per-family clip WARNING instead.
    //
    // Over the two-generation rows only: trio rows (has_g == 1) integrate p_M
    // over the grandmother-informed posterior, the switch is a no-op for them,
    // and quoting their maternal depth here would advise on a quantity this
    // code path does not use.
    if (r.maternal_marginalization) {
        double sum_inv_m2 = 0.0;
        size_t n2 = 0;
        for (const PairData& pd : data) {
            if (pd.has_g != 0) continue;
            sum_inv_m2 += 1.0 / static_cast<double>(pd.m_dp);
            ++n2;
        }
        const double harm_m2 = (n2 == 0 || !(sum_inv_m2 > 0.0)) ? 0.0 : static_cast<double>(n2) / sum_inv_m2;
        if (harm_m2 > 0.0 && harm_m2 < 60.0) {
            const std::streamsize prec = std::cerr.precision(4);
            std::cerr << "[ne-estimate] NOTE: harmonic maternal depth over the "
                      << n2 << " two-generation pair(s) is " << harm_m2
                      << ", below 60 -- the shallowest depth at which the default "
                         "marginalised fit was measured to stay comparable across "
                         "coverage (true Ne=30, d_C=2000, 240 pairs, 20 replicate "
                         "cohorts per depth: +14.1% bias at d_M=30 against <= 2.4% "
                         "at every d_M >= 60).  The Ne above is therefore likely "
                         "HIGH by roughly that margin.  Re-run with --min-depth 60 "
                         "or higher to drop the shallow maternal rows.\n";
            std::cerr.precision(prec);
        }
    }

    // Symmetric advisory for the OPT-OUT path.  The banner above states the
    // identity in the abstract; this prices it on the rows actually fitted, so
    // someone who opted out for speed sees the magnitude they accepted instead
    // of having to look it up.  Before this, opting out was a silent downgrade:
    // the only depth-specific number a run ever printed belonged to the path
    // the user had just declined.
    //
    // 1/Ne_app = 1/Ne + 1/d_M rearranges to Ne_app = Ne*d_M/(Ne+d_M), hence
    // Ne_app/d_M = Ne/(Ne+d_M) EXACTLY -- so the fitted value's own shortfall
    // is Ne_app/d_M, computable from the fit and the harmonic depth alone and
    // therefore incapable of going stale the way a hardcoded bias table would.
    //
    // Printed at every depth rather than below a threshold: the number carries
    // its own significance (0.1305% on the shipped 2000x demo cohort, where the
    // plug-in fits Ne_app=2.60917, against 50% at d_M=30 for a true Ne of 30),
    // and an arbitrary cutoff would hide the mid-depth cases where the
    // trade-off is real but not obvious.
    //
    // The REALISED shortfall tracks this without equalling it, because a finite
    // site count contributes its own upward bias that partially cancels.  Same
    // synthetic sweep as the marginalised NOTE above (true Ne=30, d_C=2000, 240
    // pairs per cohort, 20 replicate cohorts per depth): -45.0% realised at
    // d_M=30 against -50.0% from the identity, -21.2% against -23.1% at 100,
    // -7.4% against -5.7% at 500, -0.6% against -1.5% at 2000.  Quoting the
    // identity alone would overstate the deep-coverage penalty by ~2.5x, so
    // both are printed and neither is presented as the answer.
    if (plug_in_applied) {
        double sum_inv_m3 = 0.0;
        size_t n3 = 0;
        for (const PairData& pd : data) {
            if (pd.has_g != 0) continue;
            sum_inv_m3 += 1.0 / static_cast<double>(pd.m_dp);
            ++n3;
        }

        const double harm_m3 = (n3 == 0 || !(sum_inv_m3 > 0.0)) 
                             ? 0.0 : static_cast<double>(n3) / sum_inv_m3;
        if (harm_m3 > 0.0 && r.ne > 0.0 && r.ne < harm_m3) {
            const std::streamsize prec = std::cerr.precision(4);
            std::cerr << "[ne-estimate] NOTE: plug-in shortfall on this input -- the "
                         "switch reached " << n3
                      << " two-generation pair(s) at harmonic maternal depth "
                      << harm_m3 << ", and 1/Ne_app = 1/Ne + 1/d_M puts the value "
                         "above low by Ne_app/d_M = " << 100.0 * r.ne / harm_m3
                      << "% relative to the Ne marginalisation targets.  Realised "
                         "shortfall on synthetic cohorts (true Ne=30, d_C=2000, 240 "
                         "pairs, 20 replicates per depth): -45.0% at d_M=30, -21.2% "
                         "at 100, -7.4% at 500, -0.6% at 2000, against "
                         "-50.0%/-23.1%/-5.7%/-1.5% from this identity -- same sign "
                         "and magnitude, not equal, because a finite site count adds "
                         "its own upward bias.  Re-run without "
                         "--no-maternal-marginalization to marginalise; the two "
                         "fits are not comparable.\n";
            std::cerr.precision(prec);
        } else if (harm_m3 > 0.0 && r.ne >= harm_m3) {
            // Ne_app = Ne*d_M/(Ne+d_M) is strictly below d_M for every finite
            // Ne, and tends to d_M only as Ne -> infinity.  A plug-in fit at or
            // above the harmonic maternal depth therefore has NO finite
            // shortfall to quote -- Ne_app/d_M would read >= 100%.  Found by
            // execution on a 4-row d_M=30 cohort whose plug-in fit pinned at
            // --max-ne 100: the branch above stayed silent, i.e. the advisory
            // vanished precisely in the shallow, boundary-pinned regime where
            // it matters most.  Say so instead of printing nothing.
            const std::streamsize prec = std::cerr.precision(4);
            std::cerr << "[ne-estimate] NOTE: the plug-in was applied to " << n3
                      << " two-generation pair(s) at harmonic maternal depth "
                      << harm_m3 << ", but the fitted value " << r.ne
                      << " is at or above that depth.  Since "
                         "Ne_app = Ne*d_M/(Ne+d_M) stays strictly below d_M for "
                         "every finite Ne, 1/Ne_app = 1/Ne + 1/d_M has no finite "
                         "shortfall to quote here: the fit is boundary-pinned or "
                         "carries no detectable drift signal at this depth.  "
                         "Re-run without --no-maternal-marginalization, and widen "
                         "--max-ne if the CI is clipped.\n";
            std::cerr.precision(prec);
        }
    }

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
        if (!r.kimura.note.empty()) std::cerr << "  [" << r.kimura.note << "]";
        std::cerr << "\n";

        // Print trimmed Kimura result on its own line (after the main Kimura line).
        if (r.kimura.trimmed_computed) {
            if (std::isfinite(r.kimura.ne_kimura_trimmed)) {
                std::cerr << "[ne-estimate] Trimmed Kimura (trim "
                          << r.kimura.trim_frac * 100.0 << "%): b_trimmed = "
                          << r.kimura.b_trimmed
                          << ", Ne_kimura_trimmed = "
                          << r.kimura.ne_kimura_trimmed
                          << " on " << r.kimura.n_after_trim << " pairs\n";
            } else {
                std::cerr << "[ne-estimate] Trimmed Kimura (trim "
                          << r.kimura.trim_frac * 100.0
                          << "%): not computable (too few informative pairs after trim)\n";
            }
        }

        // Print top-K drift outlier diagnostic to stderr for quick inspection.
        if (!r.kimura.top_drift_outliers.empty()) {
            std::cerr << "[ne-estimate] Top-" << r.kimura.top_drift_outliers.size()
                      << " drift outlier pairs (by F_i descending):\n";
            for (size_t oi = 0; oi < r.kimura.top_drift_outliers.size(); ++oi) {
                const auto& d = r.kimura.top_drift_outliers[oi];
                std::cerr << "  #" << (oi + 1)
                          << "  pair_index=" << d.pair_index
                          << "  M_VAF=" << d.m_vaf
                          << "  C_VAF=" << d.c_vaf
                          << "  F_i=" << d.f_i
                          << "\n";
            }
        }

        // Warn when MMLE and Kimura disagree by more than 3x.
        // When trimmed Kimura is available, compare MMLE against the trimmed
        // value (which is the robust estimator); otherwise compare against
        // the untrimmed value and recommend trimming.
        const double ne_kimura_ref = (r.kimura.trimmed_computed &&
                                      std::isfinite(r.kimura.ne_kimura_trimmed))
                                     ? r.kimura.ne_kimura_trimmed
                                     : r.kimura.ne_kimura;
        const bool   has_trim      = r.kimura.trimmed_computed &&
                                     std::isfinite(r.kimura.ne_kimura_trimmed);
        if (ne_kimura_ref > 0 && std::isfinite(ne_kimura_ref)) {
            const double ratio = r.ne / ne_kimura_ref;
            if (ratio > 3.0 || ratio < 1.0/3.0) {
                if (has_trim) {
                    // Trimmed Kimura still disagrees with MMLE.
                    std::cerr << "\n*** WARNING *** Ne_MMLE (" << r.ne
                              << ") and Ne_Kimura_trimmed (" << ne_kimura_ref
                              << ") still disagree by >3x even after trimming "
                              << r.kimura.trim_frac * 100.0 << "% of high-drift pairs.\n"
                              << "    The untrimmed Ne_Kimura was " << r.kimura.ne_kimura
                              << ".\n"
                              << "    Consider increasing --kimura-trim or inspecting the top drift\n"
                              << "    outliers (--top-drift-k 20) for NUMTs / sequencing errors /\n"
                              << "    mixed populations.\n\n";
                } else {
                    std::cerr << "\n*** WARNING *** Ne_MMLE (" << r.ne
                              << ") and Ne_Kimura (" << r.kimura.ne_kimura
                              << ") disagree by >3x.\n"
                              << "    This usually indicates the data contain high-drift outlier pairs\n"
                              << "    (NUMTs / sequencing errors / mixed populations) that collapse the\n"
                              << "    variance-of-moments Kimura estimator but do not affect the MMLE.\n"
                              << "    Recommendation: re-run with --kimura-trim 0.10 --top-drift-k 20\n"
                              << "    to inspect the outlier pairs and compare the trimmed Ne_Kimura\n"
                              << "    against the MMLE.\n\n";
                }
            }
        }

        // When trim was applied and trimmed agrees with MMLE but untrimmed
        // did not, print a reconciliation note.
        if (has_trim && r.kimura.ne_kimura > 0 && std::isfinite(r.kimura.ne_kimura)) {
            const double ratio_untrimmed = r.ne / r.kimura.ne_kimura;
            const double ratio_trimmed   = r.ne / ne_kimura_ref;
            if ((ratio_untrimmed > 3.0 || ratio_untrimmed < 1.0/3.0) &&
                ratio_trimmed <= 3.0 && ratio_trimmed >= 1.0/3.0) {
                std::cerr << "[ne-estimate] NOTE: trimming "
                          << r.kimura.trim_frac * 100.0
                          << "% of high-drift pairs reconciled Ne_Kimura_trimmed ("
                          << ne_kimura_ref << ") with Ne_MMLE (" << r.ne
                          << ").\n";
            }
        }
    }

    return r;
}
