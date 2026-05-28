/**
 * @file log_factorial.cpp
 * @brief Implementation of the LogFactorial cache and PMFs.
 *
 * @author Shujia Huang (hshujia@qq.com)
 * @date 2026-05-28
 */
#include "log_factorial.h"

#include <cmath>
#include <limits>

namespace mitoquest {

LogFactorial::LogFactorial(int max_n) {
    if (max_n < 1) max_n = 1;
    _cache.assign(static_cast<size_t>(max_n) + 1, 0.0);
    // _cache[0] == log(0!) == 0; build the running sum from there.
    for (int i = 1; i <= max_n; ++i) {
        _cache[i] = _cache[i - 1] + std::log(static_cast<double>(i));
    }
}

double LogFactorial::log_fact(int n) const {
    if (n < 0) return -std::numeric_limits<double>::infinity();
    if (n <= capacity()) return _cache[n];
    // Spill over: fall back to lgamma for indices beyond the cache.
    return std::lgamma(static_cast<double>(n) + 1.0);
}

double LogFactorial::log_comb(int n, int k) const {
    if (k < 0 || k > n || n < 0) {
        return -std::numeric_limits<double>::infinity();
    }
    return log_fact(n) - log_fact(k) - log_fact(n - k);
}

double LogFactorial::log_binomial_pmf(int n, int k, double p) const {
    if (k < 0 || k > n) return -std::numeric_limits<double>::infinity();
    if (p <= 0.0) return (k == 0) ? 0.0 : -std::numeric_limits<double>::infinity();
    if (p >= 1.0) return (k == n) ? 0.0 : -std::numeric_limits<double>::infinity();
    return log_comb(n, k)
         + static_cast<double>(k)     * std::log(p)
         + static_cast<double>(n - k) * std::log1p(-p);
}

double LogFactorial::log_betabinom_pmf(int n, int k,
                                       double alpha, double beta) const {
    if (k < 0 || k > n || n < 0) return -std::numeric_limits<double>::infinity();
    if (alpha <= 0.0 || beta <= 0.0) {
        return -std::numeric_limits<double>::infinity();
    }
    // log Pr(k | n, a, b) = log C(n, k)
    //                     + lgamma(k + a) + lgamma(n - k + b) - lgamma(n + a + b)
    //                     - lgamma(a)     - lgamma(b)         + lgamma(a + b)
    const double k_d = static_cast<double>(k);
    const double n_d = static_cast<double>(n);
    return log_comb(n, k)
         + std::lgamma(k_d + alpha)
         + std::lgamma(n_d - k_d + beta)
         - std::lgamma(n_d + alpha + beta)
         - std::lgamma(alpha)
         - std::lgamma(beta)
         + std::lgamma(alpha + beta);
}

}  // namespace mitoquest
