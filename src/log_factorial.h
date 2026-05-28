/**
 * @file log_factorial.h
 * @brief Numerically stable log-factorial cache and discrete probability
 *        mass functions (Binomial, Beta-Binomial) used by likelihood
 *        models such as `mitoquest ne-estimate`.
 *
 * Design goals:
 *   * Cache log(n!) for n in [0, max_n] in O(max_n) time and memory.
 *   * Spill to std::lgamma for indices beyond the cache (slower but safe).
 *   * Boundary-safe PMFs: log(0) is returned as -INFINITY rather than NaN
 *     so that callers can use std::log_sum_exp-style accumulators.
 *
 * This header intentionally has no dependency on any caller; it lives at
 * the same layer as `algorithm.h` and can be reused by any module that
 * needs discrete-likelihood building blocks.
 *
 * @author Shujia Huang (hshujia@qq.com)
 * @date 2026-05-28
 */
#ifndef __MITOQUEST_LOG_FACTORIAL_H__
#define __MITOQUEST_LOG_FACTORIAL_H__

#include <cstddef>
#include <vector>

namespace mitoquest {

/**
 * @brief Pre-computed log(n!) cache plus log-domain Binomial / Beta-Binomial
 *        probability mass functions.
 *
 * All PMFs are returned in log space and follow a strict convention:
 *   * Impossible counts (k < 0 or k > n)            => -INFINITY.
 *   * Boundary p in {0, 1} for the Binomial         => exact 0 or -INFINITY.
 *   * Non-positive shape parameters (alpha, beta)   => -INFINITY.
 *
 * Construction is O(max_n).  Queries are O(1) up to capacity() and
 * O(1)+lgamma overhead beyond it.
 */
class LogFactorial {
public:
    /**
     * @brief Build a cache covering log(0!) ... log(max_n!).
     * @param max_n  largest index to pre-compute.  max_n < 1 is clamped to 1.
     */
    explicit LogFactorial(int max_n);

    /**
     * @brief log(n!).  Returns -INFINITY for n < 0.  For n > capacity() the
     *        result is computed via std::lgamma(n + 1) on the fly.
     */
    double log_fact(int n) const;

    /**
     * @brief log of the binomial coefficient C(n, k).
     *        Returns -INFINITY when k is out of [0, n].
     */
    double log_comb(int n, int k) const;

    /**
     * @brief log Binomial PMF: log Pr(X = k | n trials, success prob p).
     *
     * Uses log1p(-p) to avoid catastrophic cancellation when p is close
     * to 1.  Boundary handling:
     *   p == 0 : returns 0     when k == 0,            -INFINITY otherwise.
     *   p == 1 : returns 0     when k == n,            -INFINITY otherwise.
     */
    double log_binomial_pmf(int n, int k, double p) const;

    /**
     * @brief log Beta-Binomial PMF.
     *
     *   log Pr(k | n, alpha, beta) = log C(n, k)
     *                              + lgamma(k + alpha)
     *                              + lgamma(n - k + beta)
     *                              - lgamma(n + alpha + beta)
     *                              - lgamma(alpha)
     *                              - lgamma(beta)
     *                              + lgamma(alpha + beta)
     *
     * @pre alpha > 0 and beta > 0; otherwise returns -INFINITY.
     */
    double log_betabinom_pmf(int n, int k, double alpha, double beta) const;

    /// Largest n for which log_fact(n) is served from the cache.
    int capacity() const { return static_cast<int>(_cache.size()) - 1; }

private:
    std::vector<double> _cache;  // _cache[i] == log(i!)
};

}  // namespace mitoquest

#endif  // __MITOQUEST_LOG_FACTORIAL_H__
