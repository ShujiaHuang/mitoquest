/**
 * @file algorthm.h
 * 
 * @brief Algorithm functions for mtDNA variant caller
 *  
 * @author Shujia Huang
 * @date 2025-02-12
 * 
 */
#ifndef __INCLUDE_MTCALLER_ALIGORITHM_H__
#define __INCLUDE_MTCALLER_ALIGORITHM_H__

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>  // use 'log' functon
#include <numeric>

#include <htslib/kfunc.h>

template<typename ForwardIterator>
inline size_t argmin(ForwardIterator first, ForwardIterator last) {
    if (first == last) {
        throw std::invalid_argument("Empty range");
    }
    return std::distance(first, std::min_element(first, last));
}

template<typename ForwardIterator>
inline size_t argmax(ForwardIterator first, ForwardIterator last) {
    if (first == last) {
        throw std::invalid_argument("Empty range");
    }
    return std::distance(first, std::max_element(first, last));
}

// sum the value for all the data which could call '+' operator
template<typename T> // T must be a numeric type
typename std::enable_if<std::is_arithmetic<T>::value, T>::type
sum(const std::vector<T> &values) {
    if (values.empty()) {
        return T(0);
    }
    return std::accumulate(values.begin(), values.end(), T(0));
}

// mean the value for all the data which could call '+' operator
template<typename T>  // T must be a numeric type
double mean(const std::vector<T> &values) {
    if (values.empty()) {
        throw std::invalid_argument("Cannot calculate mean of empty vector");
    }
    return static_cast<double>(sum(values)) / values.size();
}

// median the value for all the data which could call '+' operator
template<typename T>  // T must be a numeric type
double median(const std::vector<T> &value) {
    size_t n = value.size();
    if (n == 0) return 0.0;

    std::vector<T> value_copy = value; // make a copy to avoid modifying the original vector
    std::sort(value_copy.begin(), value_copy.end());
    if (n % 2 == 0) {
        return (value_copy[n/2-1] + value_copy[n/2]) / 2.0;
    } else {
        return value_copy[n/2];
    }
}

// standard deviation the value for all the data which could call '+' operator
template<typename T>  // T must be a numeric type
double stddev(const std::vector<T> &values) {
    if (values.size() < 2) {
        throw std::invalid_argument("Cannot calculate standard deviation of less than 2 values");
    }
    double mean_value = mean(values);
    double sum_squared_diff = 0.0;
    for (const auto& value : values) {
        sum_squared_diff += (value - mean_value) * (value - mean_value);
    }
    return std::sqrt(sum_squared_diff / (values.size() - 1));
}

// calculate the standard error of the mean
template<typename T>  // T must be a numeric type
double standard_error(const std::vector<T> &values) {
    if (values.size() < 2) {
        throw std::invalid_argument("Cannot calculate standard error of less than 2 values");
    }
    return stddev(values) / std::sqrt(values.size());
}

// Helper function for inverse error function is used in confidence interval calculations
// Note: This is a simple approximation and may not be accurate for all values.
// For more precision, consider using a statistical library or a more accurate approximation.
double _erfinv(double x);

/**
 * Calculates confidence interval for a proportion using Wilson score or Agresti-Coull method
 * 
 * @param x Number of successes (count of positive outcomes)
 * @param n Sample size
 * @param confidence_level Desired confidence level (default 0.95 for 95% confidence)
 * @return std::pair<double, double> containing (lower bound, upper bound)
 * 
 * Wilson score interval: recommended for small sample sizes (n < 40)
 * CI = (p + z^2 / (2n) ± z * sqrt((p * (1 - p) / n) + z^2 / (4n))) / (1 + z^2 / n)
 * 
 * Agresti-Coull interval: recommended for large sample sizes (n >= 40)
 * CI = (p_tilde ± z * sqrt(p_tilde * (1 - p_tilde) / n_tilde)) / n_tilde
 * 
 */
std::pair<double, double> calculate_confidence_interval(int x, int n, double confidence_level=0.95);

// Function for chi^2 test
double chi2_test(double chi_sqrt_value, double degree_of_freedom);
double norm_dist(double x);

enum class TestSide {
    LESS,
    GREATER,
    TWO_SIDED
};
/**
 *  Perform a Fisher exact test on a 2x2 contingency table. 
 *  The null hypothesis is that the true odds ratio of the 
 *  populations underlying the observations is one.
 *  
 *  Mathmatical note: https://mathworld.wolfram.com/FishersExactTest.html
 * 
 *  @param n11 Value in cell (1,1)
 *  @param n12 Value in cell (1,2)
 *  @param n21 Value in cell (2,1)
 *  @param n22 Value in cell (2,2)
 *  @param test_side Type of test to perform (left-sided, right-sided, or two-sided)
 *  @return p-value for the specified test
 * 
 *    n11  n12  | n1_
 *    n21  n22  | n2_
 *   -----------+----
 *    n_1  n_2  | n
 * 
 *  Example: https://gatk.broadinstitute.org/hc/en-us/articles/360035532152-Fisher-s-Exact-Test
 * 
 */
double fisher_exact_test(int n11, int n12, int n21, int n22, TestSide test_side=TestSide::TWO_SIDED);

struct ContingencyTable {
    int n11, n12, n21, n22;
    ContingencyTable(int a, int b, int c, int d) : n11(a), n12(b), n21(c), n22(d) {
        if (n11 < 0 || n12 < 0 || n21 < 0 || n22 < 0) {
            throw std::invalid_argument("Contingency table entries must be non-negative");
        }
    }
};
double fisher_exact_test(const ContingencyTable& table, TestSide test_side=TestSide::TWO_SIDED);

double wilcoxon_ranksum_test(const std::vector<double>& sample1, const std::vector<double>& sample2);

/**
 * @brief Calculate the posterior probability of individual allele at each site as
 *        the bases.
 * 
 * @param obs_allele_freq        1 x _UNIQ_BASES.size() matrix.
 * @param ind_allele_likelihood  n x _UNIQ_BASES.size() matrix. n is sample size
 * @param ind_allele_post_prob   n x _UNIQ_BASES.size() matrix. allele posterior probabilty for each sample.
 * @param marginal_likelihood    n x _UNIQ_BASES.size() matrix. maginal likelihood for each sample.
 * 
 */
void e_step(const std::vector<double> &obs_allele_freq,
            const std::vector<std::vector<double>> &ind_allele_likelihood, 
            std::vector<std::vector<double>> &ind_allele_post_prob,  // return value, update value inplace 
            std::vector<double> &marginal_likelihood);               // return value, calculate value inplace

/**
 * @brief Update observed allele frequency inplace.
 * 
 * @param ind_allele_post_prob 2d-array, n x _UNIQ_BASES.size() matrix.
 * @param obs_allele_freq      1d-array, 1 x _UNIQ_BASES.size(). update the allele frequence inplace and return. 
 * 
 */
void m_step(const std::vector<std::vector<double>> &ind_allele_post_prob, std::vector<double> &obs_allele_freq);

/**
 * @brief EM algorithm
 * 
 * @param ind_allele_likelihood  n x 4 matrix, n is non-N and non-Indel's sample size (same below), 4 for [A, C, G, T]
 * @param obs_allele_freq        1 x 4 vector, Observed allele frequence for [A, C, G, T]
 * @param iter_num integer, optional.  The lager EM iteration times. default: 100
 * @param epsilon  float, optional. The threshold of likelihood different for EM 
 *                 process. default 0.001.
 * 
 */
void EM(const std::vector<std::vector<double>> &ind_allele_likelihood, // n x _UNIQ_BASES.size() matrix, do not change the raw value.
        std::vector<double> &obs_allele_freq,          // retuen value, 1 x _UNIQ_BASES.size(), expect allele frequence, it'll be update inplace here.
        std::vector<double> &log_marginal_likelihood,  // return value
        int iter_num=100, const float epsilon=0.001);

#endif
