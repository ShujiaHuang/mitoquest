/**
 * @file algorthm.h
 * 
 * 
 * @brief  Contain some main algorithms of BaseVar
 *  
 * @author Shujia Huang
 * @date 2018-08-30
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

template<class ForwardIterator>
inline size_t argmin(ForwardIterator first, ForwardIterator last) {
    return std::distance(first, std::min_element(first, last));
}

template<class ForwardIterator>
inline size_t argmax(ForwardIterator first, ForwardIterator last) {
    return std::distance(first, std::max_element(first, last));
}

// sum the value for all the data which could call '+' operator
template<typename T> // T must be a numeric type
T sum(const std::vector<T> &value) {
    T d(0);
    for (auto x: value) { d += x; }
    
    return d;
}

// mean the value for all the data which could call '+' operator
template<typename T>  // T must be a numeric type
double mean(const std::vector<T> &value) {
    return sum(value) / value.size();
}

// median the value for all the data which could call '+' operator
template<typename T>  // T must be a numeric type
double median(std::vector<T> &value) {
    size_t n = value.size();
    if (n == 0) return 0.0;

    std::sort(value.begin(), value.end());
    if (n % 2 == 0) {
        return (value[n/2-1] + value[n/2]) / 2.0;
    } else {
        return value[n/2];
    }
}

// Helper function for inverse error function used in z-score calculation
// This is a simple approximation - use a statistical library for more precision
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

/**
 *  Perform a Fisher exact test on a 2x2 contingency table. 
 *  The null hypothesis is that the true odds ratio of the 
 *  populations underlying the observations is one.
 * 
 *    n11  n12  | n1_
 *    n21  n22  | n2_
 *   -----------+----
 *    n_1  n_2  | n
 */
double fisher_exact_test(int n11, int n12, int n21, int n22, 
                         bool is_leftside=false, bool is_rightside=false, 
                         bool is_twoside=true);

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
