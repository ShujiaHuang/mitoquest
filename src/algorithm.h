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
double _erfinv(double x) {
    // Abramowitz and Stegun approximation
    const double c0 = 2.515517;
    const double c1 = 0.802853;
    const double c2 = 0.010328;
    const double d1 = 1.432788;
    const double d2 = 0.189269;
    const double d3 = 0.001308;
    
    double y = 0.5 * log(1 - x * x);
    double z = copysign(sqrt(-y), x);
    
    double numerator = c0 + c1 * z + c2 * z * z;
    double denominator = 1 + d1 * z + d2 * z * z + d3 * z * z * z;
    
    return copysign(z - numerator / denominator, x);
}

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
std::pair<double, double> calculate_confidence_interval(int x, int n, double confidence_level = 0.95) {
    // Input validation
    if (n <= 0) {
        throw std::invalid_argument("Sample size must be positive");
    }
    if (x < 0 || x > n) {
        throw std::invalid_argument("Number of x must be between 0 and sample size");
    }
    if (confidence_level <= 0 || confidence_level >= 1) {
        throw std::invalid_argument("Confidence level must be between 0 and 1");
    }
    
    // Calculate z-score for the given confidence level
    double alpha = 1.0 - confidence_level;
    double z = 0.0;
    
    // Approximate z-score based on normal distribution
    // For 95% confidence (alpha = 0.05), z ≈ 1.96
    if (confidence_level == 0.95) {
        z = 1.96;
    } else if (confidence_level == 0.99) {
        z = 2.576;
    } else if (confidence_level == 0.90) {
        z = 1.645;
    } else {
        // For other confidence levels, use this approximation
        // Note: This is an approximation for simplicity
        // A more accurate approach would use a statistical library
        z = sqrt(2.0) * _erfinv(confidence_level);
    }
    
    double p = static_cast<double>(x) / n;  // Sample proportion
    double lower = 0.0, upper = 0.0;
    // Choose method based on the sample size: n
    if (n >= 40) {  // Use Agresti-Coull for larger sample sizes
        // Agresti-Coull interval calculation
        double n_tilde = n + z * z;
        double p_tilde = (x + z * z / 2.0) / n_tilde;
        double interval_half_width = z * sqrt((p_tilde * (1.0 - p_tilde)) / n_tilde);
        
        lower = std::max(0.0, p_tilde - interval_half_width);
        upper = std::min(1.0, p_tilde + interval_half_width);
    } else {  // Use Wilson score for smaller sample sizes
        // Wilson score interval calculation
        double denominator = 1.0 + (z * z / n);
        double center = (p + (z * z) / (2 * n)) / denominator;
        double interval_half_width = z * sqrt((p * (1.0 - p) / n) + (z * z / (4 * n * n))) / denominator;

        lower = std::max(0.0, center - interval_half_width);
        upper = std::min(1.0, center + interval_half_width);
    }

    return std::make_pair(lower, upper);  // Return the confidence interval
}

// Function for chi^2 test
double chi2_test(double chi_sqrt_value, double degree_of_freedom) {
    return kf_gammaq(degree_of_freedom/2.0, chi_sqrt_value/2.0);
}

double norm_dist(double x) {
    return kf_erfc(double(x / std::sqrt(2.0))) / 2.0;
}

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
                         bool is_leftside  = false,
                         bool is_rightside = false,
                         bool is_twoside   = true)
{
    double left_pvalue, right_pvalue, twoside_pvalue;
    double v = kt_fisher_exact(n11, n12, n21, n22, 
                               &left_pvalue, 
                               &right_pvalue, 
                               &twoside_pvalue);

    return (is_twoside) ? twoside_pvalue : (is_leftside ? left_pvalue : right_pvalue);
}

double wilcoxon_ranksum_test(const std::vector<double>& sample1, const std::vector<double>& sample2) {

    size_t n1 = sample1.size(), n2 = sample2.size();

    std::vector<double> combined = sample1;
    combined.insert(combined.end(), sample2.begin(), sample2.end());

    // 排序并分配秩
    std::vector<size_t> rank_idx(combined.size());
    std::iota(rank_idx.begin(), rank_idx.end(), 0.0); // 初始化秩
    std::sort(rank_idx.begin(), rank_idx.end(), [&combined](size_t a, size_t b) { return combined[a] > combined[b]; });

    std::vector<double> rankvalues(combined.size());
    for (size_t i(0); i < rank_idx.size(); ++i) 
        rankvalues[i] = i+1; 

    // 处理重复元素
    double ranksum = 0.0, same_n = 1;
    size_t i;
    for (i = 0; i < rank_idx.size(); ++i) {
        if (i > 0 && combined[rank_idx[i]] != combined[rank_idx[i-1]]) {

            if (same_n > 1) {
                double avg_rank = ranksum / same_n; // 平均秩
                for (size_t j = i - same_n; j < i; ++j) {
                    rankvalues[j] = avg_rank; // 分配平均秩
                }
            }

            // 重置
            same_n  = 1;
            ranksum = 0;
        } else if (i > 0) {
            same_n++;
        }
        ranksum += i+1;
    }

    // 处理最后一组重复
    if (same_n > 1) {
        double avg_rank = ranksum / same_n; // 平均秩
        for (size_t j = i - same_n; j < i; ++j) {
            rankvalues[j] = avg_rank; // 分配平均秩
        }
    }

    // 计算样本1的秩和
    double smp1_ranksum = 0.0;
    for (size_t i = 0; i < rank_idx.size(); ++i) {
        if (rank_idx[i] < n1) {
            smp1_ranksum += rankvalues[i];
        }
    }

    double e = (double)(n1 * (n1 + n2 + 1)) / 2.0;
    double z = (smp1_ranksum - e) / std::sqrt(double(n1*n2*(n1+n2+1))/12.0);
    double p = 2 * norm_dist(std::abs(z));
    
    // 返回秩和检验 pvalue
    return p;
}

/**
 * @brief Calculate the posterior probability of individual allele at each site as
 *        the four A/C/G/T bases.
 * 
 * @param obs_allele_freq        1 x 4 matrix.
 * @param ind_allele_likelihood  n x 4 matrix. n is sample size
 * @param ind_allele_post_prob   n x 4 matrix. allele posterior probabilty for each sample.
 * @param marginal_likelihood    n x 1 matrix. maginal likelihood for each sample.
 * 
 */
void e_step(const std::vector<double> &obs_allele_freq,
            const std::vector<std::vector<double>> &ind_allele_likelihood, 
            std::vector<std::vector<double>> &ind_allele_post_prob,  // return value, update value inplace 
            std::vector<double> &marginal_likelihood)                // return value, calculate value inplace
{
    size_t n_sample = ind_allele_likelihood.size();
    size_t n_allele = obs_allele_freq.size();

    // 'likelihood' is as the same shape as `ind_allele_likelihood`
    std::vector<std::vector<double>> likelihood(n_sample, std::vector<double>(n_allele, 0));

    // reset the raw value to be 0
    marginal_likelihood = std::vector<double>(n_sample, 0);
    for (size_t i(0); i < n_sample; ++i) {
        for (size_t j = 0; j < n_allele; j++) {  // for [A, C, G, T]
            likelihood[i][j] = ind_allele_likelihood[i][j] * obs_allele_freq[j];
            marginal_likelihood[i] += likelihood[i][j];
        }

        // Computed the posterior probability of A/C/G/T for each individual, change the value inplace.
        for (size_t j(0); j < n_allele; ++j) { 
            // reset the posterior value
            ind_allele_post_prob[i][j] = likelihood[i][j] / marginal_likelihood[i];
        }
    }

    return;
}

/**
 * @brief Update observed allele frequency inplace.
 * 
 * @param ind_allele_post_prob 2d-array, n x 4 matrix.
 * @param obs_allele_freq      1d-array, 1 x 4. update the allele frequence inplace and return. 
 * 
 */
void m_step(const std::vector<std::vector<double>> &ind_allele_post_prob, std::vector<double> &obs_allele_freq) {
    size_t n_sample = ind_allele_post_prob.size();
    size_t n_allele = ind_allele_post_prob[0].size();

    // Reset data
    obs_allele_freq = std::vector<double>(n_allele, 0);
    for (size_t j(0); j < n_allele; ++j) {
        for (size_t i(0); i < n_sample; ++i) {
            obs_allele_freq[j] += ind_allele_post_prob[i][j];
        }
        obs_allele_freq[j] /= (double)(n_sample);  // average.
    }

    return;
}

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
void EM(const std::vector<std::vector<double>> &ind_allele_likelihood, // n x 4 matrix, do not change the raw value.
        std::vector<double> &obs_allele_freq,          // retuen value, 1 x 4, expect allele frequence, it'll be update inplace here.
        std::vector<double> &log_marginal_likelihood,  // return value
        int iter_num=100, const float epsilon=0.001)
{
    size_t n_sample = ind_allele_likelihood.size();
    size_t n_allele = obs_allele_freq.size();

    // n x 4 matrix, the same shape as 'ind_allele_likelihood'.
    std::vector<std::vector<double>> ind_allele_post_prob = std::vector<std::vector<double>>(
        n_sample, std::vector<double>(n_allele, 0));

    // It's a 1-d array (n x 1) one sample per value, n is sample size
    std::vector<double> marginal_likelihood = std::vector<double>(n_sample, 0);

    // Update the value of 'ind_allele_post_prob' and compute 'marginal_likelihood' in e_step.
    e_step(obs_allele_freq, ind_allele_likelihood, ind_allele_post_prob, marginal_likelihood);

    log_marginal_likelihood = std::vector<double>(n_sample, 0);
    for (size_t i = 0; i < marginal_likelihood.size(); i++) {
        log_marginal_likelihood[i] = log(marginal_likelihood[i]);
    }
    
    // use 'ind_allele_post_prob' to update 'obs_allele_freq'
    m_step(ind_allele_post_prob, obs_allele_freq);
    while (iter_num--) {
        // e step: update ind_allele_post_prob
        e_step(obs_allele_freq, ind_allele_likelihood, ind_allele_post_prob, marginal_likelihood);

        // m step: update the frequence of observed alleles
        m_step(ind_allele_post_prob, obs_allele_freq); 

        double delta = 0, llh;
        for (size_t i = 0; i < marginal_likelihood.size(); i++) {
            llh = log(marginal_likelihood[i]);
            delta += abs(llh - log_marginal_likelihood[i]);
            log_marginal_likelihood[i] = llh;  // update
        }
        
        // Todo: be careful here!!!
        if (delta < epsilon) break;
    }
    // update the lastest expect_allele_freq
    m_step(ind_allele_post_prob, obs_allele_freq);
    return;
}

#endif
