#include "algorithm.h"

// Calculate the inverse error function
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

std::pair<double, double> calculate_confidence_interval(int x, int n, double confidence_level) {
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

double fisher_exact_test(int n11, int n12, int n21, int n22, TestSide test_side) {
    // Input validation
    if (n11 < 0 || n12 < 0 || n21 < 0 || n22 < 0) {
        throw std::invalid_argument("Contingency table entries must be non-negative");
    }

    double left_pvalue, right_pvalue, twoside_pvalue;
    kt_fisher_exact(n11, n12, n21, n22, 
                    &left_pvalue, 
                    &right_pvalue, 
                    &twoside_pvalue);

    double pvalue = -1.0;
    switch(test_side) {
        case TestSide::LESS:      pvalue = left_pvalue;    break;
        case TestSide::GREATER:   pvalue = right_pvalue;   break;
        case TestSide::TWO_SIDED: pvalue = twoside_pvalue; break;
        default:
            throw std::invalid_argument("Invalid test side specification");
    }

    return pvalue; // -1 is an error
}

double fisher_exact_test(const ContingencyTable& table, TestSide test_side) {
    double left_pvalue, right_pvalue, twoside_pvalue;
    kt_fisher_exact(table.n11, table.n12, table.n21, table.n22,
                    &left_pvalue, &right_pvalue, &twoside_pvalue);
    
    double pvalue = -1.0;
    switch(test_side) {
        case TestSide::LESS:      pvalue = left_pvalue;    break;
        case TestSide::GREATER:   pvalue = right_pvalue;   break;
        case TestSide::TWO_SIDED: pvalue = twoside_pvalue; break;
        default:
            throw std::invalid_argument("Invalid test side specification");
    }

    return pvalue; // -1 is an error
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
        for (size_t j = 0; j < n_allele; j++) {  // for _UNIQ_BASES
            likelihood[i][j] = ind_allele_likelihood[i][j] * obs_allele_freq[j];
            marginal_likelihood[i] += likelihood[i][j];
        }

        // Computed the posterior probability of _UNIQ_BASES, change the value inplace.
        for (size_t j(0); j < n_allele; ++j) { 
            // reset the posterior value
            ind_allele_post_prob[i][j] = likelihood[i][j] / marginal_likelihood[i];
        }
    }

    return;
}

void m_step(const std::vector<std::vector<double>> &ind_allele_post_prob, std::vector<double> &obs_allele_freq) {
    size_t n_sample = ind_allele_post_prob.size();  // depth
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

void EM(const std::vector<std::vector<double>> &ind_allele_likelihood, // n x _UNIQ_BASES.size() matrix, do not change the raw value.
        std::vector<double> &obs_allele_freq,          // retuen value, 1 x _UNIQ_BASES.size(), expect allele frequence, it'll be update inplace here.
        std::vector<double> &log_marginal_likelihood,  // return value
        int iter_num, const float epsilon)
{
    if (iter_num <= 0) {
        throw std::invalid_argument("Iteration number must be positive");
    }
    if (epsilon <= 0) {
        throw std::invalid_argument("Epsilon must be positive");
    }

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
