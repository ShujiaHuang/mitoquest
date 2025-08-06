/**
 * @file basetype.cpp
 * @brief 
 * 
 * @author Shujia Huang
 * @date 2018-08-01
 * 
 */
#include "basetype.h"

//////////////////////////////////////////////////////////////
//// The codes for the member function of BaseType class /////
//////////////////////////////////////////////////////////////
BaseType::BaseType(const BatchInfo *smp_bi, double min_af) {

    // set common values
    /**
     * @brief set the unique bases and the index map
     * 
     * 1. merge the align_bases and BASIC_BASES to make sure the UNIQ_BASES 
     *    at least have the four basic bases ([A, C, G, T])
     * 2. get the unique strings vector: sort by length and then by ASCII
     * 3. set the string of UNIQ_BASES map to array index
     * 4. inital the uniq bases depth to 0
     * 
     */
    std::vector<std::string> tmp_bases(smp_bi->align_bases);
    tmp_bases.insert(tmp_bases.end(), BASIC_BASES.begin(), BASIC_BASES.end());
    _UNIQ_BASES = ngslib::get_unique_strings(tmp_bases);

    // Remove bases that are not valid (e.g., 'N', '*')
    _UNIQ_BASES.erase(
        std::remove_if(_UNIQ_BASES.begin(), _UNIQ_BASES.end(),
            [](const std::string& base) {
                return base[0] == 'N' || base[0] == 'n' || base[0] == '*';
                // Add more conditions as needed
            }
        ),
        _UNIQ_BASES.end()
    );

    // Continue with the rest of the initialization
    for (size_t i(0); i < _UNIQ_BASES.size(); ++i) {
        _B_IDX[_UNIQ_BASES[i]] = i;  // set string of UNIQ_BASES map to array index
        _depth[_UNIQ_BASES[i]] = 0;  // inital the base depth to 0.
    }

    _min_af  = min_af;
    _ref_id  = smp_bi->ref_id;
    _ref_pos = smp_bi->ref_pos;
    _total_depth = 0;
    _allele_likelihood.reserve(smp_bi->align_bases.size());
    for (size_t i(0); i < smp_bi->align_bases.size(); ++i) {

        if (_bases2ref.find(smp_bi->align_bases[i]) == _bases2ref.end()) {
            _bases2ref.insert({smp_bi->align_bases[i], smp_bi->ref_bases[i]});
        }

        double epsilon = exp((smp_bi->align_base_quals[i] - 33) * MLN10TO10); // base error probability
        if (smp_bi->align_bases[i][0] != 'N' && smp_bi->align_bases[i][0] != 'n') { 
            // ignore all the 'N' bases
            _total_depth++;
            _depth[smp_bi->align_bases[i]]++;

            // Initialized the array to {0, 0, 0, 0, ...}, which set allele likelihood for _UNIQ_BASES
            std::vector<double> allele_lh(_UNIQ_BASES.size(), 0);
            for (auto &b: _UNIQ_BASES) {
                // convert the quality phred scale to be the base confident probabilty value
                allele_lh[_B_IDX[b]] = (smp_bi->align_bases[i] == b) ? 1.0 - epsilon : epsilon / (_UNIQ_BASES.size() - 1);
            }

            // 也就是说这个 likelihood 数组将只保留有覆盖的位点信息，避免无效计算。因此 _ind_allele_likelihood 的长度较小
            // 但无所谓，因为除了该值不会传到 basetype classe 之外，仅仅只用于计算突变 
            _allele_likelihood.push_back(allele_lh);  // A 2d-array, total_depth x _UNIQ_BASES.size() matrix
        }
    }
}

BaseType::BaseType(const BaseType &b) {

    this->_UNIQ_BASES   = b._UNIQ_BASES;
    this->_ref_id       = b._ref_id;
    this->_ref_pos      = b._ref_pos;
    this->_active_bases = b._active_bases;

    this->_min_af       = b._min_af;
    this->_var_qual     = b._var_qual;
    this->_af_by_lrt    = b._af_by_lrt;
    this->_allele_likelihood = b._allele_likelihood;

    this->_depth       = b._depth;
    this->_total_depth = b._total_depth;
    this->_B_IDX       = b._B_IDX;
    this->_bases2ref   = b._bases2ref;
}

std::vector<double> BaseType::_set_initial_freq(const std::vector<std::string> &bases) {
    // bases 数组中 A,C,G,T 这四个碱基最多只能各出现一次
    // initial observed allele likelihood for _UNIQ_BASES
    std::vector<double> obs_allele_freq(_UNIQ_BASES.size(), 0);
    if (this->_total_depth > 0) {
        // computed by base count
        for (auto b: bases) {
            if (_B_IDX.find(b) == _B_IDX.end()) {
                throw std::runtime_error("The base " + b + " is not in the _UNIQ_BASES.");
            }
            obs_allele_freq[_B_IDX[b]] = this->_depth[b] / (double)(this->_total_depth);
        }
    }

    return obs_allele_freq;  // 1 x _UNIQ_BASES.size() vector. The allele frequence for _UNIQ_BASES
}

BaseType::AA BaseType::_f(const std::vector<std::string> &bases, int n) {

    AA data;
    Combinations<std::string> c(bases, n);
    std::vector<std::vector<std::string>> cbs_v = c.get();  // combination bases vector
    for (size_t i = 0; i < cbs_v.size(); i++) { // 循环该位点每一种可能的碱基组合

        std::vector<double> obs_allele_freq = this->_set_initial_freq(cbs_v[i]);
        if (sum(obs_allele_freq) == 0) // Empty coverage for this type of combination, skip.
            throw std::runtime_error("The sum of frequence of active bases must always > 0. Check: " + 
                                     ngslib::join(cbs_v[i], ",") + " - " + ngslib::join(obs_allele_freq, ","));
        
        std::vector<double> log_marginal_likelihood;
        // The value of 'obs_allele_freq' and 'log_marginal_likelihood' will be updated in EM process.
        EM(_allele_likelihood, obs_allele_freq, log_marginal_likelihood);
        double sum_log_marginal_likelihood = sum(log_marginal_likelihood);

        data.bc.push_back(cbs_v[i]);
        data.bp.push_back(obs_allele_freq);
        data.lr.push_back(sum_log_marginal_likelihood);
    }
    
    return data;
}

void BaseType::lrt(const std::vector<std::string> &specific_bases) {

    if (_total_depth == 0) return;
    
    std::vector<std::string> active_bases;
    for (auto b: specific_bases) {
        // Get active bases which count frequence > _min_af
        if (_depth[b] / _total_depth >= _min_af)
            active_bases.push_back(b);
    }

    if (active_bases.size() == 0) return;

    // init. Base combination of active_bases
    AA var = _f(active_bases, active_bases.size());  // F4

    double chi_sqrt_value = 0;
    std::vector<double> active_bases_freq = var.bp[0];
    double lr = var.lr[0];  // f4

    // Find candinate altnative alleles
    for (size_t n = active_bases.size() - 1; n > 0; --n) {
        var = _f(active_bases, n);
        std::vector<double> lrt_chivalue;
        for (size_t j(0); j < var.lr.size(); ++j) {
            lrt_chivalue.push_back(2 * (lr - var.lr[j])); // 负号代入对数方程了，所以和公式比起来就是分子分母互换
        }
        size_t i_min = argmin(lrt_chivalue.begin(), lrt_chivalue.end());  // 选出最优组合

        lr = var.lr[i_min];
        chi_sqrt_value = lrt_chivalue[i_min]; // 获得该最优组合的卡方值

        // 注意这 H0 假设是“少碱基的组合与多碱基的组合相比无显著差异”，如果是，那么选H0，也就是少碱基的组合，否则取多碱基组合(H1)。
        // 这和我的计算公式是一致的，只是反过来而已
        if (chi_sqrt_value < LRT_THRESHOLD) { 
            // Take the null hypothesis and continue
            active_bases = var.bc[i_min];
            active_bases_freq = var.bp[i_min];
        } else {
            // Take the alternate hypothesis and break the loop
            break;
        }
    }

    // Set the active bases and their frequency
    this->_active_bases.clear();
    this->_af_by_lrt.clear();
    for (auto b: active_bases) {
        this->_active_bases.push_back(b);
        this->_af_by_lrt[b] = active_bases_freq[_B_IDX[b]];
    }

    // Todo: improve the calculation method for var_qual
    if (!this->_active_bases.empty()) {

        double r = this->_depth[active_bases[0]] / (double)(this->_total_depth);
        if ((active_bases.size() == 1) && (this->_total_depth > 10) && (r > 0.5)) {
            // Hard code for 'mono-allelelic' when depth > 10 and r > 0.5
            this->_var_qual = 10000;
        } else {
            // 'chi2_test' may return nan, which is caused by 'chi_sqrt_value' <= 0 and means p value is 1.0.
            double chi_prob = chi2_test(chi_sqrt_value, 1);  // Python: chi_prob = chi2.sf(chi_sqrt_value, 1)
            if (std::isnan(chi_prob)) 
                chi_prob = 1.0;

            this->_var_qual = (chi_prob) ? -10 * log10(chi_prob) : 10000;
            // _var_qual will been setted as -0.0 instand of 0.0 if it's 0, because of the phred-scale formular
            if (this->_var_qual == -0.0) this->_var_qual = 0.0;
        }
    }

    return;
}

