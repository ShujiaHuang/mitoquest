/**
 * @file basetype.h
 * @brief 
 * 
 * @author Shujia Huang
 * @date 2025-03-01
 * 
 */
#ifndef __INCLUDE_BASETYPE_H__
#define __INCLUDE_BASETYPE_H__

#include <iostream>
#include <cstdint> // uint32_t
#include <cmath>   // use exp() function
#include <string>
#include <vector>
#include <map>

#include "algorithm.h"
#include "io/utils.h"              // join()

#include "external/combinations.h"

static const std::vector<std::string> BASIC_BASES = {"A", "C", "G", "T"}; // 预定义这个值，限定 UNIQ_BASES 数组中至少有这四个碱基
static const int LRT_THRESHOLD  = 24;             // 24 corresponding to a chi-pvalue of 10^-6
static const int QUAL_THRESHOLD = 20;             // -10 * lg(10^-2)
static const double MLN10TO10   = -0.2302585093;  // -ln(10)/10 = -0.23025850929940458，换底，把 phred-value 换成 e 为底，方便调用 exp()

// A class for calculate the base probability
class BaseType {

private:
    /**
     * @brief Define a structure for recording allele information return by BaseType::_f() 
     * in this class.
     * 
     */
    struct AA {
        std::vector<std::vector<std::string>> bc;
        std::vector<std::vector<double>> bp;
        std::vector<double> lr;
    };

    std::vector<std::string> _UNIQ_BASES;
    std::map<std::string, size_t> _B_IDX;  // A map for recroding the base in _UNIQ_BASES => index

    std::string _ref_id;
    uint32_t _ref_pos;
    std::vector<std::string> _active_bases;           // the Ref and alternative bases
    std::map<std::string, std::string> _bases2ref;  // align/alt bases => ref_base
    double _var_qual;
    double _min_af;

    int _total_depth;                      // depth on ref_pos
    std::map<std::string, double> _depth;  // bases depth, double 是为了方便做除法

    // Estimated base frequency of _UNIQ_BASES by LRT
    std::map<std::string, double> _af_by_lrt;

    // _UNIQ_BASES likelihood vector for echo base string
    std::vector<std::vector<double>> _allele_likelihood; // 2d-array, n x _UNIQ_BASES.size() matrix, n is total_depth.
    
    // init the base likelihood by input bases
    std::vector<double> _set_initial_freq(const std::vector<std::string> &bases);

    /**
     * @brief Calculate population likelihood for all the combination of bases
     * 
     * @param bases A 1-d array. An array subset of bases from _UNIQ_BASES 
     * @param n     The combination number. n must less or equal to the length of ``bases``
     * 
     * @return AA   AA.bc: An array of combinamtion bases
     *              AA.lr: Likelihood of ``bc``
     * 
     */
    AA _f(const std::vector<std::string> &bases, int n);

public:
    struct BatchInfo {
        std::string ref_id;
        uint32_t    ref_pos;

        std::vector<std::string> ref_bases;
        std::vector<std::string> align_bases;
        std::vector<char>        align_base_quals;
        std::vector<int>         base_pos_ranks;

        std::vector<int>  mapqs;
        std::vector<char> map_strands;

        BatchInfo() : ref_id(""), ref_pos(0) {};
        BatchInfo(const std::string& rid, uint32_t pos) : ref_id(rid), ref_pos(pos) {};
    };

    // Constructor
    BaseType(){};
    BaseType(const BatchInfo *smp_bi, double af);
    BaseType(const BaseType &b);  // copy constructor

    ~BaseType(){};

    /**
     * @brief The main function for likelihood ratio test
     * 
     * @param specific_bases a 1d-array. [optional]
     * Calculating LRT for specific base combination if provided.
     */
    void lrt(const std::vector<std::string> &specific_bases);
    void lrt() { /* default */ lrt(_UNIQ_BASES); }

    const std::string &get_ref_id() const { return this->_ref_id; };
    const uint32_t &get_ref_pos() const { return this->_ref_pos; };
    const std::map<std::string, std::string> &get_bases2ref() const { return this->_bases2ref; };
    const std::vector<std::string> &get_active_bases() const { return this->_active_bases; }; // candidate variant bases

    const double get_var_qual() const { return this->_var_qual; }
    const int get_total_depth() const { return this->_total_depth; }
    const double get_base_depth(const std::string &b) const {
        // operator[] doesn't have a 'const' qualifier in std::map. Use 'at' instead in C++11
        // https://stackoverflow.com/questions/42095642/error-passing-const-stdmapint-int-as-this-argument-discards-qualifiers
        double d;
        try {
            d = this->_depth.at(b);
        } catch (const std::out_of_range &ex) {
            std::string what_ex(ex.what());
            throw std::runtime_error("[ERROR] out_of_range:: " + what_ex + " '"+ b + "' not found.");
        }
        return d;
    }

    const double get_lrt_af(const std::string &b) const { 
        double lrt_af;
        try {
            lrt_af = this->_af_by_lrt.at(b);
        } catch (const std::out_of_range &ex) {
            std::string what_ex(ex.what());
            throw std::runtime_error("[ERROR] out_of_range:: " + what_ex + " '"+ b + "' not found.");
        }
        return lrt_af; 
    }
}; // BaseType class

#endif
