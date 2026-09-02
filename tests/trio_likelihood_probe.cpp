#include <algorithm>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>

#include "ne_estimate.h"

int main(int argc, char* argv[]) {
    if (argc != 8) {
        std::cerr << "usage: mitoquest_trio_probe "
                  << "<g_dp> <g_ad_alt> <m_dp> <m_ad_alt> "
                  << "<c_dp> <c_ad_alt> <Ne>\n";
        return EXIT_FAILURE;
    }

    try {
        NeEstimator::PairData pair;
        pair.g_dp = std::stoi(argv[1]);
        pair.g_ad_alt = std::stoi(argv[2]);
        pair.m_dp = std::stoi(argv[3]);
        pair.m_ad_alt = std::stoi(argv[4]);
        pair.c_dp = std::stoi(argv[5]);
        pair.c_ad_alt = std::stoi(argv[6]);
        pair.has_g = 1;

        const double ne = std::stod(argv[7]);
        const int cache_size = std::max(pair.m_dp, pair.c_dp);
        NeEstimator::LogFactorial log_factorial(cache_size);
        const double log_likelihood = NeEstimator::compute_ll_trio_continuous(pair, ne, log_factorial);
        std::cout << std::setprecision(17) << log_likelihood << "\n";

        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "mitoquest_trio_probe: " << error.what() << "\n";
        return EXIT_FAILURE;
    }
}