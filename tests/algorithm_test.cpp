#include <gtest/gtest.h>
#include <iostream>
#include <cmath>
#include <iomanip>

#include <string>
#include <vector>
#include <array>
#include "algorithm.h"

TEST(AlgorithmTest, BasicTest) {

    double chi_sqrt_value = -0.1;
    double df = 1;
    double p_value = chi2_test(chi_sqrt_value, df);
    std::cout << "chi2_test(" << chi_sqrt_value << "," << df << ") = " 
              << p_value << " - " << std::isnan(p_value) << "\n";

    // Rank-sum-test
    std::vector<double> sample1 = {1, 5, 3, 10, 3, 3, 4, 5};
    std::vector<double> sample2 = {6, 7, 2, 2, 8, 9, 10};
    double ranksum_test_p = wilcoxon_ranksum_test(sample1, sample2);
    std::cout << "wilcoxon_ranksum_test: " << ranksum_test_p << std::endl;

    std::array<int, 7> numbers{ 2, 4, 8, 0, 6, -1, 3};
	int minIndex = argmin(numbers.begin(), numbers.end());
	std::cout << "{2, 4, 8, 0, 6, -1, 3} MinIndex: " << minIndex << '\n';
	std::vector<float> prices = {12.5f, 8.9f, 100.0f, 24.5f, 30.0f};
	float maxIndex = argmax(prices.begin(), prices.end());
	std::cout << "{12.5f, 8.9f, 100.0f, 24.5f, 30.0f} MaxIndex: " <<  maxIndex << '\n' << std::endl;

    // Test for calculating confidencel interval
    struct TestCase {
        int x;
        int n;
        double confidence;
        std::string description;
    };
    
    TestCase tests[] = {
        {20, 100, 0.95, "20% proportion, n=100, 95% confidence (Agresti-Coull)"},
        {5, 30, 0.95, "16.7% proportion, n=30, 95% confidence (Wilson)"},
        {0, 50, 0.95, "0% proportion, n=50, 95% confidence"},
        {50, 50, 0.99, "100% proportion, n=50, 99% confidence"},
        {45, 50, 0.90, "90% proportion, n=50, 90% confidence"}
    };
    
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Confidence Interval Test Examples\n";
    std::cout << "================================\n\n";
    for (const auto& test : tests) {
        try {
            auto [lower, upper] = calculate_confidence_interval(test.x, test.n, test.confidence);
            
            std::cout << test.description << "\n";
            std::cout << "Successes: " << test.x << ", Sample size: " << test.n << "\n";
            std::cout << (test.confidence * 100) << "% Confidence interval: [" 
                      << lower << ", " << upper << "] or [" 
                      << (lower * 100) << "%, " << (upper * 100) << "%]\n";
            std::cout << "Interval width: " << (upper - lower) << "\n\n";
        }
        catch (const std::exception& e) {
            std::cout << "Error in test case '" << test.description << "': " 
                      << e.what() << "\n\n";
        }
    }
}

TEST(StatsTest, MeanCalculation) {
    std::vector<int> v{1, 2, 3, 4, 5};
    EXPECT_DOUBLE_EQ(mean(v), 3.0);
}

TEST(StatsTest, SingleElementTest) {
    std::vector<int> single{5};
    EXPECT_DOUBLE_EQ(mean(single), 5.0);
}

TEST(StatsTest, NegativeNumbersTest) {
    std::vector<int> negative{-1, -2, -3, -4, -5};
    EXPECT_DOUBLE_EQ(mean(negative), -3.0);
}

TEST(FisherExactTest, BasicTest) {

    try {
        // Using direct parameters
        double p_value = fisher_exact_test(12, 5, 7, 10, TestSide::TWO_SIDED);
        std::cout << "fisher_exact_test(12, 5, 7, 10, TestSide::TWO_SIDED): " << p_value << "\n";
        EXPECT_NEAR(p_value, 0.17, 0.16); // Adjust the expected value based on your test case
        
        // Using ContingencyTable
        ContingencyTable table(12, 5, 7, 10);
        double p2 = fisher_exact_test(table, TestSide::LESS);
        std::cout << "fisher_exact_test(table(12, 5, 7, 10), TestSide::LEFT_SIDED): " << p2 << "\n";

        double p3 = fisher_exact_test(table, TestSide::GREATER);
        std::cout << "fisher_exact_test(table(12, 5, 7, 10), TestSide::RIGHT_SIDED): " << p3 << "\n";

        // more Fisher exact test
        double p = fisher_exact_test(345, 455, 260, 345);
        std::cout << "Fisher exact test: " << p << "\n";
        std::cout << "Fisher exact test (8, 4, 4, 9) : " << fisher_exact_test(8, 4, 4, 9) << "\n";
        std::cout << "Fisher exact test (4, 1, 2, 3) : " << fisher_exact_test(4, 1, 2, 3) << "\n";
        std::cout << "Fisher exact test (3, 3, 2, 2) : " << fisher_exact_test(3, 3, 2, 2) << "\n";
        std::cout << "Fisher exact test (2, 4, 3, 1) : " << fisher_exact_test(2, 4, 3, 1) << "\n";
        std::cout << "Fisher exact test (5, 0, 1, 4) : " << fisher_exact_test(5, 0, 1, 4) << "\n";
        std::cout << "Fisher exact test (3, 0, 0, 3) : " << fisher_exact_test(3, 0, 0, 3) << "\n";

        std::cout << "Fisher exact test (89,353,41,0) : " << fisher_exact_test(89,353,41,0) << "\n";
        std::cout << "Fisher exact test (89,353,221,221) : " << fisher_exact_test(89,353,221,221) << "\n";
        std::cout << "Fisher exact test (221,221,89,353) : " << fisher_exact_test(89,353,221,221) << "\n";
        std::cout << "Fisher exact test (527,391,459,459) : " << fisher_exact_test(527,391,459,459) << "\n";
        
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

TEST(FisherExactTest, SymmetricTest) {
    // 测试对称性：交换行或列不应改变 p 值
    ContingencyTable table1(10, 5, 8, 4);
    ContingencyTable table2(5, 10, 4, 8);
    
    double p1 = fisher_exact_test(table1, TestSide::TWO_SIDED);
    double p2 = fisher_exact_test(table2, TestSide::TWO_SIDED);
    
    EXPECT_NEAR(p1, p2, 1e-10);
}

TEST(FisherExactTest, ExtremeCaseTest) {
    // 测试极端情况
    ContingencyTable table(100, 0, 0, 100);
    double p_value = fisher_exact_test(table, TestSide::TWO_SIDED);
    std::cout << "fisher_exact_test(100, 0, 0, 100, TestSide::TWO_SIDED): " << p_value << "\n";
    EXPECT_LT(p_value, 0.0001);  // 期望 p 值非常小
}