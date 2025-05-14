#include <gtest/gtest.h>
#include "algorithm.h"

TEST(StatsTest, MeanCalculation) {
    std::vector<int> v{1, 2, 3, 4, 5};
    EXPECT_DOUBLE_EQ(mean(v), 3.0);
}

// TEST(FisherExactTest, BasicTest) {
//     ContingencyTable table(12, 5, 7, 10);
    
//     double p_value = fisher_exact_test(table, TestSide::TWO_SIDED);
//     EXPECT_NEAR(p_value, 0.1595, 0.0001);
// }

// TEST(FisherExactTest, InvalidInput) {
//     EXPECT_THROW(fisher_exact_test(-1, 2, 3, 4), std::invalid_argument);
//     EXPECT_THROW(ContingencyTable(-1, 2, 3, 4), std::invalid_argument);
// }
