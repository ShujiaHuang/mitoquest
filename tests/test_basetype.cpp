#include <gtest/gtest.h>

#include "basetype.h"

TEST(BaseType, IgnoresAmbiguousReadBases) {
    BaseType::BatchInfo batch("chrM", 100);
    batch.ref_bases = {"A", "A"};
    batch.align_bases = {"A", "R"};
    batch.align_base_quals = {'I', 'I'};

    BaseType base_type(&batch, 0.01);
    base_type.lrt();

    EXPECT_EQ(base_type.get_total_depth(), 1);
    ASSERT_EQ(base_type.get_active_bases().size(), 1u);
    EXPECT_EQ(base_type.get_active_bases().front(), "A");
}

TEST(BaseType, AllAmbiguousCoverageHasDeterministicZeroQuality) {
    BaseType::BatchInfo batch("chrM", 100);
    batch.ref_bases = {"A", "A"};
    batch.align_bases = {"N", "N"};
    batch.align_base_quals = {'I', 'I'};

    BaseType base_type(&batch, 0.01);
    base_type.lrt();

    EXPECT_EQ(base_type.get_total_depth(), 0);
    EXPECT_TRUE(base_type.get_active_bases().empty());
    EXPECT_DOUBLE_EQ(base_type.get_var_qual(), 0.0);
}