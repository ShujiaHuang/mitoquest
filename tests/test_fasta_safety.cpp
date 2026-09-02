#include <gtest/gtest.h>

#include "io/fasta.h"

TEST(Fasta, UnloadedMetadataQueriesAreSafe) {
    ngslib::Fasta fasta;

    EXPECT_FALSE(fasta.has_seq("chrM"));
    EXPECT_EQ(fasta.nseq(), 0);
    EXPECT_TRUE(fasta.iseq_name(0).empty());
    EXPECT_EQ(fasta.seq_length("chrM"), 0u);
    EXPECT_THROW(fasta.fetch("chrM"), std::invalid_argument);
}