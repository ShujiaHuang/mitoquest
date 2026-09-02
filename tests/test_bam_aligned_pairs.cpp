#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <string>

#include <htslib/sam.h>

#include "io/bam_record.h"

TEST(BamRecord, PaddingDoesNotConsumeQueryOrReference) {
    const std::string sam_path = "/tmp/mitoquest_bam_padding_test.sam";
    {
        std::ofstream sam(sam_path);
        ASSERT_TRUE(sam.is_open());
        sam << "@HD\tVN:1.6\n"
            << "@SQ\tSN:chrM\tLN:20\n"
            << "read1\t0\tchrM\t1\t60\t2M1P2M\t*\t0\t0\tACGT\tIIII\n";
    }

    samFile* file = sam_open(sam_path.c_str(), "r");
    ASSERT_NE(file, nullptr);
    sam_hdr_t* header = sam_hdr_read(file);
    ASSERT_NE(header, nullptr);

    ngslib::BamRecord record;
    ASSERT_GE(record.load_read(file, header), 0);
    const auto pairs = record.get_aligned_pairs("ACGTACGTACGTACGTACGT");

    ASSERT_EQ(pairs.size(), 4u);
    for (size_t index = 0; index < pairs.size(); ++index) {
        EXPECT_EQ(pairs[index].op, BAM_CMATCH);
        EXPECT_EQ(pairs[index].qpos, index);
        EXPECT_EQ(pairs[index].ref_pos, static_cast<hts_pos_t>(index));
    }

    sam_hdr_destroy(header);
    sam_close(file);
    std::remove(sam_path.c_str());
}

// Multi-base insertions must carry the *mean* quality of the whole inserted
// sequence in read_qual (v1.11.3 semantics), not just the first base's
// quality.  Qualities here: 40 (I), 30 (?), 40 (I) -> mean (30 + 40) / 2.
TEST(BamRecord, InsertionReadQualIsMeanOfInsertedSequence) {
    const std::string sam_path = "/tmp/mitoquest_bam_ins_qual_test.sam";
    {
        std::ofstream sam(sam_path);
        ASSERT_TRUE(sam.is_open());
        sam << "@HD\tVN:1.6\n"
            << "@SQ\tSN:chrM\tLN:20\n"
            << "read1\t0\tchrM\t1\t60\t5M2I10M\t*\t0\t0\tAAAAACCAAAAAAAAAA\tIIIII?IIIIIIIIIII\n";
    }

    samFile* file = sam_open(sam_path.c_str(), "r");
    ASSERT_NE(file, nullptr);
    sam_hdr_t* header = sam_hdr_read(file);
    ASSERT_NE(header, nullptr);

    ngslib::BamRecord record;
    ASSERT_GE(record.load_read(file, header), 0);
    const auto pairs = record.get_aligned_pairs(std::string(20, 'A'));

    // 5 M + 1 INS + 10 M
    ASSERT_EQ(pairs.size(), 16u);
    const auto& ins = pairs[5];
    EXPECT_EQ(ins.op, BAM_CINS);
    EXPECT_EQ(ins.read_base, 'C');
    EXPECT_EQ(ins.multi_base, "C");
    // mean quality (30 + 40) / 2 = 35 -> ASCII 35 + 33 = 68 ('D')
    EXPECT_EQ(static_cast<int>(static_cast<unsigned char>(ins.read_qual)), 68);

    sam_hdr_destroy(header);
    sam_close(file);
    std::remove(sam_path.c_str());
}

// A single-base insertion keeps the quality of that one base.
TEST(BamRecord, SingleBaseInsertionKeepsItsOwnQuality) {
    const std::string sam_path = "/tmp/mitoquest_bam_ins1_qual_test.sam";
    {
        std::ofstream sam(sam_path);
        ASSERT_TRUE(sam.is_open());
        sam << "@HD\tVN:1.6\n"
            << "@SQ\tSN:chrM\tLN:20\n"
            << "read1\t0\tchrM\t1\t60\t5M1I10M\t*\t0\t0\tAAAAACAAAAAAAAAA\tIIIII?IIIIIIIIII\n";
    }

    samFile* file = sam_open(sam_path.c_str(), "r");
    ASSERT_NE(file, nullptr);
    sam_hdr_t* header = sam_hdr_read(file);
    ASSERT_NE(header, nullptr);

    ngslib::BamRecord record;
    ASSERT_GE(record.load_read(file, header), 0);
    const auto pairs = record.get_aligned_pairs(std::string(20, 'A'));

    ASSERT_EQ(pairs.size(), 16u);
    const auto& ins = pairs[5];
    EXPECT_EQ(ins.op, BAM_CINS);
    EXPECT_EQ(ins.read_base, 'C');
    EXPECT_TRUE(ins.multi_base.empty());
    // single-base quality 30 -> ASCII 30 + 33 = 63 ('?')
    EXPECT_EQ(static_cast<int>(static_cast<unsigned char>(ins.read_qual)), 63);

    sam_hdr_destroy(header);
    sam_close(file);
    std::remove(sam_path.c_str());
}