#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include <htslib/faidx.h>
#include <htslib/sam.h>

#include "io/iobgzf.h"
#include "mt_variant_caller.h"

namespace {

std::vector<std::string> split_fields(const std::string& line, char delimiter) {
    std::vector<std::string> fields;
    size_t begin = 0;
    while (begin <= line.size()) {
        const size_t end = line.find(delimiter, begin);
        fields.push_back(line.substr(begin, end - begin));
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return fields;
}

void write_indexed_bam(const std::string& sam_path, const std::string& bam_path) {
    samFile* sam = sam_open(sam_path.c_str(), "r");
    ASSERT_NE(sam, nullptr);
    sam_hdr_t* header = sam_hdr_read(sam);
    ASSERT_NE(header, nullptr);
    samFile* bam = sam_open(bam_path.c_str(), "wb");
    ASSERT_NE(bam, nullptr);
    ASSERT_EQ(sam_hdr_write(bam, header), 0);

    bam1_t* record = bam_init1();
    ASSERT_NE(record, nullptr);
    while (sam_read1(sam, header, record) >= 0) {
        ASSERT_GE(sam_write1(bam, header, record), 0);
    }

    bam_destroy1(record);
    sam_close(bam);
    sam_hdr_destroy(header);
    sam_close(sam);
    ASSERT_EQ(sam_index_build(bam_path.c_str(), 0), 0);
}

}  // namespace

TEST(MtVariantCaller, ExcludesSecondaryAndSupplementaryReadSupport) {
    const std::string prefix = "/tmp/mitoquest-caller-primary-only";
    const std::string fasta_path = prefix + ".fa";
    const std::string sam_path = prefix + ".sam";
    const std::string bam_path = prefix + ".bam";
    const std::string vcf_path = prefix + ".vcf.gz";

    {
        std::ofstream fasta(fasta_path);
        ASSERT_TRUE(fasta.is_open());
        fasta << ">chrM\n" << std::string(100, 'A') << "\n";
    }
    ASSERT_EQ(fai_build(fasta_path.c_str()), 0);

    {
        std::ofstream sam(sam_path);
        ASSERT_TRUE(sam.is_open());
        sam << "@HD\tVN:1.6\tSO:coordinate\n"
            << "@SQ\tSN:chrM\tLN:100\n"
            << "@RG\tID:rg1\tSM:S1\n";
        for (int index = 0; index < 11; ++index) {
            sam << "primary" << index << "\t0\tchrM\t1\t60\t10M\t*\t0\t0\t"
                << "CCCCCCCCCC\tIIIIIIIIII\tRG:Z:rg1\n";
        }
        sam << "secondary\t256\tchrM\t1\t60\t10M\t*\t0\t0\tCCCCCCCCCC\tIIIIIIIIII\tRG:Z:rg1\n"
            << "supplementary\t2048\tchrM\t1\t60\t10M\t*\t0\t0\tCCCCCCCCCC\tIIIIIIIIII\tRG:Z:rg1\n";
    }
    write_indexed_bam(sam_path, bam_path);

    std::vector<std::string> arguments = {
        "caller", "-f", fasta_path, "-o", vcf_path, "-Q", "0", "-q", "0",
        "-c", "100", "-t", "1", bam_path,
    };
    std::vector<char*> argv;
    argv.reserve(arguments.size());
    for (std::string& argument : arguments) argv.push_back(argument.data());
    optind = 1;
    MtVariantCaller caller(static_cast<int>(argv.size()), argv.data());
    caller.run();

    ngslib::BGZFile output(vcf_path, "r");
    std::string line;
    std::vector<std::string> first_record;
    while (output.readline(line)) {
        if (!line.empty() && line.front() != '#') {
            first_record = split_fields(line, '\t');
            break;
        }
    }
    ASSERT_GE(first_record.size(), 10u);
    EXPECT_EQ(first_record[6], ".");
    EXPECT_EQ(first_record[5], "10000");
    const auto sample_fields = split_fields(first_record[9], ':');
    ASSERT_GE(sample_fields.size(), 3u);
    EXPECT_EQ(sample_fields[2], "11");

    std::remove(fasta_path.c_str());
    std::remove((fasta_path + ".fai").c_str());
    std::remove(sam_path.c_str());
    std::remove(bam_path.c_str());
    std::remove((bam_path + ".bai").c_str());
    std::remove(vcf_path.c_str());
    std::remove((vcf_path + ".tbi").c_str());
}

TEST(MtVariantCaller, EmitsIndelAnchoredAtInternalChunkBoundary) {
    const std::string prefix = "/tmp/mitoquest-caller-indel-boundary";
    const std::string fasta_path = prefix + ".fa";
    const std::string sam_path = prefix + ".sam";
    const std::string bam_path = prefix + ".bam";
    const std::string vcf_path = prefix + ".vcf.gz";

    {
        std::ofstream fasta(fasta_path);
        ASSERT_TRUE(fasta.is_open());
        fasta << ">chrM\n" << std::string(200, 'A') << "\n";
    }
    ASSERT_EQ(fai_build(fasta_path.c_str()), 0);

    {
        std::ofstream sam(sam_path);
        ASSERT_TRUE(sam.is_open());
        sam << "@HD\tVN:1.6\tSO:coordinate\n"
            << "@SQ\tSN:chrM\tLN:200\n"
            << "@RG\tID:rg1\tSM:S1\n"
            << "boundary_ins\t0\tchrM\t91\t60\t10M1I10M\t*\t0\t0\t"
               "AAAAAAAAAACAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIII\tRG:Z:rg1\n";
    }
    write_indexed_bam(sam_path, bam_path);

    std::vector<std::string> arguments = {
        "caller", "-f", fasta_path, "-o", vcf_path, "-Q", "0", "-q", "0",
        "-c", "100", "-t", "1", bam_path,
    };
    std::vector<char*> argv;
    argv.reserve(arguments.size());
    for (std::string& argument : arguments) argv.push_back(argument.data());
    optind = 1;
    MtVariantCaller caller(static_cast<int>(argv.size()), argv.data());
    caller.run();

    ngslib::BGZFile output(vcf_path, "r");
    std::string line;
    bool found_boundary_insertion = false;
    while (output.readline(line)) {
        if (line.empty() || line.front() == '#') continue;
        const auto fields = split_fields(line, '\t');
        if (fields.size() >= 10 && fields[1] == "100") {
            found_boundary_insertion = true;
            EXPECT_EQ(fields[3], "A");
            EXPECT_EQ(fields[4], "AC");
            const auto sample_fields = split_fields(fields[9], ':');
            ASSERT_GE(sample_fields.size(), 3u);
            EXPECT_EQ(sample_fields[2], "1");
        }
    }
    EXPECT_TRUE(found_boundary_insertion);

    std::remove(fasta_path.c_str());
    std::remove((fasta_path + ".fai").c_str());
    std::remove(sam_path.c_str());
    std::remove(bam_path.c_str());
    std::remove((bam_path + ".bai").c_str());
    std::remove(vcf_path.c_str());
    std::remove((vcf_path + ".tbi").c_str());
}