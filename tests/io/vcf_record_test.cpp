#include <gtest/gtest.h>
#include <cstdio>
#include <fstream>
#include "io/vcf.h"
#include "io/vcf_header.h"
#include "io/vcf_record.h"

TEST(VCFRecordTest, TestVariablePloidy) {
    // 创建一个临时的 VCF 文件用于测试
    const char* test_vcf_content = R"(##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=.,Type=String,Description="Genotype">
##contig=<ID=chrM,length=16569,assembly=chrM_rCRS.decoy.fa.gz>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2	Sample3
chrM	1	.	A	G	.	PASS	.	GT	0/1	1/1/1	1
chrM	2	.	C	T	.	PASS	.	GT	0/1/2	1	0/1/2/3
)";

    // 写入临时文件
    std::string temp_filename = "test_variable_ploidy.vcf";
    std::ofstream temp_file(temp_filename);
    temp_file << test_vcf_content;
    temp_file.close();

    try {
        // 打开 VCF 文件
        ngslib::VCFFile reader(temp_filename);
        ngslib::VCFHeader header = reader.header();
        ngslib::VCFRecord record;

        // 读取第一条记录并检查
        ASSERT_EQ(reader.read(record), 0);
        
        std::vector<std::vector<int>> genotypes;
        int max_ploidy = record.get_genotypes(header, genotypes);
        EXPECT_EQ(max_ploidy, 3);
        std::cout << "Record 1 ploidy: " << max_ploidy << std::endl;

        // 验证第一条记录的倍性和基因型
        EXPECT_EQ(genotypes.size(), 3);     // 3个样本
        EXPECT_EQ(genotypes[0].size(), 2);  // Sample1: 二倍体
        EXPECT_EQ(genotypes[1].size(), 3);  // Sample2: 三倍体
        EXPECT_EQ(genotypes[2].size(), 1);  // Sample3: 单倍体
        
        // 检查具体的基因型值
        EXPECT_EQ(genotypes[0][0], 0);  // Sample1: 0/1
        EXPECT_EQ(genotypes[0][1], 1);
        
        EXPECT_EQ(genotypes[1][0], 1);  // Sample2: 1/1/1
        EXPECT_EQ(genotypes[1][1], 1);
        EXPECT_EQ(genotypes[1][2], 1); 

        EXPECT_EQ(genotypes[2][0], 1);  // Sample3: 1

        // 读取第二条记录并检查
        ASSERT_EQ(reader.read(record), 0);
        genotypes.clear();
        max_ploidy = record.get_genotypes(header, genotypes);
        EXPECT_EQ(max_ploidy, 4); 
        std::cout << "Record 2 ploidy: " << max_ploidy << std::endl;

        // 验证第二条记录的倍性和基因型
        EXPECT_EQ(genotypes.size(), 3);
        EXPECT_EQ(genotypes[0].size(), 3);  // Sample1: 三倍体
        EXPECT_EQ(genotypes[1].size(), 1);  // Sample2: 单倍体
        EXPECT_EQ(genotypes[2].size(), 4);  // Sample3: 四倍体

        // 清理
        std::remove(temp_filename.c_str());
    } catch (const std::exception& e) {
        std::remove(temp_filename.c_str());
        FAIL() << "Exception: " << e.what();
    }
}

// 添加一个测试用例专门检查缺失值的处理
TEST(VCFRecordTest, TestMissingGenotypes) {
    const char* test_vcf_content = R"(##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=.,Type=String,Description="Genotype">
##contig=<ID=chrM,length=16569,assembly=chrM_rCRS.decoy.fa.gz>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2	Sample3
chrM	1	.	A	G	.	PASS	.	GT	./.	./././.	.
)";

    std::string temp_filename = "test_missing_genotypes.vcf";
    std::ofstream temp_file(temp_filename);
    temp_file << test_vcf_content;
    temp_file.close();

    try {
        ngslib::VCFFile reader(temp_filename);
        ngslib::VCFHeader header = reader.header();
        ngslib::VCFRecord record;

        ASSERT_EQ(reader.read(record), 0);
        
        std::vector<std::vector<int>> genotypes;
        int max_ploidy = record.get_genotypes(header, genotypes);

        // 验证缺失值的处理
        EXPECT_EQ(genotypes.size(), 3);
        
        // 检查每个样本的缺失值
        for (const auto& sample_gt : genotypes) {
            for (int gt : sample_gt) {
                EXPECT_EQ(gt, -1);  // 所有基因型都应该是缺失值
            }
        }

        // 清理
        std::remove(temp_filename.c_str());
    } catch (const std::exception& e) {
        std::remove(temp_filename.c_str());
        FAIL() << "Exception: " << e.what();
    }
}

TEST(VCFRecordTest, RequirePassRejectsUnfilteredRecords) {
    const std::string temp_filename = "test_filter_pass_semantics.vcf";
    {
        std::ofstream file(temp_filename);
        ASSERT_TRUE(file.is_open());
        file << "##fileformat=VCFv4.2\n"
             << "##contig=<ID=chrM,length=100>\n"
             << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
             << "chrM\t1\t.\tA\tG\t.\t.\t.\n"
             << "chrM\t2\t.\tA\tG\t.\tPASS\t.\n";
    }

    ngslib::VCFFile reader(temp_filename);
    const ngslib::VCFHeader header = reader.header();
    ngslib::VCFRecord record;

    ASSERT_EQ(reader.read(record), 0);
    ASSERT_EQ(record.unpack(BCF_UN_FLT), 0);
    EXPECT_FALSE(record.passed_filters(header));

    ASSERT_EQ(reader.read(record), 0);
    ASSERT_EQ(record.unpack(BCF_UN_FLT), 0);
    EXPECT_TRUE(record.passed_filters(header));

    std::remove(temp_filename.c_str());
}

TEST(VCFRecordTest, CleanupGenotypesStandardDiploidWithAmbiguousCountsStaysHtslib) {
    // Standard third-party VCF (no MitoQuest INFO schema): diploid samples,
    // n_alt = 2, FORMAT AF declared Number=A with 2 values per sample and
    // AD declared Number=R with 3 values.  The per-sample AF value count
    // (2) equals both the declared cardinality and the GT entry count,
    // which used to misclassify the record as GT-aligned.  After trimming
    // the unused ALT2 the self path left dangling values (2 AF / 3 AD),
    // violating the declared cardinality; the record must go through the
    // htslib path and shrink to the new cardinality instead.
    const std::string temp_filename = "test_trim_alleles_diploid_standard.vcf";
    {
        std::ofstream file(temp_filename);
        ASSERT_TRUE(file.is_open());
        file << "##fileformat=VCFv4.2\n"
             << "##contig=<ID=chrM,length=100>\n"
             << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
             << "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n"
             << "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depth\">\n"
             << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfirst\tsecond\n"
             << "chrM\t1\t.\tA\tC,G\t.\tPASS\t.\tGT:AF:AD\t"
             << "0/1:0.93,0.05:10,5,0\t0/1:0.80,0.10:8,2,0\n";
    }

    ngslib::VCFFile reader(temp_filename);
    const ngslib::VCFHeader original_header = reader.header();
    ngslib::VCFRecord record;
    ASSERT_EQ(reader.read(record), 0);

    ASSERT_TRUE(record.cleanup_genotypes(original_header));
    ASSERT_EQ(record.unpack(BCF_UN_ALL), 0);

    const auto alts = record.alt();
    ASSERT_EQ(alts.size(), 1u);
    EXPECT_EQ(alts.front(), "C");

    std::vector<std::vector<int>> genotypes;
    ASSERT_EQ(record.get_genotypes(original_header, genotypes), 2);
    ASSERT_EQ(genotypes.size(), 2u);
    for (const auto& gt : genotypes) {
        ASSERT_EQ(gt.size(), 2u);
        EXPECT_EQ(gt[0], 0);
        EXPECT_EQ(gt[1], 1);
    }

    // AF must shrink from 2 values/sample to 1; AD from 3 to 2.
    std::vector<float> format_af;
    std::vector<int32_t> ad;
    EXPECT_EQ(record.get_format_float(original_header, "AF", format_af), 1);
    EXPECT_EQ(record.get_format_int(original_header, "AD", ad), 2);
    ASSERT_EQ(format_af.size(), 2u);
    ASSERT_EQ(ad.size(), 4u);
    EXPECT_FLOAT_EQ(format_af[0], 0.93f);
    EXPECT_FLOAT_EQ(format_af[1], 0.80f);
    EXPECT_EQ(ad[0], 10);
    EXPECT_EQ(ad[1], 5);
    EXPECT_EQ(ad[2], 8);
    EXPECT_EQ(ad[3], 2);

    std::remove(temp_filename.c_str());
}

TEST(VCFRecordTest, CleanupGenotypesRemapsStandardAlleleIndexedFields) {
    const std::string temp_filename = "test_trim_alleles_remap.vcf";
    {
        std::ofstream file(temp_filename);
        ASSERT_TRUE(file.is_open());
        file << "##fileformat=VCFv4.2\n"
             << "##contig=<ID=chrM,length=100>\n"
             << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
             << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n"
             << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
             << "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depth\">\n"
             << "##FORMAT=<ID=FA,Number=A,Type=Float,Description=\"Allele frequency\">\n"
             << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfirst\tsecond\n"
             << "chrM\t1\t.\tA\tC,G\t.\tPASS\tAC=3,7;AF=0.3,0.7\tGT:AD:FA\t"
             << "1:10,8,2:0.2,0.9\t2:10,1,9:0.1,0.9\n";
    }

    ngslib::VCFFile reader(temp_filename);
    const ngslib::VCFHeader original_header = reader.header();
    ngslib::VCFRecord record;
    ASSERT_EQ(reader.read(record), 0);

    const ngslib::VCFHeader subset_header = original_header.subset_samples({"second"});
    std::vector<int> sample_indices = {1};
    ngslib::VCFRecord subset = record.subset_samples(original_header, sample_indices);
    ASSERT_TRUE(subset.cleanup_genotypes(subset_header));
    ASSERT_EQ(subset.unpack(BCF_UN_ALL), 0);

    const auto alts = subset.alt();
    ASSERT_EQ(alts.size(), 1u);
    EXPECT_EQ(alts.front(), "G");
    std::vector<std::vector<int>> genotypes;
    ASSERT_EQ(subset.get_genotypes(subset_header, genotypes), 1);
    ASSERT_EQ(genotypes.size(), 1u);
    ASSERT_EQ(genotypes.front().size(), 1u);
    EXPECT_EQ(genotypes.front().front(), 1);

    std::vector<int32_t> ac;
    std::vector<float> info_af;
    EXPECT_EQ(subset.get_info_int(subset_header, "AC", ac), 1);
    EXPECT_EQ(subset.get_info_float(subset_header, "AF", info_af), 1);
    ASSERT_EQ(ac.size(), 1u);
    ASSERT_EQ(info_af.size(), 1u);
    EXPECT_EQ(ac.front(), 7);
    EXPECT_FLOAT_EQ(info_af.front(), 0.7f);

    std::vector<int32_t> ad;
    std::vector<float> format_af;
    EXPECT_EQ(subset.get_format_int(subset_header, "AD", ad), 2);
    EXPECT_EQ(subset.get_format_float(subset_header, "FA", format_af), 1);
    ASSERT_EQ(ad.size(), 2u);
    ASSERT_EQ(format_af.size(), 1u);
    EXPECT_EQ(ad[0], 10);
    EXPECT_EQ(ad[1], 9);
    EXPECT_FLOAT_EQ(format_af.front(), 0.9f);

    std::remove(temp_filename.c_str());
}