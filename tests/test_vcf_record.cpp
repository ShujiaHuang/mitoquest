#include <gtest/gtest.h>
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