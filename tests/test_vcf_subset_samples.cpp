#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "io/vcf.h"
#include "vcf_subset_samples.h"

namespace {

void write_vcf(const std::string& path, const std::string& contents) {
    std::ofstream output(path);
    ASSERT_TRUE(output.is_open());
    output << contents;
}

}  // namespace

TEST(VCFSubsetSamples, RecalculatesStandardInfoAfterKeepingNonLeadingAlt) {
    const std::string input = "/tmp/mitoquest_subsam_standard_info_input.vcf";
    const std::string output = "/tmp/mitoquest_subsam_standard_info_output.vcf";
    write_vcf(input, R"(##fileformat=VCFv4.2
##contig=<ID=chrM,length=100>
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Samples with data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AF,Number=R,Type=Float,Description="Allele fractions">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	first	second
chrM	1	.	A	C,G	.	PASS	AC=1,1;AN=2;AF=0.5,0.5;NS=2	GT:DP:AF	0:10:1,0.1,0.2	2:10:0.1,0.2,0.7
)");

    std::vector<std::string> arguments = {"subsam", "-i", input, "-o", output, "second"};
    std::vector<char*> argv;
    argv.reserve(arguments.size());
    for (std::string& argument : arguments) argv.push_back(argument.data());
    optind = 1;
    VCFSubsetSamples tool(static_cast<int>(argv.size()), argv.data());
    ASSERT_NO_THROW(tool.run());

    ngslib::VCFFile reader(output);
    const ngslib::VCFHeader header = reader.header();
    ngslib::VCFRecord record;
    ASSERT_EQ(reader.read(record), 0);
    ASSERT_EQ(record.unpack(BCF_UN_ALL), 0);

    const auto alts = record.alt();
    ASSERT_EQ(alts.size(), 1u);
    EXPECT_EQ(alts.front(), "G");
    std::vector<std::vector<int>> genotypes;
    ASSERT_EQ(record.get_genotypes(header, genotypes), 1);
    ASSERT_EQ(genotypes.size(), 1u);
    EXPECT_EQ(genotypes.front().front(), 1);

    std::vector<int32_t> ac;
    std::vector<int32_t> an;
    std::vector<int32_t> ns;
    std::vector<float> info_af;
    ASSERT_EQ(record.get_info_int(header, "AC", ac), 1);
    ASSERT_EQ(record.get_info_int(header, "AN", an), 1);
    ASSERT_EQ(record.get_info_float(header, "AF", info_af), 1);
    ASSERT_EQ(record.get_info_int(header, "NS", ns), 1);
    EXPECT_EQ(ac.front(), 1);
    EXPECT_EQ(an.front(), 1);
    EXPECT_FLOAT_EQ(info_af.front(), 1.0f);
    EXPECT_EQ(ns.front(), 1);

    std::vector<float> format_af;
    ASSERT_EQ(record.get_format_float(header, "AF", format_af), 2);
    ASSERT_EQ(format_af.size(), 2u);
    EXPECT_FLOAT_EQ(format_af[0], 0.1f);
    EXPECT_FLOAT_EQ(format_af[1], 0.7f);

    std::remove(input.c_str());
    std::remove(output.c_str());
}

// MitoQuest schema VCF：INFO 含 PT/REF_N/HET_N/HOM_N/VAF_MEAN 时，AN 保留
// MitoQuest 自定义"样本数"语义，同时镜像写入标准 NS tag。
TEST(VCFSubsetSamples, WritesNsAlongsideMitoquestAnAfterSubsetting) {
    const std::string input = "/tmp/mitoquest_subsam_ns_input.vcf";
    const std::string output = "/tmp/mitoquest_subsam_ns_output.vcf";
    write_vcf(input, R"(##fileformat=VCFv4.2
##contig=<ID=chrM,length=100>
##INFO=<ID=AN,Number=1,Type=Integer,Description="MitoQuest-specific sample count">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with non-missing GT">
##INFO=<ID=PT,Number=1,Type=String,Description="Plasmic type">
##INFO=<ID=REF_N,Number=1,Type=Integer,Description="Ref individuals">
##INFO=<ID=HET_N,Number=1,Type=Integer,Description="Het individuals">
##INFO=<ID=HOM_N,Number=1,Type=Integer,Description="Hom individuals">
##INFO=<ID=VAF_MEAN,Number=A,Type=Float,Description="Mean VAF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	first	second
chrM	1	.	A	C	.	PASS	AN=2;NS=2;PT=Mixed;REF_N=0;HET_N=1;HOM_N=1;VAF_MEAN=0.5	GT	0/1	1
)");

    std::vector<std::string> arguments = {"subsam", "-i", input, "-o", output, "second"};
    std::vector<char*> argv;
    argv.reserve(arguments.size());
    for (std::string& argument : arguments) argv.push_back(argument.data());
    optind = 1;
    VCFSubsetSamples tool(static_cast<int>(argv.size()), argv.data());
    ASSERT_NO_THROW(tool.run());

    ngslib::VCFFile reader(output);
    const ngslib::VCFHeader header = reader.header();
    ngslib::VCFRecord record;
    ASSERT_EQ(reader.read(record), 0);
    ASSERT_EQ(record.unpack(BCF_UN_ALL), 0);

    std::vector<int32_t> an;
    std::vector<int32_t> ns;
    ASSERT_EQ(record.get_info_int(header, "AN", an), 1);
    ASSERT_EQ(record.get_info_int(header, "NS", ns), 1);
    EXPECT_EQ(an.front(), 1);
    EXPECT_EQ(ns.front(), 1);

    std::remove(input.c_str());
    std::remove(output.c_str());
}

TEST(VCFSubsetSamples, DropsReferenceOnlySiteWhenInfoUpdateIsDisabled) {
    const std::string input = "/tmp/mitoquest_subsam_reference_only_input.vcf";
    const std::string output = "/tmp/mitoquest_subsam_reference_only_output.vcf";
    write_vcf(input, R"(##fileformat=VCFv4.2
##contig=<ID=chrM,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ref_sample	alt_sample
chrM	1	.	A	G	.	PASS	.	GT	0	1
)");

    std::vector<std::string> arguments = {
        "subsam", "-i", input, "-o", output, "--no-update-info", "ref_sample",
    };
    std::vector<char*> argv;
    argv.reserve(arguments.size());
    for (std::string& argument : arguments) argv.push_back(argument.data());
    optind = 1;
    VCFSubsetSamples tool(static_cast<int>(argv.size()), argv.data());
    ASSERT_NO_THROW(tool.run());

    ngslib::VCFFile reader(output);
    const ngslib::VCFHeader header = reader.header();
    ngslib::VCFRecord record;
    EXPECT_LT(reader.read(record), 0);
    EXPECT_EQ(header.n_samples(), 1);

    std::remove(input.c_str());
    std::remove(output.c_str());
}

// MitoQuest caller 的 AF/AQ 声明 Number=A 但按 GT 位置对齐（每样本值数
// == GT 条目数，可变）。bcf_trim_alleles 会因值数与 Number=A 基数不符
// 报错退出（mini_crash 布局），或在值数恰好等于基数时按标准布局误删
// 值；修复后应不崩溃、保留 GT 对齐布局、并归一化 REF/ALT。
// 输入带 caller 的 INFO schema（PT/REF_N/HET_N/HOM_N/VAF_MEAN）：值数
// 恰好等于基数时的 GT-aligned 歧义判定仅对 MitoQuest schema VCF 生效，
// 标准第三方 VCF 同布局走 htslib 路径（见 CleanupGenotypesStandard-
// DiploidWithAmbiguousCountsStaysHtslib）。
TEST(VCFSubsetSamples, KeepsGtAlignedFormatValuesAfterTrimmingUnusedAlt) {
    const std::string input = "/tmp/mitoquest_subsam_gt_aligned_input.vcf";
    const std::string output = "/tmp/mitoquest_subsam_gt_aligned_output.vcf";
    write_vcf(input, R"(##fileformat=VCFv4.2
##contig=<ID=chrM,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction in GT order">
##FORMAT=<ID=AQ,Number=A,Type=Integer,Description="Allele quality in GT order">
##INFO=<ID=PT,Number=1,Type=String,Description="Plasmic type">
##INFO=<ID=REF_N,Number=1,Type=Integer,Description="Ref count">
##INFO=<ID=HET_N,Number=1,Type=Integer,Description="Het count">
##INFO=<ID=HOM_N,Number=1,Type=Integer,Description="Hom count">
##INFO=<ID=VAF_MEAN,Number=A,Type=Float,Description="Mean VAF">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chrM	1	.	AC	CC,AT	.	PASS	.	GT:AF:AQ	0/1:0.9,0.1:30,25	0/2:0.8,0.2:31,26
)");

    // 抽 S1：ALT2 未用被删，GT 重映射，REF=AC/ALT=CC 应归一化为
    // REF=A/ALT=C（共同后缀修剪）；S1 的 GT-aligned AF/AQ 保持原样。
    std::vector<std::string> arguments = {"subsam", "-i", input, "-o", output,
                                          "--no-update-info", "S1"};
    std::vector<char*> argv;
    argv.reserve(arguments.size());
    for (std::string& argument : arguments) argv.push_back(argument.data());
    optind = 1;
    VCFSubsetSamples tool(static_cast<int>(argv.size()), argv.data());
    ASSERT_NO_THROW(tool.run());

    ngslib::VCFFile reader(output);
    const ngslib::VCFHeader header = reader.header();
    ngslib::VCFRecord record;
    ASSERT_EQ(reader.read(record), 0);
    ASSERT_EQ(record.unpack(BCF_UN_ALL), 0);

    EXPECT_EQ(record.ref(), "A");
    const auto alts = record.alt();
    ASSERT_EQ(alts.size(), 1u);
    EXPECT_EQ(alts.front(), "C");

    std::vector<std::vector<int>> genotypes;
    ASSERT_EQ(record.get_genotypes(header, genotypes), 2);
    ASSERT_EQ(genotypes.size(), 1u);
    EXPECT_EQ(genotypes.front().size(), 2u);
    EXPECT_EQ(genotypes.front()[0], 0);
    EXPECT_EQ(genotypes.front()[1], 1);

    std::vector<float> af;
    ASSERT_EQ(record.get_format_float(header, "AF", af), 2);
    ASSERT_EQ(af.size(), 2u);
    EXPECT_FLOAT_EQ(af[0], 0.9f);
    EXPECT_FLOAT_EQ(af[1], 0.1f);

    std::vector<int32_t> aq;
    ASSERT_EQ(record.get_format_int(header, "AQ", aq), 2);
    ASSERT_EQ(aq.size(), 2u);
    EXPECT_EQ(aq[0], 30);
    EXPECT_EQ(aq[1], 25);

    std::remove(input.c_str());
    std::remove(output.c_str());
}

// GT-aligned 值数可变（mini_crash 布局：S1 1 值、S2 2 值）+ 未用 ALT 清理：
// bcf_trim_alleles 严格校验报错退出，修复后应不崩溃。抽 S1（纯 REF）并
// 加 --keep-all-site 保留记录，验证 AF 单值完整输出。当所有 ALT 都未被
// 引用时（new_alts 为空），沿用 v1.11.3 语义：ALT 保留、GT 重映射。
TEST(VCFSubsetSamples, TrimsUnusedAltWithoutCrashingOnVariableGtAlignedValues) {
    const std::string input = "/tmp/mitoquest_subsam_var_gt_aligned_input.vcf";
    const std::string output = "/tmp/mitoquest_subsam_var_gt_aligned_output.vcf";
    write_vcf(input, R"(##fileformat=VCFv4.2
##contig=<ID=chrM,length=100>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction in GT order">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chrM	1	.	G	C	.	PASS	.	GT:AF	0:0.97	0/1:0.93,0.05
)");

    std::vector<std::string> arguments = {"subsam", "-i", input, "-o", output,
                                          "--no-update-info", "--keep-all-site", "S1"};
    std::vector<char*> argv;
    argv.reserve(arguments.size());
    for (std::string& argument : arguments) argv.push_back(argument.data());
    optind = 1;
    VCFSubsetSamples tool(static_cast<int>(argv.size()), argv.data());
    ASSERT_NO_THROW(tool.run());

    ngslib::VCFFile reader(output);
    const ngslib::VCFHeader header = reader.header();
    ngslib::VCFRecord record;
    ASSERT_EQ(reader.read(record), 0);
    ASSERT_EQ(record.unpack(BCF_UN_ALL), 0);

    EXPECT_EQ(record.ref(), "G");
    // v1.11.3 语义：所有 ALT 均未被引用时保留 ALT（new_alts 为空则跳过
    // update_alleles），仅重映射 GT；A/B 基准输出为 "G C" GT=0。
    EXPECT_EQ(record.n_alt(), 1);
    const auto alts = record.alt();
    ASSERT_EQ(alts.size(), 1u);
    EXPECT_EQ(alts.front(), "C");

    std::vector<std::vector<int>> genotypes;
    ASSERT_EQ(record.get_genotypes(header, genotypes), 1);
    ASSERT_EQ(genotypes.size(), 1u);
    EXPECT_EQ(genotypes.front().front(), 0);

    std::vector<float> af;
    ASSERT_EQ(record.get_format_float(header, "AF", af), 1);
    ASSERT_EQ(af.size(), 1u);
    EXPECT_FLOAT_EQ(af[0], 0.97f);

    std::remove(input.c_str());
    std::remove(output.c_str());
}