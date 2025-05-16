#include <gtest/gtest.h>
#include <iostream>
#include <string>

#include "io/iobgzf.h"
#include "io/vcf.h"
#include "io/vcf_header.h"


// Example usage
TEST(VCFTest, BasicTest) {

    // 创建一个临时的 VCF 文件用于测试
    const char* test_vcf_content = R"(##fileformat=VCFv4.0
##FORMAT=<ID=GT,Number=.,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Reads covering the REF position">
##FORMAT=<ID=HF,Number=.,Type=Float,Description="Heteroplasmy Frequency of variant allele">
##FORMAT=<ID=CILOW,Number=.,Type=Float,Description="Value defining the lower limit of the confidence interval of the heteroplasmy fraction">
##FORMAT=<ID=CIUP,Number=.,Type=Float,Description="Value defining the upper limit of the confidence interval of the heteroplasmy fraction">
##FORMAT=<ID=SDP,Number=.,Type=String,Description="Strand-specific read depth of the ALT allele">
##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##contig=<ID=seq,length=1000,assembly=seq.fa.gz>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample
seq	601	.	T	C	.	PASS	AC=1;AN=1	GT:DP:HF:CILOW:CIUP:SDP	1:24:1.0:0.862:1.0:24;0
seq	714	.	AC	A	.	PASS	AC=1;AN=1	GT:DP:HF:CILOW:CIUP:SDP	1:24:1.0:0.862:1.0:24;0
)";
    // 写入临时文件
    std::string temp_filename = "test_vcffile.vcf.gz";
    ngslib::BGZFile OUT(temp_filename, "wb");
    OUT << test_vcf_content << "\n";
    OUT.close();

    ngslib::VCFFile reader(temp_filename);
    ngslib::VCFHeader hdr = reader.header(); // Get a copy we can query
    std::cout << reader.header() << std::endl;
    std::cout << "\n";
    std::cout << "Header: " << hdr << "\nSample num: " << hdr.n_samples() << std::endl;

    ngslib::VCFRecord rec;
    
    reader.read(rec);
    std::cout << "Record: " << rec << "\t" << rec.chrom(hdr) << std::endl;
    
    reader.read(rec);
    std::cout << "Record: " << rec << std::endl;
    
    // // Error
    // reader.read(rec);    
    // std::cout << "Record: " << rec << std::endl;

    std::remove(temp_filename.c_str());
}

