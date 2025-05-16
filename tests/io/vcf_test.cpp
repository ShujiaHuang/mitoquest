#include <iostream>
#include "io/vcf.h"
#include "io/vcf_header.h"


// Example usage
void example() {
    ngslib::VCFFile reader("data/test.vcf.gz");
    // ngslib::VCFFile reader("data/ttt.gz");
    ngslib::VCFHeader hdr = reader.header(); // Get a copy we can query
    std::cout << reader.header() << std::endl;
    std::cout << "\n";
    std::cout << "Header: " << hdr << "\nSample num: " << hdr.n_samples() << std::endl;

    ngslib::VCFRecord rec;
    reader.read(rec);
    std::cout << "Record: " << rec << "\t" << rec.chrom(hdr) << std::endl;
    reader.read(rec);
    std::cout << "Record: " << rec << std::endl;
    reader.read(rec);
    std::cout << "Record: " << rec << std::endl;
}

int main(int argc, char* argv[]) {
    example();
    return 0;
}

