#include <iostream>
#include "io/hts_utils.h"

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1;
    }

    const std::string filename(argv[1]);
    ngslib::Format format = ngslib::get_format(filename);
    
    std::cout << "File: " << filename << std::endl;
    std::cout << "Format: ";
    
    switch (format) {
        case ngslib::Format::bam:
            std::cout << "BAM";
            break;
        case ngslib::Format::cram:
            std::cout << "CRAM";
            break;
        case ngslib::Format::sam:
            std::cout << "SAM";
            break;
        case ngslib::Format::vcf:
            std::cout << "VCF";
            break;
        case ngslib::Format::bcf:
            std::cout << "BCF";
            break;
        case ngslib::Format::text:
            std::cout << "text";
            break;
        case ngslib::Format::unknown:
            std::cout << "Unknown";
            break;
        default:
            std::cout << "Other";
    }
    std::cout << std::endl;

    if (ngslib::is_cram(filename)) {
        std::cout << "This is a CRAM file" << std::endl;
    }

    return 0;
}

