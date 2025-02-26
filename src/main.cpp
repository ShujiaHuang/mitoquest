/**
 * The main program entry point for the mtvariantcaller command-line tool.
 * 
 * @brief 
 * @file main.cpp
 * @author Shujia Huang (hshujia@qq.com)
 * @date 2025-02-11
 * 
 */
#include <getopt.h>
#include <iostream>
#include <cstdlib>

#include "mt_variant_caller.h"

void print_usage(const MtVariantCaller::Config &config) {
    std::cout << "Usage: hmtk [options] -f ref.fa -o output.vcf.gz in1.bam [in2.bam ...]\n"
              << "Required options:\n"
              << "  -f, --reference FILE       Reference FASTA file\n"
              << "  -o, --output FILE          Output VCF file\n\n"

              << "Optional options:\n"
              << "  -b, --bam-list FILE        list of input BAM/CRAM filenames, one per line.\n"
              << "  -r, --regions REG[,...]    Comma separated list of regions in which to process (default: entire genome).\n"
              << "                             REG format: chr:start-end (e.g.: chrM or chrM:1-1000,chrM:8000-8200)\n"
              << "  -q, --min-MQ INT           skip alignments with mapQ smaller than INT (default: "    << config.min_mapq  << ")\n"
              << "  -Q, --min-BQ INT           skip bases with base quality smaller than INT (default: " << config.min_baseq << ")\n"
              << "  -p, --pairs-map-only       Only use the paired reads which mapped to the some chromosome.\n"
              << "  -P, --proper-pairs-only    Only use properly paired reads.\n"
              << "  --filename-has-samplename  If the name of bamfile is something like 'SampleID.xxxx.bam', set this\n"
              << "                             argrument could save a lot of time during get the sample id from BAMfile.\n"
              << "  -j, --threshold FLOAT      Heteroplasmy threshold (default: " << config.heteroplasmy_threshold << ")\n"
              << "  -c, --chunk INT            Chunk size for parallel processing (default: " << config.chunk_size << ")\n"
              << "  -t, --threads INT          Number of threads (default: " << config.thread_count << ")\n"
              << "  -h, --help                 Print this help message.\n\n";
}

int main(int argc, char* argv[]) {
    
    MtVariantCaller::Config config;

    // Set default values
    config.min_mapq                = 0;
    config.min_baseq               = 20;
    config.pairs_map_only          = false;
    config.proper_pairs_only       = false;
    config.filename_has_samplename = false;
    config.heteroplasmy_threshold  = 0.2;
    config.thread_count            = 1;
    config.chunk_size              = 1000;

    static const struct option MT_CMDLINE_LOPTS[] = {
        {"reference",   required_argument, 0, 'R'},
        {"output",      required_argument, 0, 'o'},

        {"bam-list",    optional_argument, 0, 'b'},
        {"regions",     optional_argument, 0, 'r'},
        {"min-MQ",      optional_argument, 0, 'q'},
        {"min-BQ",      optional_argument, 0, 'Q'},
        {"threads",     optional_argument, 0, 't'},
        {"threshold",   optional_argument, 0, 'j'},
        {"chunk",       optional_argument, 0, 'c'},

        {"pairs-map-only",          no_argument, 0,    'p'}, // 小写 p
        {"proper-pairs-only",       no_argument, 0,    'P'},
        {"filename-has-samplename", no_argument, NULL, '1'},
        {"help",                    no_argument, 0,    'h'},

        // must set this value, to get the correct value from getopt_long
        {0, 0, 0, 0}
    };

    int opt;
    std::vector<std::string> bam_filelist;
    while ((opt = getopt_long(argc, argv, "f:b:o:r:q:Q:c:j:t:pPh", MT_CMDLINE_LOPTS, NULL)) != -1) {
        switch (opt) {
            case 'f': config.reference_file = optarg;                    break;
            case 'o': config.output_file    = optarg;                    break;
            case 'b': 
                bam_filelist = ngslib::get_firstcolumn_from_file(optarg);
                config.bam_files.insert(
                    config.bam_files.end(), 
                    bam_filelist.begin(), 
                    bam_filelist.end()
                );
                break;

            case 'r': config.calling_regions         = optarg;            break;
            case 'q': config.min_mapq                = std::atoi(optarg); break;
            case 'Q': config.min_baseq               = std::atoi(optarg); break;
            case 'p': config.pairs_map_only          = true;              break;
            case 'P': config.proper_pairs_only       = true;              break;
            case '1': config.filename_has_samplename = true;              break;
            case 'c': config.chunk_size              = std::atoi(optarg); break;
            case 'j': config.heteroplasmy_threshold  = std::atof(optarg); break;
            case 't': config.thread_count            = std::atoi(optarg); break;
            case 'h': print_usage(config); return 0;
            
            default:
                std::cerr << "Unknown argument: " << opt << std::endl;
                return 1;
        }
    }

    // Collect BAM/CRAM files
    while (optind < argc) {
        config.bam_files.push_back(argv[optind++]);
    }

    /* Make sure we set valid arguments */
    if (config.reference_file.empty() || config.bam_files.empty() || config.output_file.empty()) {
        std::cerr << "Error: Missing required arguments\n";
        print_usage(config);
        return 1;
    }

    if (config.min_mapq < 0) {
        std::cerr << "Error: Quality score must be non-negative\n";
        return 1;
    }

    if (config.heteroplasmy_threshold <= 0.0 || config.heteroplasmy_threshold > 1.0) {
        std::cerr << "Error: Heteroplasmy threshold must be between 0 and 1\n";
        return 1;
    }

    if (config.thread_count < 1) {
        std::cerr << "Error: Thread count must be at least 1\n";
        return 1;
    }

    if (config.chunk_size < 100) {
        std::cerr << "Error: Chunk size must be at least 100\n";
        return 1;
    }

    try {
        MtVariantCaller caller(config);
        if (!caller.run()) {
            std::cerr << "Error processing variants\n";
            return 1;
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
