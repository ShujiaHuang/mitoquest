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
#include <ctime>

#include "version.h"
#include "mt_variant_caller.h"
#include "vcf_subset_samples.h"

static int usage() {
    std::cout << MITOQUEST_DESCRIPTION << "\n"
              << "Version: " << MITOQUEST_VERSION << "\n\n"
              << "Usage: mitoquest <command> [options]\n"
                 "Commands:\n"
                 "  caller    Mitochondrial variants and heteroplasmy/homoplasmy caller.\n"
                 "  subsam    Extract mitochondrial variants for specified samples from VCF files and output a new VCF file.\n"
              << "\n" << std::endl;
    return 1;
}

int main(int argc, char* argv[]) {
    time_t real_start_time = time(0);
    clock_t cpu_start_time = clock();

    if (argc < 2) {
        return usage();
    }

    // Save the complete command line with quotes for arguments containing spaces
    std::string cmdline = argv[0];  // start with program name
    for (int i = 1; i < argc; ++i) {
        cmdline += " " + std::string(argv[i]);
    }

    std::string cmd(argv[1]);
    if(cmd != "-h" && cmd != "--help") std::cout << "Commandline: " + cmdline + "\n" << std::endl;
    if (cmd == "caller") {
        try {
            MtVariantCaller caller(argc-1, argv+1);
            caller.run();
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << '\n';
            return 1;
        }

    } else if (cmd == "subsam") {
        try {
            // For subsampling, we need at least 3 arguments: input file, output file, 
            // and at least one sample name.
            VCFSubsetSamples subsam(argc-1, argv+1);
            subsam.run();

        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << '\n';
            return 1;
        }
    
    } else if (cmd == "-h" || cmd == "--help") {
        return usage();
        
    } else {
        std::cout << "\nError: Unrecognizable option: " + cmd << std::endl;
        return 1;
    }

    // Time information
    time_t now = time(0);
    std::string ct(ctime(&now));
    ct.pop_back();  // rm the trailing '\n' put by `asctime`
    std::cout << "\n[INFO] " + ct + ". Processes are all done, "
              << difftime(now, real_start_time) << " (CPU time: "
              << std::round((clock() - cpu_start_time) / CLOCKS_PER_SEC) 
              << ") seconds elapsed.\n" << std::endl;

    return 0;
}
