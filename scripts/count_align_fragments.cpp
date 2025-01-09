/*
 * Author: Shujia Huang
 * Date: 2025-01-02
 *
 * g++ -o count_align_fragments count_align_fragments.cpp -lhts -lz -lcurl -pthread -O3
 * 
 * 
 **/
#include <iostream>
#include <getopt.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include <cmath>
#include <numeric>
#include <regex>

#include <iomanip>
#include <thread>
#include <mutex>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>

// Enumeration for sequencing type
enum class SeqType {
    AUTO,   // Automatically detect from BAM/CRAM flags
    PE,     // Paired-end sequencing
    SE      // Single-end sequencing
};

// Structure to hold statistics
struct Statistics {
    double mean;
    double confidence_interval_lower;
    double confidence_interval_upper;
    
    Statistics() : mean(0.0), confidence_interval_lower(0.0), confidence_interval_upper(0.0) {}
};

// Structure to hold chromosome-specific data
struct ChromosomeData {
    std::string name;
    int count;
    uint32_t length;
    double normalized_ratio;
    double gc_content;  // record GC content for whole chromosomes
    std::unordered_map<std::string, bool> processed_reads;
    Statistics CN_stats;  // record the chromosome copy number stats 
    
    ChromosomeData(const std::string& n, uint32_t len) 
        : name(n), count(0), length(len), normalized_ratio(0.0), gc_content(0.0) {}
};

// Function to check if a chromosome is autosomal
bool is_autosomal(const std::string& chrom_name) {
    // Match chr1-chr22 or just 1-22
    std::regex autosomal_regex("^(chr)?([1-9]|1[0-9]|2[0-2])$");
    return std::regex_match(chrom_name, autosomal_regex);
}

// Function to calculate mean
double calculate_mean(const std::vector<double>& values) {
    if (values.empty()) return 0.0;
    return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

// Function to calculate standard error of the mean
double calculate_standard_error(const std::vector<double>& values, double mean) {
    if (values.size() < 2) return 0.0;
    
    double sum_squared_diff = 0.0;
    for (double value : values) {
        sum_squared_diff += std::pow(value - mean, 2);
    }
    
    double variance = sum_squared_diff / (values.size() - 1);
    return std::sqrt(variance / values.size());
}

// Function to calculate statistics including 95% confidence interval
Statistics calculate_statistics(const std::vector<double>& values) {
    Statistics stats;
    if (values.empty()) return stats;

    // Calculate mean
    stats.mean = calculate_mean(values);

    // Calculate standard error
    double standard_error = calculate_standard_error(values, stats.mean);

    // Calculate 95% confidence interval (using t-distribution critical value 1.96)
    double confidence_interval = 1.96 * standard_error;
    stats.confidence_interval_lower = stats.mean - confidence_interval;
    stats.confidence_interval_upper = stats.mean + confidence_interval;

    return stats;
}

// Function to detect sequencing type from BAM/CRAM file
SeqType detect_seq_type(const char* input_path) {
    samFile* input = sam_open(input_path, "r");
    if (!input) {
        std::cerr << "Error: Could not open input file for type detection\n";
        return SeqType::PE; // Default to PE if cannot detect
    }

    bam_hdr_t* header = sam_hdr_read(input);
    if (!header) {
        sam_close(input);
        return SeqType::PE;
    }

    bam1_t* record = bam_init1();
    SeqType detected_type = SeqType::SE;
    int read_count = 0;
    const int READ_LIMIT = 1000; // Check first 1000 reads

    while (sam_read1(input, header, record) >= 0 && read_count < READ_LIMIT) {
        if (!(record->core.flag & BAM_FUNMAP)) {
            if (record->core.flag & BAM_FPAIRED) {
                detected_type = SeqType::PE;
                break;
            }
            read_count++;
        }
    }

    bam_destroy1(record);
    bam_hdr_destroy(header);
    sam_close(input);
    return detected_type;
}

// Mutex for thread-safe output
std::mutex output_mutex;

// Add this helper function for thread-safe output
void thread_safe_output(const std::string& message) {
    std::lock_guard<std::mutex> lock(output_mutex);
    std::cerr << message;
}

/*** Function to calculate the GC content ***/ 

// Structure for genomic interval
struct GenomicInterval {
    std::string chrom;
    int64_t start;  // 0-based start position
    int64_t end;    // 0-based end position (exclusive)
    
    GenomicInterval(const std::string& c, int64_t s, int64_t e) 
        : chrom(c), start(s), end(e) {}
};

// Function to calculate GC content for a specific interval
double calculate_interval_gc_content(faidx_t* fai, const GenomicInterval& interval) {
    if (!fai) {
        thread_safe_output("Error: Invalid FASTA index\n");
        return -1.0;
    }

    // Fetch sequence for the interval
    int seq_len;
    char* sequence = faidx_fetch_seq(fai, 
                                     interval.chrom.c_str(), 
                                     interval.start, 
                                     interval.end - 1,  // faidx_fetch_seq is inclusive of end
                                     &seq_len);
    
    if (!sequence) {
        thread_safe_output("Error: Could not fetch sequence for interval " + 
                           interval.chrom + ":" + 
                           std::to_string(interval.start) + "-" + 
                           std::to_string(interval.end) + "\n");
        return -1.0;
    }

    // Count GC bases
    int64_t gc_count = 0;
    int64_t total_count = 0;

    for (int i = 0; i < seq_len; i++) {
        char base = std::toupper(sequence[i]);
        if (base == 'G' || base == 'C') {
            gc_count++;
        }
        if (base != 'N') {  // Only count non-N bases in total
            total_count++;
        }
    }

    free(sequence);
    return total_count > 0 ? (double)gc_count / total_count : 0.0;
}

// Function to calculate GC content for whole chromosomes or specific intervals
class GCContentCalculator {
private:
    faidx_t* fai;
    std::mutex gc_mutex;

public:
    GCContentCalculator(const char* reference_path) {
        fai = fai_load(reference_path);
        if (!fai) {
            throw std::runtime_error("Could not load reference genome index");
        }
    }

    ~GCContentCalculator() {
        if (fai) {
            fai_destroy(fai);
        }
    }

    // Calculate GC content for a whole chromosome
    double calculate_chromosome_gc(const std::string& chrom_name) {
        int64_t length = faidx_seq_len(fai, chrom_name.c_str());
        if (length < 0) {
            thread_safe_output("Error: Could not get length for chromosome " + 
                             chrom_name + "\n");
            return -1.0;
        }
        
        GenomicInterval interval(chrom_name, 0, length);
        return calculate_interval_gc_content(fai, interval);
    }

    // Calculate GC content for a specific interval
    double calculate_interval_gc(const GenomicInterval& interval) {
        return calculate_interval_gc_content(fai, interval);
    }

    // Calculate GC content for multiple intervals in parallel
    void calculate_intervals_gc(std::vector<std::pair<GenomicInterval, double>>& intervals_gc) {
        std::vector<std::thread> threads;
        
        auto gc_worker = [this](std::pair<GenomicInterval, double>& interval_gc) {
            double gc = calculate_interval_gc_content(this->fai, interval_gc.first);
            
            std::lock_guard<std::mutex> lock(gc_mutex);
            interval_gc.second = gc;
        };

        // Process intervals in parallel
        for (auto& interval_gc : intervals_gc) {
            threads.emplace_back(gc_worker, std::ref(interval_gc));
            
            if (threads.size() >= std::thread::hardware_concurrency()) {
                for (auto& t : threads) {
                    t.join();
                }
                threads.clear();
            }
        }

        // Wait for remaining threads
        for (auto& t : threads) {
            t.join();
        }
    }

    // Parse interval string (e.g., "chr1:1000-2000")
    static GenomicInterval parse_interval(const std::string& interval_str) {
        std::regex interval_regex("([^:]+):([0-9]+)-([0-9]+)");
        std::smatch matches;
        
        if (!std::regex_match(interval_str, matches, interval_regex)) {
            throw std::runtime_error("Invalid interval format: " + interval_str);
        }
        
        std::string chrom = matches[1];
        int64_t start = std::stoll(matches[2]) - 1;  // Convert to 0-based
        int64_t end = std::stoll(matches[3]);
        
        return GenomicInterval(chrom, start, end);
    }
};

// Function to process reads for a specific chromosome
void process_chromosome(const char* input_path, 
                        const char* reference_path, 
                        ChromosomeData& chrom_data, 
                        int min_mapq,
                        const std::string& target_chrom, 
                        SeqType seq_type) {

    samFile* input = sam_open(input_path, "r");
    if (!input) {
        std::cerr << "Error: Could not open input file " << input_path << std::endl;
        return;
    }

    if (reference_path) {
        hts_set_fai_filename(input, reference_path);
    }

    bam_hdr_t* header = sam_hdr_read(input);
    if (!header) {
        std::cerr << "Error: Could not read header\n";
        return;
    }

    hts_idx_t* idx = sam_index_load(input, input_path);
    if (!idx) {
        std::cerr << "Error: Could not load index for " << input_path << std::endl;
        return;
    }

    int tid = bam_name2id(header, target_chrom.c_str());
    if (tid < 0) {
        std::cerr << "Warning: Chromosome " << target_chrom << " not found in BAM/CRAM file\n";
        return;
    }

    hts_itr_t* iter = sam_itr_queryi(idx, tid, 0, header->target_len[tid]);
    if (!iter) {
        std::cerr << "Error: Could not create iterator for " << target_chrom << std::endl;
        return;
    }

    bam1_t* record = bam_init1();
    while (sam_itr_next(input, iter, record) >= 0) {
        // Skip unmapped, secondary, and supplementary alignments
        if (record->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
            continue;
        }

        // Check mapping quality
        if (record->core.qual < min_mapq) {
            continue;
        }

        std::string read_name = bam_get_qname(record);
        if (seq_type == SeqType::SE) {
            // For single-end, count every read
            chrom_data.count++;
        } else {
            // For paired-end, count only the first read in pair
            if (record->core.flag & BAM_FREAD1) {
                if (chrom_data.processed_reads.find(read_name) == chrom_data.processed_reads.end()) {
                    chrom_data.count++;
                    chrom_data.processed_reads[read_name] = true;
                }
            }
        }
    }

    bam_destroy1(record);
    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    bam_hdr_destroy(header);
    sam_close(input);
}

void calculate_normalized_ratios(std::vector<ChromosomeData>& chromosomes) {
    int64_t total_fragments = 0;
    int64_t total_genome_length = 0;
    
    for (const auto& chrom : chromosomes) {
        total_fragments += chrom.count;
        total_genome_length += chrom.length;
    }

    for (auto& chrom : chromosomes) {
        double fragment_ratio  = static_cast<double>(chrom.count) / total_fragments;
        double length_ratio    = static_cast<double>(chrom.length) / total_genome_length;
        chrom.normalized_ratio = fragment_ratio / length_ratio;
    }

    // Calculate ratio of each chromosome's normalized ratio to each autosomal
    for (auto& chrom : chromosomes) {
        std::vector<double> ratio_ratios;  // Store ratios of normalized ratios
        for (auto& chrtmp : chromosomes) {
            // Ignore non Autosomal
            if (is_autosomal(chrtmp.name)) {
                ratio_ratios.push_back(chrom.normalized_ratio / chrtmp.normalized_ratio);
// std::cout << chrtmp.name << " : " << chrom.normalized_ratio << "\t" << chrtmp.normacatlized_ratio << "\t" << chrom.normalized_ratio / chrtmp.normalized_ratio << "\n";
            }
        }
        chrom.CN_stats = calculate_statistics(ratio_ratios);
    }
}

void print_usage() {
    std::cerr << "Usage: count_align_fragments [options] <input.bam/cram>\n"
              << "Options:\n"
              << "  -r, --reference   Reference genome file (required for CRAM)\n"
              << "  -q, --mapq        Minimum mapping quality score [0]\n"
              << "  -t, --threads     Number of threads [auto]\n"
              << "  -s, --seqtype     Sequencing type: auto|pe|se [auto]\n"
              << "  -h, --help        Print this help message\n";
    exit(1);
}

int main(int argc, char *argv[]) {
    // Default parameters
    int min_mapq = 0;
    const char* reference_path = nullptr;
    int num_threads  = std::thread::hardware_concurrency();
    SeqType seq_type = SeqType::AUTO;
    
    static struct option long_options[] = {
        {"reference", required_argument, 0, 'r'},
        {"mapq",      required_argument, 0, 'q'},
        {"threads",   required_argument, 0, 't'},
        {"seqtype",   required_argument, 0, 's'},
        {"help",      no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "r:q:t:s:h", long_options, &option_index)) != -1) {
        switch (c) {
            case 'r': reference_path = optarg;         break;
            case 'q': min_mapq = std::stoi(optarg);    break;
            case 't': num_threads = std::stoi(optarg); break;
            case 's':
                if (std::string(optarg)      ==   "pe") seq_type = SeqType::PE;
                else if (std::string(optarg) ==   "se") seq_type = SeqType::SE;
                else if (std::string(optarg) == "auto") seq_type = SeqType::AUTO;
                else {
                    std::cerr << "Error: Invalid sequencing type. Use 'auto', 'pe', or 'se'\n";
                    return 1;
                }
                break;
            case 'h': print_usage(); break;
            default:
                print_usage();
        }
    }

    if (optind >= argc) {
        std::cerr << "Error: Input file is required.\n";
        print_usage();
    }

    const char* input_path = argv[optind];

    // Auto-detect sequencing type if needed
    if (seq_type == SeqType::AUTO) {
        seq_type = detect_seq_type(input_path);
        std::cout << "#Auto-detected sequencing type: " 
                  << (seq_type == SeqType::PE ? "paired-end" : "single-end")
                  << std::endl;
    }

    // Open input file to get chromosome list
    samFile* input = sam_open(input_path, "r");
    if (!input) {
        std::cerr << "Error: Could not open input file " << input_path << std::endl;
        return 1;
    }

    bam_hdr_t* header = sam_hdr_read(input);
    if (!header) {
        std::cerr << "Error: Could not read header\n";
        return 1;
    }

    std::vector<ChromosomeData> chromosomes;
    for (int i = 0; i < header->n_targets; ++i) {
        chromosomes.emplace_back(header->target_name[i], header->target_len[i]);
    }

    bam_hdr_destroy(header);
    sam_close(input);

    // Process chromosomes in parallel
    std::vector<std::thread> threads;
    size_t chrom_index = 0;
    
    while (chrom_index < chromosomes.size()) {
        if (threads.size() < static_cast<size_t>(num_threads)) {
            threads.emplace_back(process_chromosome, input_path, reference_path,
                                 std::ref(chromosomes[chrom_index]), min_mapq,
                                 chromosomes[chrom_index].name, seq_type);
            chrom_index++;
        } else {
            for (auto& thread : threads) {
                if (thread.joinable()) {
                    thread.join();
                }
            }
            threads.clear();
        }
    }

    for (auto& thread : threads) {
        if (thread.joinable()) {
            thread.join();
        }
    }

    calculate_normalized_ratios(chromosomes);

    try {
        GCContentCalculator gc_calc(reference_path);

        // Calculate GC content for whole chromosomes
        for (auto& chrom : chromosomes) {
            chrom.gc_content = gc_calc.calculate_chromosome_gc(chrom.name);
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    // Output results
    std::cout << "#Sequencing type: " << (seq_type == SeqType::PE ? "paired-end" : "single-end") << "\n";
    std::cout << "#Chromosome\tFragments\tChrom_Length\tGC_Content\tFragment_Normalized_Ratio\t" 
              << "CopyNum\tCopyNum-CI95-Lower\tCopyNum-CI95-Upper\n";
    for (const auto& chrom : chromosomes) {
        std::cout << chrom.name   << "\t" 
                  << chrom.count  << "\t"
                  << chrom.length << "\t"
                  << std::fixed   << std::setprecision(4) << chrom.gc_content                         << "\t"
                  << std::fixed   << std::setprecision(4) << chrom.normalized_ratio                   << "\t"
                  << std::fixed   << std::setprecision(4) << chrom.CN_stats.mean                      << "\t"
                  << std::fixed   << std::setprecision(4) << chrom.CN_stats.confidence_interval_lower << "\t"
                  << std::fixed   << std::setprecision(4) << chrom.CN_stats.confidence_interval_upper << "\n"; 
    }

    return 0;
}
