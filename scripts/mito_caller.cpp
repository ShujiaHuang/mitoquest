/*
 * Call variants and heterophrasmy from mtDNA
 * Author: Shujia Huang
 * Date: 2025-01-03
 * 
 * g++ -o mito_caller mito_caller.cpp -lhts -lz -lpthread
 * 
 **/
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <getopt.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/vcf.h>

struct Options {
    std::string bam_file;
    std::string ref_file;
    std::string out_file;
    int min_base_qual = 20;
    int min_mapping_qual = 20;
    float min_vaf = 0.01;     // 最小变异频率
    int min_depth = 10;       // 最小测序深度
    float het_threshold = 0.95; // 异质性阈值
};

struct PileupData {
    int depth = 0;
    std::map<char, int> base_counts;
    std::vector<int> base_quals;
};

class MitoVariantCaller {
private:
    Options opts;
    samFile* bam = nullptr;
    faidx_t* fai = nullptr;
    bam_hdr_t* header = nullptr;
    htsFile* out_vcf = nullptr;
    bcf_hdr_t* vcf_header = nullptr;

    void initialize_vcf_header() {
        vcf_header = bcf_hdr_init("w");
        bcf_hdr_append(vcf_header, "##fileformat=VCFv4.2");
        bcf_hdr_append(vcf_header, "##source=MitoVariantCaller");
        bcf_hdr_append(vcf_header, "##reference=chrM");
        bcf_hdr_append(vcf_header, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
        bcf_hdr_append(vcf_header, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">");
        bcf_hdr_append(vcf_header, "##INFO=<ID=HET,Number=0,Type=Flag,Description=\"Heteroplasmy\">");
        bcf_hdr_append(vcf_header, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        bcf_hdr_append(vcf_header, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic Depths\">");
        bcf_hdr_append(vcf_header, "##FORMAT=<ID=VAF,Number=A,Type=Float,Description=\"Variant Allele Frequency\">");
        
        bcf_hdr_add_sample(vcf_header, header->target_name[0]);
        bcf_hdr_write(out_vcf, vcf_header);
    }

public:
    MitoVariantCaller(const Options& options) : opts(options) {
        // 打开输入文件
        bam = sam_open(opts.bam_file.c_str(), "r");
        header = sam_hdr_read(bam);
        fai = fai_load(opts.ref_file.c_str());
        
        // 打开输出VCF文件
        out_vcf = bcf_open(opts.out_file.c_str(), "w");
        initialize_vcf_header();
    }

    ~MitoVariantCaller() {
        if (bam) sam_close(bam);
        if (header) bam_hdr_destroy(header);
        if (fai) fai_destroy(fai);
        if (out_vcf) bcf_close(out_vcf);
        if (vcf_header) bcf_hdr_destroy(vcf_header);
    }

    // 计算碱基质量的似然值
    double calc_base_likelihood(int qual, char obs, char exp) {
        double error_prob = pow(10, -qual/10.0);
        if (obs == exp) return 1.0 - error_prob;
        return error_prob/3.0;
    }

    // 检测变异和异质性
    void detect_variants() {
        int tid = bam_name2id(header, "chrM");
        if (tid < 0) {
            std::cerr << "Error: chrM not found in BAM file" << std::endl;
            return;
        }

        int ref_len;
        char* ref = fai_fetch(fai, "chrM", &ref_len);
        
        // 设置pileup参数
        bam_plp_t plp = bam_plp_init(nullptr, nullptr);
        bam_plp_set_maxcnt(plp, 1000000);  // 设置最大覆盖度

        // 对每个位点进行遍历
        bam1_t* b = bam_init1();
        int pos;
        int tid_p;
        const bam_pileup1_t* pl;
        int n_plp;

        while ((pl = bam_plp_auto(plp, &tid_p, &pos, &n_plp)) != nullptr) {
            if (tid_p != tid) continue;
            
            PileupData pdata;
            
            // 统计每个碱基的数量和质量
            for (int i = 0; i < n_plp; ++i) {
                if (pl[i].is_del || pl[i].is_refskip) continue;
                if (pl[i].qpos < pl[i].b->core.l_qseq) {
                    int qual = bam_get_qual(pl[i].b)[pl[i].qpos];
                    if (qual >= opts.min_base_qual && pl[i].b->core.qual >= opts.min_mapping_qual) {
                        char base = seq_nt16_str[bam_seqi(bam_get_seq(pl[i].b), pl[i].qpos)];
                        pdata.base_counts[base]++;
                        pdata.base_quals.push_back(qual);
                        pdata.depth++;
                    }
                }
            }

            // 如果深度足够，进行变异检测
            if (pdata.depth >= opts.min_depth) {
                char ref_base = ref[pos];
                std::vector<std::pair<char, int>> variants;
                
                for (const auto& bc : pdata.base_counts) {
                    float vaf = (float)bc.second / pdata.depth;
                    if (bc.first != ref_base && vaf >= opts.min_vaf) {
                        variants.push_back(std::make_pair(bc.first, bc.second));
                    }
                }

                // 如果检测到变异，输出到VCF文件
                if (!variants.empty()) {
                    bcf1_t* rec = bcf_init();
                    rec->rid = tid;
                    rec->pos = pos;
                    
                    std::string ref_allele(1, ref_base);
                    std::string alt_allele(1, variants[0].first);
                    
                    bcf_update_alleles_str(vcf_header, rec, 
                        (ref_allele + "," + alt_allele).c_str());
                    
                    // 计算VAF
                    float vaf = (float)variants[0].second / pdata.depth;
                    bcf_update_info_float(vcf_header, rec, "AF", &vaf, 1);
                    
                    // 设置深度信息
                    bcf_update_info_int32(vcf_header, rec, "DP", &pdata.depth, 1);
                    
                    // 检测异质性
                    if (vaf < opts.het_threshold) {
                        bcf_update_info_flag(vcf_header, rec, "HET", nullptr, 1);
                    }
                    
                    // 设置基因型
                    std::vector<int32_t> ad = {pdata.base_counts[ref_base], 
                                             variants[0].second};
                    bcf_update_format_int32(vcf_header, rec, "AD", ad.data(), 2);
                    bcf_update_format_float(vcf_header, rec, "VAF", &vaf, 1);
                    
                    // 写入VCF记录
                    bcf_write(out_vcf, vcf_header, rec);
                    bcf_destroy(rec);
                }
            }
        }

        bam_plp_destroy(plp);
        bam_destroy1(b);
        free(ref);
    }
};

void print_usage() {
    std::cerr << "Usage: mito_caller [options] -b <bam_file> -r <ref_file> -o <out_file>\n"
              << "Options:\n"
              << "  -b FILE    Input BAM/CRAM file\n"
              << "  -r FILE    Reference FASTA file\n"
              << "  -o FILE    Output VCF file\n"
              << "  -q INT     Minimum base quality [20]\n"
              << "  -Q INT     Minimum mapping quality [20]\n"
              << "  -v FLOAT   Minimum variant allele frequency [0.01]\n"
              << "  -d INT     Minimum read depth [10]\n"
              << "  -t FLOAT   Heteroplasmy threshold [0.95]\n";
}

int main(int argc, char** argv) {
    Options opts;
    int c;
    while ((c = getopt(argc, argv, "b:r:o:q:Q:v:d:t:h")) >= 0) {
        switch (c) {
            case 'b': opts.bam_file = optarg; break;
            case 'r': opts.ref_file = optarg; break;
            case 'o': opts.out_file = optarg; break;
            case 'q': opts.min_base_qual = atoi(optarg); break;
            case 'Q': opts.min_mapping_qual = atoi(optarg); break;
            case 'v': opts.min_vaf = atof(optarg); break;
            case 'd': opts.min_depth = atoi(optarg); break;
            case 't': opts.het_threshold = atof(optarg); break;
            case 'h': print_usage(); return 0;
            default: print_usage(); return 1;
        }
    }

    if (opts.bam_file.empty() || opts.ref_file.empty() || opts.out_file.empty()) {
        print_usage();
        return 1;
    }

    try {
        MitoVariantCaller caller(opts);
        caller.detect_variants();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}