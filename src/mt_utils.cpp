#include "mt_utils.h"

std::string format_double(double value, int precision) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(precision) << value;
    return ss.str();
}

int get_total_depth(const AlignInfo &align_infor) {
    int depth(0);
    for (auto &ab: align_infor.align_bases) {
        depth += ab.read_base.size();
    }  

    return depth;
}

StrandBiasInfo strand_bias(const std::string &major_base, 
                           const std::string &alt_base,
                           const std::vector<std::string> &bases,
                           const std::vector<char> &strands)  // strands 和 bases 是配对的，一一对应 
{
    int maj_fwd = 0, maj_rev = 0;
    int alt_fwd = 0, alt_rev = 0;
    for (size_t i(0); i < bases.size(); ++i) {
        if (bases[i][0] == 'N' || bases[i][0] == 'n') continue;

        if (strands[i] == '+') {
            if (bases[i] == major_base) {
                ++maj_fwd;
            } else if (alt_base == bases[i]) {
                ++alt_fwd;
            }

        } else if (strands[i] == '-') {
            if (bases[i] == major_base) {
                ++maj_rev;
            } else if (alt_base == bases[i]) {
                ++alt_rev;
            }

        } else {
            throw std::runtime_error("[ERROR] Get strange strand symbol: " + bases[i] + " " + std::to_string(strands[i]));
        }
    }

    if (alt_base == major_base) { 
        // 如果 alt_base 刚好是 major_base 那么就按照 50% 的比例构造 alt_fwd 和 alt_rev 的理论值
        alt_fwd = alt_rev = int(std::round(0.5 * (maj_fwd + maj_rev)));
    }

    double fs = -10 * log10(fisher_exact_test(maj_fwd, maj_rev, alt_fwd, alt_rev));
    // 应对'全 major allele' 或者是 '全 alt allele'
    if (std::isinf(fs)) {
        fs = 10000;
    } else if (fs == 0) {
        fs = 0.0;
    }

    // Strand bias estimated by the Symmetric Odds Ratio test
    // https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_annotator_StrandOddsRatio.php
    double sor = (maj_rev * alt_fwd > 0) ? (double)(maj_fwd * alt_rev) / (double)(maj_rev * alt_fwd): 10000;

    StrandBiasInfo sbi;
    sbi.fwd = (alt_base == major_base) ? maj_fwd : alt_fwd;
    sbi.rev = (alt_base == major_base) ? maj_rev : alt_rev;
    sbi.fs  = fs; 
    sbi.sor = sor;

    return sbi;
}

AlleleInfo collect_allele_info(std::vector<VariantInfo>& variants) {
    // Find longest REF
    std::string shared_ref;
    for (const auto& smp_var : variants) {
        for (const auto& ref : smp_var.ref_bases) {
            if (ref.length() > shared_ref.length()) {
                shared_ref = ref;
            }
        }
    }
    AlleleInfo ai; 
    ai.ref = shared_ref;  // set raw ref and should be upper base already

    // Normalize ALTs and update variants
    std::set<std::string> unique_alts;
    for (auto& smp_var : variants) {  // Loop all samples in the position, collect and normalized the ALT information
        for (size_t j = 0; j < smp_var.alt_bases.size(); j++) {
            std::string alt = smp_var.alt_bases[j];
            std::string ref = smp_var.ref_bases[j];

            // rebase if Indels
            if (alt[0] == '-') {
                alt = ref[0];               // replace the bases by the first ref base
                smp_var.alt_bases[j] = alt; // rewrite deletion seq
            } else if (alt[0] == '+') {
                alt = ref + alt.substr(1);  // replace the first base('+') by ref bases
                smp_var.alt_bases[j] = alt; // rewrite insertion seq
            }
            
            if (ref != shared_ref && shared_ref.length() > ref.length()) {
                alt += shared_ref.substr(ref.length());
                smp_var.alt_bases[j] = alt;  // Update the ALT sequence
            }
            unique_alts.insert(alt);
        }
    }
    unique_alts.erase(shared_ref);  // Remove REF from ALTs

    // Set ALT field: Unique and sorted ALT (non-reference) sequences by length and then by ASCII
    ai.alts = ngslib::get_unique_strings(
        std::vector<std::string>(unique_alts.begin(), unique_alts.end())
    );

    for (const auto& alt : ai.alts) {
        ai.alt_all_freqs[alt] = std::vector<double>();
        ai.alt_het_freqs[alt] = std::vector<double>();
    }
    
    return ai;
}

VCFSampleAnnotation process_sample_variant(const VariantInfo& var_info,
                                           const std::vector<std::string>& ref_alt_order,
                                           double hf_cutoff) {
    // 计算并返回单个样本的信息
    VCFSampleAnnotation sa;
    
    int exp_major_count = (1 - hf_cutoff) * var_info.total_depth;
    int exp_minor_count = hf_cutoff * var_info.total_depth;
    
    for (size_t i = 0; i < ref_alt_order.size(); i++) {  // i == 0 represents the REF GT
        const auto& alt = ref_alt_order[i];
        for (size_t j = 0; j < var_info.alt_bases.size(); j++) {
            if (var_info.alt_bases[j] == alt) {
                sa.gtcode.push_back(i);
                sa.sample_alts.push_back(alt);
                sa.allele_depths.push_back(var_info.depths[j]);
                
                // double h = var_info.freqs[j]; // 这里不要用 lrt 计算出来的 allele frequency，因为可能不知为何会有负数（极少情况下）
                // calculate the allele frequency by allele_depth/total_depth
                double h = double(var_info.depths[j]) / double(var_info.total_depth);
                sa.allele_freqs.push_back(h);
                sa.logit_af.push_back(log(h/(1-h)));
                
                sa.ci_strings.push_back(format_double(var_info.ci[j].first, 4) + "," + 
                                        format_double(var_info.ci[j].second, 4));
                
                sa.sb_strings.push_back(std::to_string(var_info.strand_bias[j].fwd) + "," + 
                                        std::to_string(var_info.strand_bias[j].rev));
                
                sa.fs_strings.push_back(var_info.strand_bias[j].fs != 10000 ?
                                        format_double(var_info.strand_bias[j].fs, 3) : "10000"); // it's a phred-scale score
                sa.sor_strings.push_back(var_info.strand_bias[j].sor != 10000 ?
                                         format_double(var_info.strand_bias[j].sor, 3) : "10000");
                
                sa.var_types.push_back(var_info.var_types[j]);

                /**
                 * @brief determine if the rate of heteroplasmy is significantly greater than user defined cutoff.
                 * 
                 *  H0 (Null hypothesis): LESS
                 *  H1 (Alternative hypothesis): Equal or Greater
                 * 
                 *            major   minor
                 *  observed    n11     n12 | n1p
                 *  expected    n21     n22 | n2p
                 *          -----------------
                 *              np1     np2   npp
                 * 
                 *  where n11 and n12 are observed depth of major and minor alleles,
                 *  `n21 = (n11+n12) * (1 - hf_cutoff)` in which hf_cutoff is defined by `--het-threshold` 
                 *  `n22 = (n11+n12) * hf_cutoff` in which hf_cutoff is defined by `--het-threshold`
                 * 
                 */
                int obs_major_count = var_info.depths[var_info.major_allele_idx];
                int obs_minor_count = var_info.depths[j];
                double aq = -10 * log10(fisher_exact_test(obs_major_count, obs_minor_count,
                                                          exp_major_count, exp_minor_count,
                                                          TestSide::LESS));
                if (std::isinf(aq)) {
                    aq = 10000;
                }
                sa.aq.push_back(int(aq));
            }
        }
    }
    
    return sa;
}

std::string format_sample_string(const VCFSampleAnnotation& sa, const VariantInfo& var_info) {
    
    if (sa.sample_alts.empty()) {
        return ".:0:" + std::to_string(var_info.total_depth);  // No variants found
    }

    std::string sample_info = ngslib::join(sa.gtcode, "/")         + ":" +  // GT, genotype
                              std::to_string(int(var_info.qual))   + ":" +  // GQ, genotype quality
                              std::to_string(var_info.total_depth) + ":" +  // DP, total depth
                              ngslib::join(sa.allele_depths, ",")  + ":" +  // AD, active allele depth, so sum(AD) <= PD
                              ngslib::join(sa.allele_freqs, ",")   + ":" +  // AF, allele frequency
                              ngslib::join(sa.ci_strings, ";")     + ":" +  // CI, confidence interval
                              ngslib::join(sa.aq, ",")             + ":" +  // AQ, allele quality score
                              ngslib::join(sa.logit_af, ",")       + ":" +  // LAF,Transformed allele frequency by `logit`
                              ngslib::join(sa.sb_strings, ";")     + ":" +  // SB, strand bias
                              ngslib::join(sa.fs_strings, ",")     + ":" +  // FS, fisher-exact-test strand bias
                              ngslib::join(sa.sor_strings, ",")    + ":" +  // SOR,Strand odds ratio
                              ngslib::join(sa.var_types, ",");              // VT, Variant type

    return sample_info;
}

std::string vcf_header_define(const std::string &ref_file_path, const std::vector<std::string> &samples, const std::string other_comment) {
    std::vector<std::string> header = {
        "##fileformat=VCFv4.2",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth on the REF position\">",
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depth for each allele, in the order listed by GT\">",
        "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele fraction for each allele, in the order listed by GT. Fraction "
            "of non-reference allele corresponds to the variant allele fraction(VAF)\">",
        "##FORMAT=<ID=CI,Number=1,Type=String,Description=\"95\% confidence interval around the estimated allele fraction for "
            "the allele in the order listed by GT. format: ci_low,ci_up;ci_low,ci_up;...\">",
        "##FORMAT=<ID=AQ,Number=A,Type=Integer,Description=\"Allele quality, phred quality scores of pvalue of one-tail Fisher exact test "
            "to determine if the rate of allele is significantly greater than user defined cutoff (-j), in the order listed by GT. "
            "[CAUTION] In most cases, the minor allele corresponds to the heteroplasmic allele; therefore, the AQ at the minor allele position "
            "reflects the quality value of heterozygous allele mostly\">",
        "##FORMAT=<ID=LAF,Number=A,Type=Float,Description=\"Transformed AF: The logit of the allele fraction (AF) is "
            "computed as logit(AF) = ln(AF/(1-AF)) for each allele, in the order listed by GT\">",
        "##FORMAT=<ID=SB,Number=1,Type=String,Description=\"Allele-specific forward/reverse read counts for strand bias tests for the alleles, in "
            "the order listed by GT, separated by ';'. Format: fwd,rev;fwd,rev;...\">",
        "##FORMAT=<ID=FS,Number=A,Type=Float,Description=\"An ordered, comma delimited list of phred-scaled p-value using Fisher's exact test to detect strand bias\">",
        "##FORMAT=<ID=SOR,Number=A,Type=Float,Description=\"An ordered, comma delimited list of strand bias estimated by the Symmetric Odds Ratio test\">",
        "##FORMAT=<ID=VT,Number=1,Type=String,Description=\"An ordered, comma delimited list of variant type: REF, SNV, INS, DEL, or MNV\">",
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Number of samples with non-missing GT at this site\">",
        "##INFO=<ID=REF_N,Number=1,Type=Integer,Description=\"Total number of individuals exhibiting the reference state in the population\">",
        "##INFO=<ID=HET_N,Number=1,Type=Integer,Description=\"Total number of individuals exhibiting the heteroplasmic state in the population\">",
        "##INFO=<ID=HOM_N,Number=1,Type=Integer,Description=\"Total number of individuals exhibiting the homoplasmic state in the population\">",
        "##INFO=<ID=DP_MEAN,Number=1,Type=Float,Description=\"Mean mitochondrial sequencing depth across samples contributing to AN\">",
        "##INFO=<ID=DP_MEDIAN,Number=1,Type=Integer,Description=\"Median mitochondrial sequencing depth across samples contributing to AN\">",
        "##INFO=<ID=VAF_MEAN,Number=A,Type=Float,Description=\"Mean mitochondrial variant allele fraction(VAF) across all samples contributing to AN, "
            "with VAF=0 assigned to samples without detectable variant\">",
        "##INFO=<ID=VAF_MEDIAN,Number=A,Type=Float,Description=\"Median mitochondrial VAF across all samples contributing to AN\">",
        "##INFO=<ID=VAF_MEAN_HET,Number=A,Type=Float,Description=\"Mean mitochondrial VAF among heteroplasmic samples only\">",
        "##INFO=<ID=VAF_MEDIAN_HET,Number=A,Type=Float,Description=\"Median mitochondrial VAF among heteroplasmic samples only\">",
        "##INFO=<ID=PT,Number=1,Type=String,Description=\"Type of plasmicity observed in population: Ref, Hom, Het, or Mixed(Hom and Het)\">"


        // "##INFO=<ID=HOM_PF,Number=1,Type=Float,Description=\"Total frequency of individuals exhibiting the homoplasmic state for the non-reference allele in the population\">",
        // "##INFO=<ID=HET_PF,Number=1,Type=Float,Description=\"Total frequency of individuals exhibiting the heteroplasmic state for the non-reference allele in the population\">",
    };  // initial by common information of header

    ngslib::Fasta fa = ref_file_path;
    std::vector<std::string> contigs;
    for (size_t i(0); i < fa.nseq(); ++i) {
        std::string seqname = fa.iseq_name(i);
        uint32_t seqlen = fa.seq_length(seqname);
        contigs.push_back("##contig=<ID=" + seqname + ",length=" + std::to_string(seqlen) + ",assembly=" + ref_file_path + ">");
    }
    header.insert(header.end(), contigs.begin(), contigs.end());
    header.push_back("##reference=file://" + ngslib::abspath(ref_file_path));
    if (!other_comment.empty()) header.push_back(other_comment);
    header.push_back("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + ngslib::join(samples, "\t"));

    return ngslib::join(header, "\n");
}

void merge_file_by_line(const std::vector<std::string> &infiles, const std::string &outfile, 
                        std::string header, bool is_remove_tempfile) 
{
    if (infiles.empty()) return;

    bool is_compress = (ngslib::suffix_name(outfile) == ".gz") ? true : false;
    ngslib::BGZFile OUT(outfile, is_compress ? "wb" : "uw"); 
    OUT << header << "\n";

    /* Merge all files here */
    for (auto fn: infiles) {
        ngslib::BGZFile f(fn, "r");
        std::string line;

        while (f.readline(line)) {
            if (line[0] == '#') continue;
            OUT << line << "\n";
        }
        OUT.flush(); // 确保数据被写入

        if (is_remove_tempfile) ngslib::safe_remove(fn);
    }

    OUT.close();
    return;
}
