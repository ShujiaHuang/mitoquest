#include "mt_caller_utils.h"

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
            throw std::runtime_error("[ERROR] Get strange strand symbol: " + std::to_string(strands[i]));
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
    ai.ref = shared_ref;  // set raw ref
    std::transform(shared_ref.begin(), shared_ref.end(), shared_ref.begin(), ::toupper);

    // Normalize ALTs and update variants
    std::set<std::string> unique_alts;
    for (auto& smp_var : variants) {  // Loop all samples in the position, collect and normalized the ALT information
        for (size_t j = 0; j < smp_var.alt_bases.size(); j++) {
            std::string alt = smp_var.alt_bases[j];
            std::string ref = smp_var.ref_bases[j];
            std::transform(ref.begin(), ref.end(), ref.begin(), ::toupper);

            // rebase if Indels
            if (alt[0] == '-') {
                alt = ref[0];               // replace by the first ref bases for DEL seq
                smp_var.alt_bases[j] = alt; // Update the ALT sequence
            } else if (alt[0] == '+') {
                alt = ref + alt.substr(1);  // replace the first base('+') by ref bases
                smp_var.alt_bases[j] = alt; // rewrite deletion seq
            }
            
            if (ref != shared_ref && shared_ref.length() > ref.length()) {
                alt += shared_ref.substr(ref.length());
                smp_var.alt_bases[j] = alt;  // Update the ALT sequence
            }
            unique_alts.insert(alt);
        }
    }
    unique_alts.erase(shared_ref);  // Remove REF from ALTs

    // Set ALT field: Unique and sorted ALT sequences by length and then by ASCII
    ai.alts = ngslib::get_unique_strings(
        std::vector<std::string>(unique_alts.begin(), unique_alts.end())
    );

    for (const auto& alt : ai.alts) {
        ai.allele_counts[alt] = 0;
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
                sa.gt_indices.push_back(i);
                sa.sample_alts.push_back(alt);
                sa.allele_depths.push_back(var_info.depths[j]);
                
                // allele_freqs.push_back(var_info.freqs[j]); // 这里不要用 lrt 计算出来的 allele frequency，因为可能不知为何会有负数（极少情况下）
                 // calculate the allele frequency by allele_depth/total_depth
                double h = double(var_info.depths[j]) / double(var_info.total_depth);
                sa.allele_freqs.push_back(h);
                sa.logit_hf.push_back(log(h/(1-h)));
                
                sa.ci_strings.push_back(format_double(var_info.ci[j].first) + "," + 
                                        format_double(var_info.ci[j].second));
                
                sa.sb_strings.push_back(std::to_string(var_info.strand_bias[j].fwd) + "," + 
                                        std::to_string(var_info.strand_bias[j].rev));
                
                sa.fs_strings.push_back(var_info.strand_bias[j].fs != 10000 ?
                                        format_double(var_info.strand_bias[j].fs) : "10000"); // it's a phred-scale score
                
                sa.sor_strings.push_back(var_info.strand_bias[j].sor != 10000 ?
                                         format_double(var_info.strand_bias[j].sor) : "10000");
                
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
                double hq = -10 * log10(fisher_exact_test(obs_major_count, obs_minor_count,
                                                          exp_major_count, exp_minor_count,
                                                          TestSide::LESS));
                if (std::isinf(hq)) {
                    hq = 10000;
                }
                sa.hq.push_back(int(hq));
            }
        }
    }
    
    return sa;
}

std::string format_sample_string(const VCFSampleAnnotation& sa, const VariantInfo& var_info) {
    
    if (sa.sample_alts.empty()) {
        return ".:0:" + std::to_string(var_info.total_depth);  // No variants found
    }

    std::string sample_info = ngslib::join(sa.gt_indices, "/")     + ":" +  // GT, genotype
                              std::to_string(int(var_info.qual))   + ":" +  // GQ, genotype quality (Variant quality)
                              std::to_string(var_info.total_depth) + ":" +  // DP, total depth
                              ngslib::join(sa.allele_depths, ",")  + ":" +  // AD, active allele depth, so sum(AD) <= PD
                              ngslib::join(sa.allele_freqs, ",")   + ":" +  // HF, allele frequency, homo-/hetero-plasmy
                              ngslib::join(sa.ci_strings, ";")     + ":" +  // CI, confidence interval
                              ngslib::join(sa.hq, ",")             + ":" +  // HQ, homo-/hetero-plasmy quality score
                              ngslib::join(sa.logit_hf, ",")       + ":" +  // LHF,Transformed homo-/hetero-plasmy by `logit`
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
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">",
        "##FORMAT=<ID=HF,Number=A,Type=Float,Description=\"Heteroplasmy/Homoplasmy fraction for the ref and alt alleles in the order listed\">",
        "##FORMAT=<ID=CI,Number=1,Type=String,Description=\"95\% confidence interval around the estimated homoplasmy/heteroplasmy fraction for "
        "the GT alleles in the order listed. format: ci_low,ci_up;ci_low,ci_up;...\">",
        "##FORMAT=<ID=HQ,Number=A,Type=Integer,Description=\"Heteroplasmy/Homoplasmy Quality, phred quality scores of pvalue of one-tail Fisher exact test "
        "to determine if the rate of allele is significantly greater than user defined cutoff (-j), in the order listed by GT. "
        "[CAUTION] In most cases, the minor allele corresponds to the heteroplasmic allele; therefore, the HQ at the minor allele position "
        "reflects the quality value of heterozygous allele mostly.\">",
        "##FORMAT=<ID=LHF,Number=A,Type=Float,Description=\"Transformed heteroplasmy/homoplasmy: The logit of the heteroplasmy/homoplasmy fraction (HF) is "
        "computed as logit(HF) = ln(HF/(1-HF)) for each allele, in the order listed by GT.\">",
        "##FORMAT=<ID=SB,Number=1,Type=String,Description=\"Allele-specific forward/reverse read counts for strand bias tests for the alleles, in "
        "the order listed by GT, separated by ';'. Format: fwd,rev;fwd,rev;...\">",
        "##FORMAT=<ID=FS,Number=A,Type=Float,Description=\"An ordered, comma delimited list of phred-scaled p-value using Fisher's exact test to detect strand bias\">",
        "##FORMAT=<ID=SOR,Number=A,Type=Float,Description=\"An ordered, comma delimited list of strand bias estimated by the Symmetric Odds Ratio test\">",
        "##FORMAT=<ID=VT,Number=1,Type=String,Description=\"An ordered, comma delimited list of variant type: REF, SNV, INS, DEL, or MNV\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"An ordered, comma delimited list of non-reference allele frequencies\">",
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"An ordered, comma delimited list of non-reference allele count in genotypes\">",
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of ref and non-ref allele in called genotypes\">",
        "##INFO=<ID=HOM_N,Number=1,Type=Integer,Description=\"Total number of individuals exhibiting the homoplasmic state for the non-reference allele in the population.\">",
        "##INFO=<ID=HET_N,Number=1,Type=Integer,Description=\"Total number of individuals exhibiting the heteroplasmic state for the non-reference allele in the population.\">",
        "##INFO=<ID=Total_N,Number=1,Type=Integer,Description=\"Available sample size in this record.\">",
        "##INFO=<ID=HOM_AF,Number=1,Type=Float,Description=\"Total frequency of individuals exhibiting the homoplasmic state for the non-reference allele in the population.\">",
        "##INFO=<ID=HET_AF,Number=1,Type=Float,Description=\"Total frequency of individuals exhibiting the heteroplasmic state for the non-reference allele in the population.\">",
        "##INFO=<ID=SUM_AF,Number=1,Type=Float,Description=\"The frequency of HOM_AF+HET_AF.\">",
        "##INFO=<ID=PT,Number=1,Type=String,Description=\"Type of plasmicity observed in population: Hom_only, Het_only, or Both\">"
    };  // initial by common information of header

    ngslib::Fasta fa = ref_file_path;
    std::vector<std::string> contigs;
    for (size_t i(0); i < fa.nseq(); ++i) {
        std::string seqname = fa.iseq_name(i);
        uint32_t seqlen = fa.seq_length(seqname);
        contigs.push_back("##contig=<ID=" + seqname + ",length=" + std::to_string(seqlen) + 
                          ",assembly=" + ref_file_path + ">");
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

        while (f.getline(line)) {
            if (line[0] == '#') continue;
            OUT << line << "\n";
        }
        OUT.flush(); // 确保数据被写入

        if (is_remove_tempfile) ngslib::safe_remove(fn);
    }

    OUT.close();
    return;
}
