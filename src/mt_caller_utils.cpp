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

double ref_vs_alt_ranksumtest(const char ref_base, 
                              const std::string alt_bases_string,
                              const std::vector<char> &bases,
                              const std::vector<int> &values)  // values 和 bases 的值是配对的，一一对应 
{
    std::vector<double> ref, alt;
    ref.reserve(values.size());  // change capacity and save time
    alt.reserve(values.size());  // change capacity and save time

    for (size_t i = 0; i < bases.size(); i++) {
        if (bases[i] == ref_base) {
            ref.push_back(values[i]);
        } else if (alt_bases_string.find(bases[i]) != std::string::npos) {
            alt.push_back(values[i]);
        }
    }
    
    double p_phred_scale_value;
    if (ref.size() > 0 && alt.size() > 0) { // not empty
        double p_value = wilcoxon_ranksum_test(ref, alt);
        p_phred_scale_value = -10 * log10(p_value);
        if (std::isinf(p_phred_scale_value)) {
            p_phred_scale_value = 10000;
        }
    } else {
        // set a big enough phred-scale value if only get REF or ALT base on this postion.
        p_phred_scale_value = 10000;
    }
    return p_phred_scale_value;
}

double ref_vs_alt_ranksumtest(const char ref_base, 
                              const std::string alt_bases_string,
                              const std::vector<char> &bases,
                              const std::vector<char> &values) 
{
    std::vector<int> v(values.begin(), values.end()); 
    return ref_vs_alt_ranksumtest(ref_base, alt_bases_string, bases, v);
}

std::string vcf_header_define(const std::string &ref_file_path, const std::vector<std::string> &samples, const std::string other_comment) {
    std::vector<std::string> header = {
        "##fileformat=VCFv4.2",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth on the REF position\">",
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">",
        "##FORMAT=<ID=HF,Number=A,Type=Float,Description=\"Homoplasmy/Heteroplasmy fraction for the ref and alt alleles in the order listed\">",
        "##FORMAT=<ID=CI,Number=1,Type=String,Description=\"95\% confidence interval around the estimated homoplasmy/heteroplasmy fraction for "
        "the GT alleles in the order listed. format: ci_low,ci_up;ci_low,ci_up;...\">",
        "##FORMAT=<ID=HQ,Number=A,Type=Integer,Description=\"Heteroplasmy Quality, phred quality scores of pvalue of one-tail Fisher exact test "
        "to determine if the rate of heteroplasmy is significantly greater than user defined cutoff (-j), an ordered list of the GT alleles. "
        "[CAUTION] In most cases, the minor allele corresponds to the heteroplasmic allele; therefore, the HQ at the minor allele position "
        "reflects the quality value of heterozygous allele mostly.\">",
        "##FORMAT=<ID=LHF,Number=A,Type=Float,Description=\"Transformed heteroplasmy: The logit of the heteroplasmy fraction (HF) is computed as logit(HF) = ln(HF / (1 - HF)) for each heteroplasmic allele, in the order listed by GT.\">",
        "##FORMAT=<ID=SB,Number=1,Type=String,Description=\"Allele-specific forward/reverse read counts for strand bias tests for the GT alleles in "
        "the order listed, separated by ';'. Format: fwd,rev;fwd,rev;...\">",
        "##FORMAT=<ID=FS,Number=A,Type=Float,Description=\"An ordered, comma delimited list of phred-scaled p-value using Fisher's exact test to detect strand bias\">",
        "##FORMAT=<ID=SOR,Number=A,Type=Float,Description=\"An ordered, comma delimited list of strand bias estimated by the Symmetric Odds Ratio test\">",
        "##FORMAT=<ID=VT,Number=1,Type=String,Description=\"An ordered, comma delimited list of variant type: REF, SNV, INS, DEL, or MNV\">",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"An ordered, comma delimited list of non-reference allele frequencies\">",
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"An ordered, comma delimited list of non-reference allele count in genotypes\">",
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of ref and non-ref allele in called genotypes\">",
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
