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

StrandBiasInfo strand_bias(const std::string &ref_base, 
                           const std::string &alt_bases_string,
                           const std::vector<std::string> &bases,
                           const std::vector<char> &strands)  // strands 和 bases 是配对的，一一对应 
{
    int ref_fwd = 0, ref_rev = 0;
    int alt_fwd = 0, alt_rev = 0;

    for (size_t i(0); i < bases.size(); ++i) {
        if (strands[i] == '+') {
            if (bases[i] == ref_base) {
                ++ref_fwd;
            } else if (alt_bases_string == bases[i]) {
                ++alt_fwd;
            }

        } else if (strands[i] == '-') {
            if (bases[i] == ref_base) {
                ++ref_rev;
            } else if (alt_bases_string == bases[i]) {
                ++alt_rev;
            }

        } else {
            throw std::runtime_error("[ERROR] Get strange strand symbol: " + std::to_string(strands[i]));
        }
    }

    // 如果是'全 Ref' 或者是 '全 ALT' 怎么办？
    double fs = -10 * log10(fisher_exact_test(ref_fwd, ref_rev, alt_fwd, alt_rev));
    if (std::isinf(fs)) {
        fs = 10000;
    } else if (fs == 0) {
        fs = 0.0;
    }

    // Strand bias estimated by the Symmetric Odds Ratio test
    // https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_annotator_StrandOddsRatio.php
    double sor = (ref_rev * alt_fwd > 0) ? (double)(ref_fwd * alt_rev) / (double)(ref_rev * alt_fwd): 10000;

    StrandBiasInfo sbi;
    sbi.ref_fwd = ref_fwd; sbi.ref_rev = ref_rev; 
    sbi.alt_fwd = alt_fwd; sbi.alt_rev = alt_rev;
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

std::string vcf_header_define(const std::string &ref_file_path, const std::vector<std::string> &samples) {
    std::vector<std::string> header = {
        "##fileformat=VCFv4.2",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth on the REF position\">",
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">",
        "##FORMAT=<ID=HF,Number=A,Type=Float,Description=\"Heteroplasmy fraction for the ref and alt alleles in the order listed\">",
        "##FORMAT=<ID=CI,Number=1,Type=String,Description=\"An ordered, '|' delimited list of the 95\% confidence interval around the estimated heteroplasmy fraction for the ref and alt alleles in the order listed. format: ci_low,ci_up|ci_low,ci_up|...\">",
        "##FORMAT=<ID=SB,Number=1,Type=String,Description=\"Read number of mapping strand orientation for ref and alt alleles in the order listed: ref_fwd,ref_rev,alt_fwd,alt_rev|...\">",
        "##FORMAT=<ID=FS,Number=A,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test to detect strand bias\">",
        "##FORMAT=<ID=SOR,Number=A,Type=Float,Description=\"Strand bias estimated by the Symmetric Odds Ratio test\">",
        "##FORMAT=<ID=VT,Number=1,Type=String,Description=\"Variant type: REF, SNV, INS, DEL, or MNV\">",

        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"An ordered, comma delimited list of allele frequencies base\">",
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"An ordered, comma delimited allele count in genotypes\">",
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">",
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
    header.push_back("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + ngslib::join(samples, "\t"));

    return ngslib::join(header, "\n");
}

void merge_file_by_line(const std::vector<std::string> &infiles, const std::string &outfile, 
                        std::string header, bool is_remove_tempfile) 
{
    if (infiles.empty()) return;

    bool is_compress_out = (ngslib::suffix_name(outfile) == ".gz") ? true : false;
    BGZF *OUT = bgzf_open(outfile.c_str(), is_compress_out ? "w" : "uw");  // output file
    if (!OUT) throw std::runtime_error("[ERROR] " + outfile + " open failure.");

    header += "\n";
    if (bgzf_write(OUT, header.c_str(), header.length()) != header.length())
        throw std::runtime_error("[ERROR] fail to write data");

    /* Merge all files here */
    for (auto fn: infiles) {
        BGZF *f = bgzf_open(fn.c_str(), "r");
        kstring_t s; s.s = NULL; s.l = s.m = 0;
        while (bgzf_getline(f, '\n', &s) >= 0) {
            if (s.s[0] == '#') continue;  // ignore the header of subfiles.
            std::string out(s.s); out += "\n";

            if (bgzf_write(OUT, out.c_str(), out.length()) != out.length())
                throw std::runtime_error("[ERROR] fail to write data");
        }
        bgzf_close(f);

        if (is_remove_tempfile) ngslib::safe_remove(fn);
    }
    
    int is_cl = bgzf_close(OUT);
    if (is_cl < 0) throw std::runtime_error("[ERROR] " + outfile + " fail close.");
    
    return;
}
