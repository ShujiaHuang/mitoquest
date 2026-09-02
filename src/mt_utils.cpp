#include "mt_utils.h"

#include "io/fasta.h"
#include "io/iobgzf.h"
#include "algorithm.h"

std::string format_double(double value, int precision) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(precision) << value;
    return ss.str();
}

namespace {

std::string format_logit_af(double allele_frequency) {
    if (!(allele_frequency > 0.0 && allele_frequency < 1.0)
        || !std::isfinite(allele_frequency)) 
    {
        return ".";
    }
    return format_double(std::log(allele_frequency / (1.0 - allele_frequency)), 6);
}

double binomial_survival_probability(int successes, int trials, double probability) {
    if (successes <= 0) return 1.0;
    if (successes > trials || trials <= 0 || probability <= 0.0) return 0.0;
    if (probability >= 1.0) return 1.0;

    // Pr(X >= successes) for X ~ Binomial(trials, probability) is the
    // regularized incomplete beta I_probability(successes, trials-successes+1).
    const double tail = kf_betai(static_cast<double>(successes),
                                 static_cast<double>(trials - successes + 1),
                                 probability);
    return std::max(0.0, std::min(1.0, tail));
}

std::string format_strand_metric(double value, int precision) {
    if (!std::isfinite(value)) return ".";
    if (value == 10000.0) return "10000";
    return format_double(value, precision);
}

// Phred-scale a p-value, capping the underflow at the 10,000 sentinel.
double phred_scale_pvalue(double p_value) {
    const double fs = -10.0 * std::log10(p_value);
    if (!std::isfinite(fs)) return 10000.0;
    if (fs == 0.0) return 0.0;
    return fs;
}

}  // namespace

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

    const bool monomorphic = (alt_base == major_base);

    double fs(0.0);
    double sor(0.0);
    if (monomorphic) {
        // Monomorphic genotype: no second allele to contrast against.  Test
        // the absolute strand balance of the sole allele against the 50:50
        // expectation with a two-sided binomial test, and report the
        // symmetric log ratio ln(x + 1/x) of its pseudocount-corrected
        // forward/reverse counts.
        const int n_trials = maj_fwd + maj_rev;
        const double p_two_sided = std::min(
            1.0, 2.0 * binomial_survival_probability(std::max(maj_fwd, maj_rev), n_trials, 0.5)
        );
        fs = phred_scale_pvalue(p_two_sided);

        const double ratio = static_cast<double>(maj_fwd + 1) / static_cast<double>(maj_rev + 1);
        sor = std::log(ratio + 1.0 / ratio);
    } else {
        const double fisher_pvalue = fisher_exact_test(maj_fwd, maj_rev, alt_fwd, alt_rev);
        fs = phred_scale_pvalue(fisher_pvalue);

        // GATK's symmetric odds ratio: add a pseudocount to avoid singular
        // tables, combine the odds ratio with its reciprocal, and normalize
        // for strand imbalance shared by the major and selected alleles.
        const double major_forward = static_cast<double>(maj_fwd) + 1.0;
        const double major_reverse = static_cast<double>(maj_rev) + 1.0;
        const double alt_forward = static_cast<double>(alt_fwd) + 1.0;
        const double alt_reverse = static_cast<double>(alt_rev) + 1.0;

        const double odds_ratio = (major_forward / major_reverse) * (alt_reverse / alt_forward);
        const double symmetric_ratio = odds_ratio + 1.0 / odds_ratio;

        const double major_ratio = std::min(major_forward, major_reverse) / std::max(major_forward, major_reverse);
        const double alt_ratio = std::min(alt_forward, alt_reverse) / std::max(alt_forward, alt_reverse);
        sor = std::log(symmetric_ratio) + std::log(major_ratio) - std::log(alt_ratio);
    }

    StrandBiasInfo sbi;
    sbi.fwd = monomorphic ? maj_fwd : alt_fwd;
    sbi.rev = monomorphic ? maj_rev : alt_rev;
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
    
    for (size_t i = 0; i < ref_alt_order.size(); i++) {  // i == 0 represents the REF GT
        const auto& alt = ref_alt_order[i];
        for (size_t j = 0; j < var_info.alt_bases.size(); j++) {
            if (var_info.alt_bases[j] == alt) {
                sa.gtcode.push_back(i);
                sa.sample_alts.push_back(alt);
                sa.allele_depths.push_back(var_info.depths[j]);
                
                // double h = var_info.freqs[j]; // 这里不要用 lrt 计算出来的 allele frequency，因为可能不知为何会有负数（极少情况下）
                // calculate the allele frequency by allele_depth/total_depth
                const double h = double(var_info.depths[j]) / double(var_info.total_depth);
                sa.allele_freqs.push_back(h);
                sa.logit_af.push_back(format_logit_af(h));
                
                sa.ci_strings.push_back(format_double(var_info.ci[j].first, 4) + "," + 
                                        format_double(var_info.ci[j].second, 4));
                
                sa.sb_strings.push_back(std::to_string(var_info.strand_bias[j].fwd) + "," + 
                                        std::to_string(var_info.strand_bias[j].rev));
                
                sa.fs_strings.push_back(format_strand_metric(var_info.strand_bias[j].fs, 3));
                sa.sor_strings.push_back(format_strand_metric(var_info.strand_bias[j].sor, 3));
                
                sa.var_types.push_back(var_info.var_types[j]);

                const double aq_pvalue = binomial_survival_probability(
                    var_info.depths[j], var_info.total_depth, hf_cutoff);
                const double aq = (aq_pvalue > 0.0) ? -10.0 * std::log10(aq_pvalue) : 10000.0;
                sa.aq.push_back(std::isfinite(aq) ? static_cast<int>(aq) : 10000);
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
        // ---- Per-sample variable-length FORMAT fields --------------------------
        // The following fields all share `mitoquest caller`'s GT-aligned
        // per-sample layout: each sample emits one value per allele present
        // in that sample's GT, in GT order.  The number of values therefore
        // varies per sample (1 for homoplasmic calls, 2 for heteroplasmic,
        // 3 for tri-allelic, etc.) and does NOT match the standard
        // `Number=A` (= n_alts, fixed across samples) cardinality.  Per the
        // VCFv4.2 spec the correct cardinality for per-sample variable-
        // length fields is `Number=.`.
        "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depth for each allele, in the order listed by GT\">",
        "##FORMAT=<ID=AF,Number=.,Type=Float,Description=\"Allele fraction for ref- and alt-alleles, in the order listed by GT\">",
        "##FORMAT=<ID=CI,Number=1,Type=String,Description=\"95\% confidence interval around the estimated allele fraction for "
            "the allele in the order listed by GT. format: ci_low,ci_up;ci_low,ci_up;...\">",
        "##FORMAT=<ID=AQ,Number=.,Type=Integer,Description=\"Allele quality: Phred-scaled exact one-sided Binomial p-value for observing at least AD reads under allele fraction equal to the user cutoff (-j), in the order listed by GT\">",
        "##FORMAT=<ID=LAF,Number=.,Type=Float,Description=\"Transformed AF: logit(AF) = ln(AF/(1-AF)) for each allele in GT order; missing at AF=0 or AF=1\">",
        "##FORMAT=<ID=SB,Number=1,Type=String,Description=\"Allele-specific forward/reverse read counts for strand bias tests for the alleles, in "
            "the order listed by GT, separated by ';'. Format: fwd,rev;fwd,rev;...\">",
        "##FORMAT=<ID=FS,Number=.,Type=Float,Description=\"An ordered, comma delimited list of phred-scaled p-values of the strand-bias test, in the order listed by GT: Fisher's exact test against the major allele for heteroplasmic calls, or a two-sided Binomial test against the 50:50 expectation for the sole allele of a monomorphic call\">",
        "##FORMAT=<ID=SOR,Number=.,Type=Float,Description=\"An ordered, comma delimited list of GATK-style pseudocount-corrected symmetric strand odds ratios, in the order listed by GT; for the sole allele of a monomorphic call, the symmetric log ratio ln(x+1/x) of its pseudocount-corrected forward/reverse counts\">",
        "##FORMAT=<ID=VT,Number=1,Type=String,Description=\"An ordered, comma delimited list of variant type: REF, SNV, INS, DEL, or MNV\">",
        "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with non-missing GT at this site\">",
        "##INFO=<ID=REF_N,Number=1,Type=Integer,Description=\"Total number of individuals exhibiting the reference state in the population\">",
        "##INFO=<ID=HET_N,Number=1,Type=Integer,Description=\"Total number of individuals exhibiting the heteroplasmic state in the population\">",
        "##INFO=<ID=HOM_N,Number=1,Type=Integer,Description=\"Total number of individuals exhibiting the homoplasmic state in the population\">",
        "##INFO=<ID=DP_MEAN,Number=1,Type=Integer,Description=\"Mean mitochondrial sequencing depth across samples contributing to NS\">",
        "##INFO=<ID=DP_MEDIAN,Number=1,Type=Integer,Description=\"Median mitochondrial sequencing depth across samples contributing to NS\">",
        "##INFO=<ID=VAF_MEAN,Number=A,Type=Float,Description=\"Mean mitochondrial variant allele (non-reference alleles) fraction(VAF) across all samples "
            "contributing to NS, with VAF=0 assigned to samples without detectable variant, denominator is the count of individuals with non-missing genotype\">",
        "##INFO=<ID=VAF_MEAN_HET,Number=A,Type=Float,Description=\"Mean mitochondrial VAF among heteroplasmic samples only, denominator is the count of heteroplasmic individuals\">",
        "##INFO=<ID=PT,Number=1,Type=String,Description=\"Type of plasmicity observed in population: Ref, Hom, Het, or Mixed(Hom and Het)\">"
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
