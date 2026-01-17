#!/usr/bin/env python3

import sys
import gzip

def rewrite_vcf(
    in_vcf,
    out_vcf
):
    NEW_HEADER = [
        "##fileformat=VCFv4.2",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth on the REF position\">",
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depth for each allele, in the order listed by GT\">",
        "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele fraction for each allele, in the order listed by GT. Fraction of non-reference allele corresponds to the variant allele fraction(VAF)\">",
        "##FORMAT=<ID=CI,Number=1,Type=String,Description=\"95% confidence interval around the estimated allele fraction for the allele in the order listed by GT. format: ci_low,ci_up;ci_low,ci_up;...\">",
        "##FORMAT=<ID=AQ,Number=A,Type=Integer,Description=\"Allele quality, phred quality scores of pvalue of one-tail Fisher exact test to determine if the rate of allele is significantly greater than user defined cutoff (-j), in the order listed by GT. [CAUTION] In most cases, the minor allele corresponds to the heteroplasmic allele; therefore, the AQ at the minor allele position reflects the quality value of heterozygous allele mostly.\">",
        "##FORMAT=<ID=LAF,Number=A,Type=Float,Description=\"Transformed AF: The logit of the allele fraction (AF) is computed as logit(AF) = ln(AF/(1-AF)) for each allele, in the order listed by GT.\">",
        "##FORMAT=<ID=SB,Number=1,Type=String,Description=\"Allele-specific forward/reverse read counts for strand bias tests for the alleles, in the order listed by GT, separated by ';'. Format: fwd,rev;fwd,rev;...\">",
        "##FORMAT=<ID=FS,Number=A,Type=Float,Description=\"An ordered, comma delimited list of phred-scaled p-value using Fisher's exact test to detect strand bias\">",
        "##FORMAT=<ID=SOR,Number=A,Type=Float,Description=\"An ordered, comma delimited list of strand bias estimated by the Symmetric Odds Ratio test\">",
        "##FORMAT=<ID=VT,Number=1,Type=String,Description=\"An ordered, comma delimited list of variant type: REF, SNV, INS, DEL, or MNV\">",
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Number of samples with non-missing GT at this site\">",
        "##INFO=<ID=REF_N,Number=1,Type=Integer,Description=\"Total number of individuals exhibiting the reference state in the population.\">",
        "##INFO=<ID=HOM_N,Number=1,Type=Integer,Description=\"Total number of individuals exhibiting the homoplasmic state in the population.\">",
        "##INFO=<ID=HET_N,Number=1,Type=Integer,Description=\"Total number of individuals exhibiting the heteroplasmic state in the population.\">",
        "##INFO=<ID=DP_MEAN,Number=1,Type=Float,Description=\"Mean mitochondrial sequencing depth across samples contributing to AN\">",
        "##INFO=<ID=DP_MEDIAN,Number=1,Type=Integer,Description=\"Median mitochondrial sequencing depth across samples contributing to AN\">",
        "##INFO=<ID=VAF_MEAN,Number=A,Type=Float,Description=\"Mean mitochondrial variant allele fraction(VAF) across all samples contributing to AN, with VAF=0 assigned to samples without detectable variant\">",
        "##INFO=<ID=VAF_MEDIAN,Number=A,Type=Float,Description=\"Median mitochondrial VAF across all samples contributing to AN\">",
        "##INFO=<ID=VAF_MEAN_HET,Number=A,Type=Float,Description=\"Mean mitochondrial VAF among heteroplasmic samples only\">",
        "##INFO=<ID=VAF_MEDIAN_HET,Number=A,Type=Float,Description=\"Median mitochondrial VAF among heteroplasmic samples only\">",
        "##INFO=<ID=PT,Number=1,Type=String,Description=\"Type of plasmicity observed in population: Ref, Hom, Het, or Mixed(Hom and Het)\">",
        "##contig=<ID=chrM,length=16569,assembly=chrM_rCRS.decoy.fa.gz>",
        "##contig=<ID=NUMT_JoinedSequences_gaps1000N_decoy,length=64266,assembly=chrM_rCRS.decoy.fa.gz>",
        "##reference=file:///Users/huangshujia/Projects/mitoquest/tests/data/chrM_rCRS.decoy.fa.gz>",
    ]

    INFO_TEMPLATE = (
        "AN=0;"
        "REF_N=0;"
        "HOM_N=0;"
        "HET_N=0;"
        "DP_MEAN=0;"
        "DP_MEDIAN=0;"
        "VAF_MEAN=0;"
        "VAF_MEDIAN=0;"
        "VAF_MEAN_HET=0;"
        "VAF_MEDIAN_HET=0;"
        "PT=Ref"
    )

    with gzip.open(in_vcf, 'rt') if in_vcf.endswith(".gz") else open(in_vcf) as fin, open(out_vcf, "w") as fout:
        for line in fin:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                for h in NEW_HEADER:
                    fout.write(h + "\n")
                fout.write(line)
                continue

            fields = line.rstrip().split("\t")

            # INFO column
            fields[7] = INFO_TEMPLATE

            # FORMAT column
            if fields[8] == "GT:GQ:DP:AD:HF:CI:HQ:LHF:SB:FS:SOR:VT":
                fields[8] = "GT:GQ:DP:AD:AF:CI:AQ:LAF:SB:FS:SOR:VT"

            fout.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} input.vcf output.vcf")

    rewrite_vcf(sys.argv[1], sys.argv[2])

