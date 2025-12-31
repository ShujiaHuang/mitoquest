import sys
import gzip
"""
This script is used to parse vcf file and extract information of each sample.
"""
if len(sys.argv) != 3:
    print("Usage: python parse_vcf.py vcf_file parse_file")
    sys.exit(1)

vcf_file = sys.argv[1]
parse_file = sys.argv[2]

def open_file(file_path):
    """
    Open a file, regardless of whether it is gzipped or not.
    """
    if file_path.endswith('.gz'):
        file = gzip.open(file_path, 'rt')
    else:
        file = open(file_path, 'r')
    return file

def remove_common_suffix(s1, s2):
    """
    Remove the common suffix of two strings.
    """
    i = 0
    while i < min(len(s1), len(s2)) and s1[-1 - i] == s2[-1 - i]:
        i += 1

    if i > 0:
        if s1[:-i] and s2[:-i]:
            new_s1 = s1[:-i]
            new_s2 = s2[:-i]
        else:
            new_s1 = s1[:-i+1]
            new_s2 = s2[:-i+1]
    else:
        new_s1 = s1
        new_s2 = s2
    return (new_s1, new_s2)

inputf = open_file(vcf_file)
with open(parse_file, 'w') as outputf:
    print("Sample_name","allGT", "snv_type", "CHROM", "POS", "REF", "ALT", "GT", "AD", "HF", "HQ", "LHF", "SB", "VT", sep='\t', file=outputf)
    for line in inputf:
        if line.startswith('##'):
            continue
        elif line.startswith('#CHROM'):
            _CHROM,_POS,_ID,_REF,_ALT,_QUAL,_FILTER,_INFO,_FORMAT,*SAMPLES = line.strip().split('\t')
        else:
            CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,*samples = line.strip().split('\t')
            for sample_name, sample_data in zip(SAMPLES, samples):
                if sample_data.startswith('.'):
                    continue
                GT,GQ,DP,AD,HF,CI,HQ,LHF,SB,FS,SOR,VT = sample_data.split(':')
                GT_list = GT.split('/')
                if len(GT_list) == 1 and GT != '0':
                    snv_type = 'SNP'
                else:
                    snv_type = 'HF'
                for i in range(len(GT_list)):
                    if GT_list[i] == '0':
                        continue
                    gt = GT_list[i]
                    alt = ALT.split(',')[int(gt)-1]
                    ad = int(AD.split(',')[i])
                    hf = float(HF.split(',')[i])
                    hq = HQ.split(',')[i]
                    lhf = float(LHF.split(',')[i])
                    sb = SB.split(';')[i]
                    vt = VT.split(',')[i]
                    
                    REF_simplified, alt_simplified = remove_common_suffix(REF, alt)
                    print(sample_name,GT, snv_type,CHROM,POS,REF_simplified,alt_simplified,gt,ad,hf,hq,lhf,sb,vt, sep='\t', file=outputf)

inputf.close()