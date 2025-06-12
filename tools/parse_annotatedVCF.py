import argparse
import gzip
import re
import logging
import warnings
from multiprocessing import Pool
import pandas as pd
import os

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", message="divide by zero encountered in divide")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

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

def parse_line(line, SAMPLES):
    """Parse a line of VCF file and return a dictionary of values."""
    line_list = []
    CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,*samples = line.strip().split('\t')
    D_loop_region = [i for i in range(1,577)] + [i for i in range(16024,16570)]
    for sample_name, sample_data in zip(SAMPLES, samples):
        trinucleotide  = re.findall(r';?trinucleotide=([^;]+)', INFO)[0]
        try:
            var_type  = re.findall(r'VEP_CSQ=.\|([^|]+)\|', INFO)[0]
        except:
            if len(REF) > 1 or len(ALT) > 1:
                var_type = ''
            else:
                raise ValueError(f"{CHROM}:{POS} {REF} {ALT} var_type not found in INFO field.")
        try:
            gene_symbol = re.findall(r';?mitomap_locus=([^;]+)', INFO)[0]
        except:
            gene_symbol = ''
        # gene_symbol = re.findall(r';?mitomap_locus=([^;]+)', INFO)[0]
        if gene_symbol =='':
            if int(POS) in D_loop_region:
                gene_symbol = 'D-loop'
            else:
                raise ValueError(f"{CHROM}:{POS} Gene symbol not found in INFO field.")
        try:
            # GT,GQ,DP,AD,HF,CI,HQ,LHF,SB,FS,SOR,VT,LODR = sample_data.split(':')
            GT,GQ,DP,AD,HF,*_ = sample_data.split(':')
        except:
            if sample_data.startswith('.'):
                continue
            else:
                raise ValueError(f"{CHROM}:{POS} {sample_data} sample_data is error in VCF file.")
        GT_list = GT.split('/')
        if len(GT_list) == 1:
            continue
        for i in range(len(GT_list)):
            if GT_list[i] == '0' or GT_list[i] == '.' or GT_list[i] == '-65':
                continue
            else:
                gt = GT_list[i]
                try:
                    alt = ALT.split(',')[int(gt)-1]
                except:
                    raise ValueError(f"{ALT, sample_data} ALT field not found in VCF file.")
            hf = float(HF.split(',')[i])
            newREF, newALT = remove_common_suffix(REF, alt)
            line_list.append([sample_name,CHROM,POS,newREF, newALT,hf,trinucleotide,var_type, gene_symbol])
    return line_list


def load_data(input_vcf_file, output_prefix, threads):
    """Read vcf file and convert some information to pandas dataframe."""
    
    with gzip.open(input_vcf_file, "rt") if input_vcf_file.endswith(".gz") else open(input_vcf_file, "r") as IN_VCF, \
        open(output_prefix + ".tmp", "w") as tmpf:
        colname = ["Sample_name","CHROM", "POS", "REF", "ALT", "HF", 'Trinucleotide','Consequence', 'Gene symbol']
        print('\t'.join(colname), file=tmpf)
        
        vcf_list = []
        pool = Pool(processes= threads)
        for line in IN_VCF:
            if line.startswith('##'):
                continue
            
            if line.startswith('#CHROM'):
                SAMPLES = line.strip().split('\t')[9:]
                continue
            vcf_list.append(pool.apply_async(parse_line, args=(line,SAMPLES)))
            
        for line_list in vcf_list:
            for line in line_list.get():
                print(*line, sep='\t', file=tmpf)
        pool.close()
        pool.join()

def merge_info(output_file, information_file):
    """Merge information file with sample information."""
    sample_info = pd.read_csv(information_file) if information_file.endswith('.csv') else pd.read_table(information_file, sep='\t')
    df = pd.read_csv(output_file+".tmp", sep='\t')
    used_col = [x for x in ['Family_Index', 'FamilyID', 'Sample_name', 'Tissue', 'Member'] if x in sample_info.columns]
    sample_info = sample_info[used_col]
    df = pd.merge(sample_info, df,  on='Sample_name', how='right')
    df.sort_values(['Family_Index'], inplace=True)
    df.to_csv(output_file, sep='\t', index=False)
    os.remove(output_file+".tmp")
    
def main():
    parser = argparse.ArgumentParser(description='Parse annotated VCF file and extract HF information.')
    parser.add_argument('-i', '--input_vcf_file', required=True, help='Input VCF file.')
    parser.add_argument('-o', '--output_file', required=True, help='Output file.')
    parser.add_argument('-t', '--threads', default=1, type=int, help='Number of threads to use.')
    parser.add_argument('-f', '--information_file', type=str, default=None, help='the information file Path.')
    args = parser.parse_args()
    
    if args.information_file:
        load_data(args.input_vcf_file, args.output_file, args.threads)
        merge_info(args.output_file, args.information_file)
    else:
        load_data(args.input_vcf_file, args.output_file, args.threads)
        os.rename(args.output_file+".tmp", args.output_file)

if __name__ == '__main__':
    main()

