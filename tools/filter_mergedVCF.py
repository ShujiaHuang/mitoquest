from collections import Counter
import argparse
import subprocess

def rm_region(POS, region_list=[[299,317], [511,524],[564, 571],[952,955],[3106, 3108],[8268, 8279],[13646,13650], [16179,16183]]):
    """
    Site Blocklist(artifact-prone-sites):  299–317, 511–524, 564-571, 952-955(specific), 3106–3108(3107 is a N base), 
    8268-8279, 13646-13650(specific) and 16179–16183 (May cause by Ploy-C tracts). Total length:  70bp
    """
    # chrom, pos, id, ref, alt, qual, filter, info, format, sample = line.strip().split('\t')
    poses = []
    for list in region_list:
        for i in range(list[0], list[1]+1):
            poses.append(i)
    if int(POS) in poses:
        return False
    else:
        return True


def modify_single_sample(sample, HF_cutoff, HQ_cutoff, AD_cutoff):
    """
    Modify the Format value of sample after filtering.
    """
    try:
        GT,GQ,DP,AD,HF,CI,HQ,LHF,SB,FS,SOR,VT = sample.split(':')
    except:
        GT = ''
        return sample, GT
    GT_list = GT.split('/')
    AD_list = AD.split(',')
    HF_list = HF.split(',')
    CI_list = CI.split(';')
    HQ_list = HQ.split(',')
    LHF_list = LHF.split(',')
    SB_list = SB.split(';')
    FS_list = FS.split(',')
    SOR_list = SOR.split(',')
    VT_list = VT.split(',')
    index_list = []
    for index in range(len(GT.split('/'))):
        if AD_list[index] == '.' or HF_list[index] == '.' or HQ_list[index] == '.':
            continue
        if float(AD_list[index]) >= AD_cutoff and float(HF_list[index]) >= HF_cutoff and float(HQ_list[index]) >= HQ_cutoff:
            index_list.append(index)
    GT = '/'.join([str(GT_list[i]) for i in index_list])
    AD = ','.join([str(AD_list[i]) for i in index_list])
    HF = ','.join([str(HF_list[i]) for i in index_list])
    CI = ';'.join([CI_list[i] for i in index_list])
    HQ = ','.join([str(HQ_list[i]) for i in index_list])
    LHF = ','.join([str(LHF_list[i]) for i in index_list])
    SB = ';'.join([SB_list[i] for i in index_list])
    FS = ','.join([FS_list[i] for i in index_list])
    SOR = ','.join([SOR_list[i] for i in index_list])
    VT = ','.join([VT_list[i] for i in index_list])
    # if GT and AD and HF and CI and HQ and LHF and SB and FS and SOR and VT:
    #     GT = AD = HF = CI = HQ = LHF = SB = FS = SOR = VT = '.'
    sample = ':'.join([GT, GQ, DP, AD, HF, CI, HQ, LHF, SB, FS, SOR, VT])
    return sample, GT

def modify_ALT_INFO(ALT, GT_list):
    """
    Modify the ALT and INFO value of sample after filtering.
    """
    all_alleles_list = [int(x) for GT in GT_list for x in GT.split('/')]
    AN = len(all_alleles_list)
    ALT_list = []
    AF_list = []
    AC_list = []
    for allele in sorted(Counter(all_alleles_list)):
        if allele != 0:
            ALT_list.append(ALT.split(',')[allele-1])
            AC_list.append(Counter(all_alleles_list)[allele])
            AF_list.append(Counter(all_alleles_list)[allele]/AN)
    ALT_str = ','.join(ALT_list)
    AF = ','.join([f"{num:.4f}" for num in AF_list])
    AC = ','.join([f"{num}" for num in AC_list])
    INFO =  'AF=' + AF + ';AC=' + AC + ';AN=' + str(AN)
    return ALT_str, INFO

def update_gt(gt, index_map):
    """
    modify the GT in per-sample format accordinly to the newALT and oldALT after filtering.
    """
    parts = gt.split('/')
    updated_parts = []
    for part in parts:
        if part == '.' or part == '0' or part == '':
            updated_parts.append(part)
        else:
            # parse old index and map to new index
            # try:
            old_alt_idx = int(part)
            # if old_alt_idx in index_map:
            updated_parts.append(str(index_map[old_alt_idx]))
            # else:
            #     raise ValueError  # invalid index
            # except ValueError:
            #     updated_parts.append(part)  # invalid index（如 '.'）
    newgt = '/'.join(updated_parts)
    if newgt:
        return newgt
    else:
        return '.'

def modify_line(line, HF_cutoff, HQ_cutoff, AD_cutoff):
    """
    Modify the INFO value of sample after filtering.
    """
    GT_list = []
    samples_list = []
    CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,*samples = line.strip().split('\t')
    for sample in samples:
        sample,GT = modify_single_sample(sample, HF_cutoff, HQ_cutoff, AD_cutoff)
        if GT:
            GT_list.append(GT)
        samples_list.append(sample)
    newALT, newINFO = modify_ALT_INFO(ALT, GT_list)
    # try to modify the GT in per-sample format accordinly to the newALT and oldALT after filtering
    old_alt_items = ALT.split(',')
    new_alt_items = newALT.split(',')
    index_map = {}
    new_idx = 0
    for old_idx, item in enumerate(old_alt_items):
        if item in new_alt_items:
            index_map[old_idx + 1] = new_idx + 1
            new_idx += 1
    new_samples_list = []
    for sample in samples_list:
        GT,*_ = sample.split(':')
        new_GT = update_gt(GT, index_map)
        new_samples_list.append(':'.join([new_GT,*_]))
    line = '\t'.join([CHROM,POS,ID,REF,newALT,QUAL,FILTER,INFO,FORMAT]+new_samples_list)
    return line, newALT

def main():
    import argparse
    import gzip
    parser = argparse.ArgumentParser(description='Filter variants based on the following 3 criteria: 1. Non blocklisted regions; 2. HF, AD and LHF cutoffs; 3. ALT number need lt 3 after filtering by AD, HF, and LHF.')
    parser.add_argument('-i', '--input_vcf', type=str, required=True, help='Input vcf file')
    parser.add_argument('-f', '--HF_cutoff',  type=float, required=True, help='HF cutoff')
    parser.add_argument('-q', '--HQ_cutoff', type=float, required=True, help='LHF cutoff')
    parser.add_argument('-d', '--AD_cutoff',  type=float, required=True, help='AD cutoff')
    parser.add_argument('-o', '--output_vcf', type=str, required=True, help='Output vcf file')
    args = parser.parse_args()
    with gzip.open(args.input_vcf, 'rt') if args.input_vcf.endswith('.gz') else open(args.input_vcf, 'r') as inputf, open(args.output_vcf, 'w') as outf:
        total_line = 0
        block_poses = []
        filtered_line = []
        filtered_ALT = []
        for line in inputf:
            if line.startswith('##'):
                outf.write(line)
            elif line.startswith('#CHROM'):
                print(f'##filtering_command=python filtering_mergedVCF.py -i {args.input_vcf} -f {args.HF_cutoff} -q {args.HQ_cutoff} -d {args.AD_cutoff} -o {args.output_vcf}', file=outf)
                outf.write(line)
            else:
                total_line += 1
                if rm_region(line.strip().split('\t')[1]):
                    line, newALT = modify_line(line, args.HF_cutoff, args.HQ_cutoff, args.AD_cutoff)
                    if newALT:
                        if len(newALT.strip().split(',')) <3:
                            print(line.strip(), file=outf)
                        else:
                            filtered_ALT.append(line.strip().split('\t')[1])
                    else:
                        filtered_line.append(line.strip().split('\t')[1])
                else:
                    block_poses.append(line.strip().split('\t')[1])
        print(f"Total sites: {total_line}")
        print(f"Blocked sites: {len(block_poses)}:\n{block_poses}")
        print(f"Filtered sites by AD-{args.AD_cutoff}, HF-{args.HF_cutoff}, and HQ-{args.HQ_cutoff} cutoffs: {len(filtered_line)}:\n{filtered_line}")
        print(f"Filtered sites by ALT (ALT need lt 3): {len(filtered_ALT)}:\n{filtered_ALT}")
        print(f'Total remained sites: {total_line - len(block_poses) - len(filtered_line) - len(filtered_ALT)}')
    subprocess.run(["bcftools", "view", "-Oz", "-o", args.output_vcf+".gz", args.output_vcf])
    subprocess.run(["bcftools", "index", args.output_vcf+".gz"]) 
if __name__ == '__main__':
    main()
    
