import gzip
from collections import defaultdict
import argparse

def open_file(file_path):
    """
    Open a file, regardless of whether it is gzipped or not.
    """
    if file_path.endswith('.gz'):
        file = gzip.open(file_path, 'rt')
    else:
        file = open(file_path, 'r')
    return file

def vcf_header(vcf_file):
    """
    Get the header of a VCF file.
    """
    vcf = open_file(vcf_file)
    header = [lines for lines in vcf if lines.startswith('#')]
    vcf.close()
    return header

def get_region_vcf(vcf_file, region, CR=False):
    """
    Get the variation site of the corresponding region for NCR or CR.
    """
    vcf = open_file(vcf_file)
    with open(region, 'r') as f:
        for line in f:
            start = line.strip().split('\t')[1]
            end = line.strip().split('\t')[2]
            region_list = [pos for pos in range(int(start), int(end)+1)]
    vcf_dict = defaultdict()
    for line in vcf:
        if line.startswith('#'):
            continue
        else:
            CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,*_sample = line.strip().split('\t')
            if int(POS) in region_list:
                if CR:
                    if int(POS) <= 8569:
                        POS = str(int(POS) + 8000)
                        _line = '\t'.join([CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,*_sample])
                        vcf_dict[int(POS)] = _line.strip()
                    else:
                        POS = str(int(POS) - 8569)
                        _line = '\t'.join([CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,*_sample])
                        vcf_dict[int(POS)] = _line.strip()
                else:
                    vcf_dict[int(POS)] = line.strip()
    vcf.close()
    return vcf_dict

def main():
    parser = argparse.ArgumentParser(description='merge raw and shifting vcf files')
    parser.add_argument('-v1', '--raw_vcf', help='the raw vcf file', required=True)
    parser.add_argument('-v2', '--shift_vcf', help='the shifted vcf file', required=True)
    parser.add_argument('-r1', '--raw_region', help='the raw region file', required=True)
    parser.add_argument('-r2', '--shift_region', help='the shifted region file', required=True)
    parser.add_argument('-o', '--output', help='the output file', required=True)
    argparse_args = parser.parse_args()

    raw_vcf = argparse_args.raw_vcf
    shift_vcf = argparse_args.shift_vcf
    raw_region = argparse_args.raw_region
    shift_region = argparse_args.shift_region
    output = argparse_args.output

    raw_header = vcf_header(raw_vcf)
    ncr_vcf_dict = get_region_vcf(raw_vcf, raw_region)
    cr_vcf_dict = get_region_vcf(shift_vcf, shift_region, CR=True)
    ncr_vcf_dict.update(cr_vcf_dict) # merge two vcf files
    sorted_dict = dict(sorted(ncr_vcf_dict.items())) # sort the vcf file by position
    
    with open(output, 'w') as f:
        f.write(''.join(raw_header))
        for pos in sorted_dict:
            f.write(sorted_dict[pos] + '\n')

if __name__ == '__main__':
    main()