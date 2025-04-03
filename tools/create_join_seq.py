"""
Author: Shujia Huang
Date: 2024-12-25
"""
import pysam
import argparse
import sys
import os


def read_bed_file(bed_file):
    """Read BED file and return a sorted list of intervals"""
    intervals = []
    with open(bed_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):

                # Parse BED format (chromosome, start, end)
                parts = line.strip().split()
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                intervals.append((chrom, start, end))

    return sorted(intervals, key=lambda x: (x[0], x[1]))


def extract_and_join_sequences(fasta_file, bed_file, output_file, n_length=1000, line_width=70, 
                               output_name="combined_sequences", verbose=False):
    """
    Extract sequences from reference genome based on BED intervals and combine them with N spacers
    
    Args:
        fasta_file: Input reference genome file
        bed_file: Input BED file
        output_file: Output FASTA file path
        n_length: Length of N spacer between intervals
        line_width: Number of bases per line in output FASTA
        output_name: Name for the output sequence
        verbose: Whether to print verbose output
    """
    # Open reference genome file
    try:
        fasta = pysam.FastaFile(fasta_file)
    except Exception as e:
        print(f"Error opening reference file: {e}")
        return

    # Read intervals from BED file
    intervals = read_bed_file(bed_file)
    if not intervals:
        print("No valid intervals found in BED file")
        return
    
    if verbose:
        print(f"Found {len(intervals)} intervals in BED file")
    
    # Extract sequences and join them with N gaps
    joined_seq = ''
    for i, (chrom, start, end) in enumerate(intervals):
        # Extract sequence for current interval
        try:
            seq = fasta.fetch(chrom, start, end)
        except Exception as e:
            print(f"Warning: Could not extract sequence for interval: {chrom}:{start}-{end}")
            print(f"Error: {str(e)}")
            continue
            
        # Add N gap if not first sequence
        if i > 0:
            joined_seq += 'N' * n_length
        
        joined_seq += seq
    
    # Close reference file
    fasta.close()

    seq_len = len(joined_seq)
    # Write output in FASTA format
    with open(output_file, 'w') as out:
        # Write header
        out.write(f">{output_name}\n")
        
        # Write sequence in lines of {line_width} characters
        for i in range(0, seq_len, line_width):
            out.write(joined_seq[i:i+line_width] + '\n')

    if verbose:
        print(f"Total sequence length: {seq_len} bp")
        print(f"Combined sequence written to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract sequences from a reference genome based on BED intervals "
                    "and combine them with N spacers",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('-r', '--reference', required=True,
                        help='Input reference genome in FASTA format')
    parser.add_argument('-b', '--bed', required=True,
                        help='Input BED file with intervals')
    parser.add_argument('-o', '--output', required=True,
                        help='Output FASTA file')
    
    # Optional arguments
    parser.add_argument('-n', '--n-length', type=int, default=1000,
                        help='Length of N spacer between intervals')
    parser.add_argument('-w', '--line-width', type=int, default=70,
                        help='Number of bases per line in output FASTA')
    parser.add_argument('-s', '--sequence-name', default="combined_sequences",
                        help='Name for the output sequence')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print verbose output')
    
    args = parser.parse_args()
    
    # Check if input files exist
    if not os.path.exists(args.reference):
        print(f"Error: Reference file '{args.reference}' does not exist")
        sys.exit(1)
    if not os.path.exists(args.bed):
        print(f"Error: BED file '{args.bed}' does not exist")
        sys.exit(1)
    
    # Check and create output directory if needed
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        print(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir)
    
    extract_and_join_sequences(
        args.reference,
        args.bed,
        args.output,
        n_length=args.n_length,
        line_width=args.line_width,
        output_name=args.sequence_name,
        verbose=args.verbose
    )


if __name__ == "__main__":
    main()
