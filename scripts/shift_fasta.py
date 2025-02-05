"""Perform shift operation on FASTA format sequence file

Author: Shujia Huang
Date: 2024-12-25
"""
import argparse


def shift_fasta_sequence(input_file, output_file, shift_size, line_width):
    # Read input file
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Get header and sequence content
    header = lines[0].strip()
    sequence = ''.join(line.strip() for line in lines[1:])
    
    # Calculate sequence length and ensure shift size doesn't exceed sequence length
    seq_length = len(sequence)
    shift_size = shift_size % seq_length
    
    # Perform sequence shift from the beginning
    shifted_sequence = sequence[shift_size:] + sequence[:shift_size]

    # Write shifted sequence to output file
    with open(output_file, 'w') as f:
        f.write(header + '\n')
        # Output bases with specified line width
        for i in range(0, len(shifted_sequence), line_width):
            f.write(shifted_sequence[i:i+line_width] + '\n')


def main():
    # Create command line argument parser
    parser = argparse.ArgumentParser(description='Perform shift operation on FASTA format sequence file')

    # Add command line arguments
    parser.add_argument('-i', '--input', required=True,
                        help='Input FASTA file path')
    parser.add_argument('-o', '--output', required=True,
                        help='Output FASTA file path')
    parser.add_argument('-s', '--shift', type=int, required=True,
                        help='Number of bases to shift')
    parser.add_argument('-w', '--width', type=int, default=70,
                        help='Number of bases per line in output file '
                             '(default: 70)')
                             
    # Parse command line arguments
    args = parser.parse_args()

    try:
        shift_fasta_sequence(args.input, args.output, args.shift, args.width)
        print(f"Sequence successfully shifted by {args.shift} bp and saved to {args.output}")
        print(f"Output file contains {args.width} bases per line")
    except Exception as e:
        print(f"Error occurred during processing: {str(e)}")


if __name__ == "__main__":
    main()
