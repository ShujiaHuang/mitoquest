#!/usr/bin/env python3
"""
Convert mtDNA VCF into tidy (long-format) table.

Each row corresponds to:
    one sample × one mtDNA position × one ALT allele

Output columns:
    sample_id, pos, ref, alt, vaf, depth, GT

"""
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Optional

import pysam


@dataclass
class VariantRecord:
    """Represents a single variant for one sample and one alt allele."""
    sample_id: str
    pos: int
    ref: str
    alt: str
    vaf: float
    depth: int
    genotype: str

    def to_tsv_line(self) -> str:
        """Convert record to TSV line."""
        return "\t".join([
            self.sample_id,
            str(self.pos),
            self.ref,
            self.alt,
            f"{self.vaf:.6f}",
            str(self.depth),
            self.genotype,
        ])


class VCFProcessor:
    """Processor for VCF files using pysam."""
    HEADER_COLUMNS = ["sample_id", "pos", "ref", "alt", "vaf", "depth", "GT"]
    def __init__(self, vcf_path: str):
        """
        Initialize VCF processor.
        
        Parameters
        ----------
        vcf_path : str
            Path to input VCF file
        """
        self.vcf_path = Path(vcf_path)
        self._validate_file()
    
    def _validate_file(self) -> None:
        """Validate that VCF file exists and is readable."""
        if not self.vcf_path.exists():
            raise FileNotFoundError(f"VCF file not found: {self.vcf_path}")
        
        if not self.vcf_path.is_file():
            raise ValueError(f"Not a file: {self.vcf_path}")
    
    def process(self) -> Iterator[VariantRecord]:
        """
        Process VCF file and yield variant records.
        
        Yields
        ------
        VariantRecord
            Tidy variant records
        
        Raises
        ------
        ValueError
            If VCF file has no samples
        """
        try:
            with pysam.VariantFile(str(self.vcf_path)) as vcf:
                samples = list(vcf.header.samples)
                
                if not samples:
                    raise ValueError("No samples found in VCF file")
                
                for record in vcf:
                    yield from self._process_record(record, samples)
        
        except (OSError, IOError) as e:
            raise IOError(f"Error reading VCF file: {e}")
        except Exception as e:
            raise RuntimeError(f"Error processing VCF: {e}")
    
    def _process_record(
        self, 
        record: pysam.VariantRecord, 
        samples: List[str]
    ) -> Iterator[VariantRecord]:
        """
        Process a single VCF record and yield variant records.
        
        Parameters
        ----------
        record : pysam.VariantRecord
            VCF record from pysam
        samples : List[str]
            List of sample names
        
        Yields
        ------
        VariantRecord
            Variant records for each sample and alt allele
        """
        pos = record.pos
        ref = record.ref
        alts = record.alts if record.alts else []
        
        for sample_id in samples:
            sample = record.samples[sample_id]
            
            # Skip missing genotypes
            if self._is_missing_genotype(sample):
                continue
            
            genotype = self._extract_genotype(sample)
            depth = self._extract_depth(sample)
            vaf_values = self._extract_vaf_values(sample, len(alts))
            
            for alt_idx, alt in enumerate(alts):
                vaf = vaf_values[alt_idx] if alt_idx < len(vaf_values) else 0.0
                
                yield VariantRecord(
                    sample_id=sample_id,
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    vaf=vaf,
                    depth=depth,
                    genotype=genotype,
                )
    
    @staticmethod
    def _is_missing_genotype(sample: pysam.VariantRecordSample) -> bool:
        """
        Check if sample has missing genotype.
        
        Parameters
        ----------
        sample : pysam.VariantRecordSample
            Sample from VCF record
        
        Returns
        -------
        bool
            True if genotype is missing or None
        """
        try:
            gt = sample.get('GT')
            if gt is None:
                return True
            
            # Check if all alleles are None (missing)
            if isinstance(gt, tuple):
                return all(allele is None for allele in gt)
            
            return False
        except (KeyError, AttributeError):
            return True
    
    @staticmethod
    def _extract_genotype(sample: pysam.VariantRecordSample) -> str:
        """
        Extract genotype string from sample.
        
        Parameters
        ----------
        sample : pysam.VariantRecordSample
            Sample from VCF record
        
        Returns
        -------
        str
            Genotype string (e.g., '0/1', '1/1') or '.' if missing
        """
        try:
            gt = sample.get('GT')
            if gt is None:
                return '.'
            
            # pysam returns GT as tuple of allele indices
            if isinstance(gt, tuple):
                # Get phasing character
                phased = sample.phased
                sep = '|' if phased else '/'
                
                # Convert None to '.' for missing alleles
                alleles = [str(a) if a is not None else '.' for a in gt]
                return sep.join(alleles)
            
            return str(gt)
        except (KeyError, AttributeError):
            return '.'
    
    @staticmethod
    def _extract_depth(sample: pysam.VariantRecordSample) -> int:
        """
        Extract depth from sample.
        
        Parameters
        ----------
        sample : pysam.VariantRecordSample
            Sample from VCF record
        
        Returns
        -------
        int
            Read depth, or 0 if missing
        """
        try:
            dp = sample.get('DP')
            if dp is None:
                return 0
            return int(dp)
        except (KeyError, ValueError, AttributeError):
            return 0
    
    @staticmethod
    def _extract_vaf_values(
        sample: pysam.VariantRecordSample, 
        num_alts: int
    ) -> List[float]:
        """
        Extract VAF values from sample.
        
        Parameters
        ----------
        sample : pysam.VariantRecordSample
            Sample from VCF record
        num_alts : int
            Number of alternate alleles
        
        Returns
        -------
        List[float]
            List of VAF values for each alt allele
        """
        try:
            af = sample.get('AF')
            
            if af is None:
                return [0.0] * num_alts
            
            # AF can be a single value or tuple
            if isinstance(af, (tuple, list)):
                return [float(v) if v is not None else 0.0 for v in af]
            else:
                # Single value case
                return [float(af) if af is not None else 0.0]
        
        except (KeyError, ValueError, AttributeError):
            return [0.0] * num_alts


def write_tidy_table(
    records: Iterator[VariantRecord], 
    output_path: str
) -> int:
    """
    Write variant records to tidy TSV file.
    
    Parameters
    ----------
    records : Iterator[VariantRecord]
        Iterator of variant records
    output_path : str
        Path to output TSV file
    
    Returns
    -------
    int
        Number of records written
    
    Raises
    ------
    PermissionError
        If output file cannot be written
    OSError
        If there's an error writing the file
    """
    output_file = Path(output_path)
    count = 0
    
    try:
        with open(output_file, 'w') as f:
            # Write header
            f.write("\t".join(VCFProcessor.HEADER_COLUMNS) + "\n")
            
            # Write records
            for record in records:
                if record.genotype == "0": 
                    continue
                
                f.write(record.to_tsv_line() + "\n")
                count += 1
        
        return count
    
    except PermissionError:
        raise PermissionError(f"Permission denied: {output_file}")
    except OSError as e:
        raise OSError(f"Error writing to {output_file}: {e}")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert mtDNA VCF into tidy long-format table using pysam.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input-vcf",
        required=True,
        help="Input mtDNA VCF file (can be .vcf, .vcf.gz, or .bcf)",
    )
    parser.add_argument(
        "-o",
        "--output-tsv",
        required=True,
        help="Output tidy TSV file",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print progress information",
    )
    return parser.parse_args()


def main() -> None:
    """Main entry point."""
    args = parse_args()
    try:
        if args.verbose:
            print(f"Processing VCF file: {args.input_vcf}", file=sys.stderr)
        
        processor = VCFProcessor(args.input_vcf)
        records = processor.process()
        count = write_tidy_table(records, args.output_tsv)
        
        if args.verbose:
            print(
                f"Successfully wrote {count} records to {args.output_tsv}", 
                file=sys.stderr
            )
    
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except PermissionError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except (IOError, OSError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nInterrupted by user", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
