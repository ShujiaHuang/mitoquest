#!/usr/bin/env python3
"""
Convert mtDNA VCF into tidy (long-format) table.

Each row corresponds to:
    one sample x one mtDNA position x one ALT allele

Output columns:
    sample_id, pos, ref, alt, vaf, depth, GT

Authors: Shujia Huang
Date: 2026-01-24
"""
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Optional

import numpy as np
import pysam

@dataclass
class VariantRecord:
    """Represents a single variant for one sample and one alt allele."""
    sample_id: str
    chrom: str
    pos: int
    rsid: str
    ref: str
    alt: str
    vaf: float
    srf: float
    depth: int
    genotype: str
    var_type: str # The type of variant: SNV, DEL, INS
    status: str   # HET or HOM

    def to_tsv_line(self) -> str:
        """Convert record to TSV line."""
        return "\t".join([
            self.sample_id,
            self.chrom,
            str(self.pos),
            self.rsid,
            self.ref,
            self.alt,
            f"{self.vaf:.6f}",
            f"{self.srf:.6f}",
            str(self.depth),
            self.genotype,
            self.var_type,
            self.status
        ])


class VCFProcessor:
    """Processor for VCF files using pysam."""
    HEADER_COLUMNS = ["Sample_id", "Chrom", "Pos", "ID", "REF", "ALT", "VAF", "SRF", 
                      "Depth", "GT", "Variant", "Status"]
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
    
    def remove_common_suffix(self, s1, s2):
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
        chrom = record.chrom
        pos   = record.pos
        ref   = record.ref if record.ref is not None else '.'
        alts  = record.alts if record.alts else []
        rsid  = record.id if record.id is not None else '.'
        for sample_id in samples:
            sample = record.samples[sample_id]
            
            # Skip missing genotypes
            if self._is_missing_genotype(sample):
                continue

            # Coverts string to boolean, default to False if not present 
            is_good_call = sample.get('GOOD_CALL', 'False') == 'True'
            if not is_good_call:
                continue
            
            # pysam doesn't support directly getting a tuple of tuples for the specific label of SB (format problem), 
            # so we need to parse it from the string format. Set back to be a string of format in 'a1_fwd,a1_reverse;a2_forward,a2_reverse;...' 
            # for each sample, where a1, a2, ... are ref and alt alleles in order.
            sb_s = ','.join(sample.get('SB')) if sample.get('SB') else '.'
            sb_split = sb_s.split(';') if sb_s != '.' else []
            
            # A list of tuples [(a1_forward, a1_reverse), (a2_forward, a2_reverse), ...]
            sbs = [tuple(map(int, sb.split(','))) for sb in sb_split] # could be empthy if SB is missing
            
            # Strand ratio factor can be calculated as np.min([fwd, rev])/np.max([fwd, rev]) for each allele.
            # The range of SRF is [0, 1], where values close to 0 indicate strong strand bias and values 
            # close to 1 indicate balanced strand representation.
            srfs = [(np.min([fwd, rev]) + 1e-10)/(np.max([fwd, rev]) + 1e-10) for fwd, rev in sbs]
            
            gts = sample.get('GT')
            gts_non_missing = [gt for gt in gts if gt is not None] if gts is not None else []
            vaf_values = sample.get('AF')
            depth = self._extract_depth(sample)
            for i, (gt, vaf) in enumerate(zip(gts, vaf_values)):
                if gt is None:
                    var_type = "MISSING"
                    new_ref = ref
                    alt_seq = '.'
                elif gt == 0:
                    var_type = "REF"
                    new_ref = ref
                    alt_seq = ref
                else:
                    alt_seq  = alts[gt - 1]
                    if len(ref) == len(alt_seq):
                        var_type = "SNV"
                    elif len(ref) > len(alt_seq):
                        var_type = "DEL"
                    elif len(ref) < len(alt_seq):
                        var_type = "INS"
                    else:
                        var_type = "UNK"
                        
                    new_ref, alt_seq = self.remove_common_suffix(ref, alt_seq)
                yield VariantRecord(
                    sample_id=sample_id,
                    chrom=chrom,
                    pos=pos,
                    rsid=rsid,
                    ref=new_ref,
                    alt=alt_seq,
                    vaf=vaf,
                    srf=srfs[i] if i < len(srfs) else 1.0,  # Default to 1.0 (no bias) if SB info is missing
                    depth=depth,
                    genotype="/".join(map(str, gts)),
                    var_type=var_type,
                    status="HET" if len(gts_non_missing) > 1 else "HOM",
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
