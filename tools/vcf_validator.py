#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VCF Format Validator.

This module provides a comprehensive validator for VCF (Variant Call Format) files.
It checks the validity of headers, INFO fields, FORMAT fields, and their
relationships in the VCF file.

Author: Shujia Huang
Date: 2026-01-15
"""
import sys
import re
import argparse
import gzip
from typing import Dict, List, Set, Tuple, Optional


class VCFValidator:
    """
    A comprehensive validator for VCF files.
    
    This class validates VCF file format including:
    - File format version
    - Header lines
    - INFO field definitions and usage
    - FORMAT field definitions and usage
    - Column structure
    - Data type consistency
    """
    
    def __init__(self, filepath: str, verbose: bool = False):
        """
        Initialize the VCF validator.
        
        Args:
            filepath: Path to the VCF file to validate
            verbose: Whether to print verbose output
        """
        self.filepath = filepath
        self.verbose = verbose
        self.errors: List[str] = []
        self.warnings: List[str] = []
        
        # Store header definitions
        self.info_definitions: Dict[str, Dict] = {}
        self.format_definitions: Dict[str, Dict] = {}
        self.filter_definitions: Set[str] = set()
        self.contig_definitions: Set[str] = set()
        self.alt_definitions: Dict[str, str] = {}
        
        # Column headers
        self.column_headers: List[str] = []
        self.sample_names: List[str] = []
        
        # VCF version
        self.vcf_version: Optional[str] = None
        
        # Valid VCF data types
        self.valid_types = {'Integer', 'Float', 'String', 'Character', 'Flag'}
        self.valid_numbers = {'A', 'R', 'G', '.'}  # Or positive integers
        
    def add_error(self, line_num: int, message: str) -> None:
        """Add an error message."""
        self.errors.append(f"Line {line_num}: ERROR - {message}")
        
    def add_warning(self, line_num: int, message: str) -> None:
        """Add a warning message."""
        self.warnings.append(f"Line {line_num}: WARNING - {message}")
        
    def log(self, message: str) -> None:
        """Print message if verbose mode is enabled."""
        if self.verbose:
            print(f"[INFO] {message}")
            
    def parse_meta_info(self, line: str) -> Optional[Tuple[str, Dict]]:
        """
        Parse meta-information lines like INFO, FORMAT, FILTER, etc.
        
        Args:
            line: The meta-information line
            
        Returns:
            Tuple of (field_id, metadata_dict) or None if parsing fails
        """
        # Pattern to match key=value pairs within <>
        pattern = r'<(.+)>'
        match = re.search(pattern, line)
        
        if not match:
            return None
            
        content = match.group(1)
        metadata = {}
        
        # Parse key=value pairs, handling quoted values
        current_key = None
        current_value = []
        in_quotes = False
        i = 0
        
        while i < len(content):
            char = content[i]
            
            if char == '"':
                in_quotes = not in_quotes
                i += 1
                continue
                
            if char == '=' and not in_quotes and current_key is None:
                # Found key
                key = ''.join(current_value).strip()
                current_key = key
                current_value = []
                i += 1
                continue
                
            if char == ',' and not in_quotes and current_key is not None:
                # Found end of value
                value = ''.join(current_value).strip()
                metadata[current_key] = value
                current_key = None
                current_value = []
                i += 1
                continue
                
            current_value.append(char)
            i += 1
            
        # Handle last key-value pair
        if current_key is not None:
            value = ''.join(current_value).strip()
            metadata[current_key] = value
            
        # Ensure we return either a (str, dict) tuple or None.
        field_id = metadata.get('ID')
        if field_id is None:
            return None

        return field_id, metadata
        
    def validate_header_line(self, line_num: int, line: str) -> None:
        """
        Validate a header line.
        
        Args:
            line_num: Line number in the file
            line: The header line content
        """
        line = line.strip()
        
        # Check for fileformat
        if line.startswith('##fileformat='):
            self.vcf_version = line.split('=', 1)[1]
            self.log(f"VCF version: {self.vcf_version}")
            
            if not self.vcf_version.startswith('VCFv'):
                self.add_error(line_num, 
                             f"Invalid fileformat: {self.vcf_version}")
                             
        # Parse INFO definitions
        elif line.startswith('##INFO='):
            result = self.parse_meta_info(line)
            if result:
                field_id, metadata = result
                
                # Validate required fields
                if 'Number' not in metadata:
                    self.add_error(line_num, 
                                 f"INFO {field_id} missing 'Number' field")
                elif not self._is_valid_number(metadata['Number']):
                    self.add_error(line_num,
                                 f"INFO {field_id} has invalid Number: "
                                 f"{metadata['Number']}")
                    
                if 'Type' not in metadata:
                    self.add_error(line_num,
                                 f"INFO {field_id} missing 'Type' field")
                elif metadata['Type'] not in self.valid_types:
                    self.add_error(line_num,
                                 f"INFO {field_id} has invalid Type: "
                                 f"{metadata['Type']}")
                    
                if 'Description' not in metadata:
                    self.add_warning(line_num,
                                   f"INFO {field_id} missing 'Description' field")
                    
                self.info_definitions[field_id] = metadata
                self.log(f"Registered INFO field: {field_id}")
                
        # Parse FORMAT definitions
        elif line.startswith('##FORMAT='):
            result = self.parse_meta_info(line)
            if result:
                field_id, metadata = result
                
                # Validate required fields
                if 'Number' not in metadata:
                    self.add_error(line_num,
                                 f"FORMAT {field_id} missing 'Number' field")
                elif not self._is_valid_number(metadata['Number']):
                    self.add_error(line_num,
                                 f"FORMAT {field_id} has invalid Number: "
                                 f"{metadata['Number']}")
                    
                if 'Type' not in metadata:
                    self.add_error(line_num,
                                 f"FORMAT {field_id} missing 'Type' field")
                elif metadata['Type'] not in self.valid_types:
                    self.add_error(line_num,
                                 f"FORMAT {field_id} has invalid Type: "
                                 f"{metadata['Type']}")
                    
                if 'Description' not in metadata:
                    self.add_warning(line_num,
                                   f"FORMAT {field_id} missing 'Description' field")
                    
                self.format_definitions[field_id] = metadata
                self.log(f"Registered FORMAT field: {field_id}")
                
        # Parse FILTER definitions
        elif line.startswith('##FILTER='):
            result = self.parse_meta_info(line)
            if result:
                field_id, metadata = result
                self.filter_definitions.add(field_id)
                self.log(f"Registered FILTER: {field_id}")
                
        # Parse contig definitions
        elif line.startswith('##contig='):
            result = self.parse_meta_info(line)
            if result:
                field_id, metadata = result
                self.contig_definitions.add(field_id)
                self.log(f"Registered contig: {field_id}")
                
        # Parse ALT definitions
        elif line.startswith('##ALT='):
            result = self.parse_meta_info(line)
            if result:
                field_id, metadata = result
                self.alt_definitions[field_id] = metadata.get('Description', '')
                self.log(f"Registered ALT: {field_id}")
                
    def _is_valid_number(self, number: str) -> bool:
        """
        Check if a Number field value is valid.
        
        Args:
            number: The Number field value
            
        Returns:
            True if valid, False otherwise
        """
        if number in self.valid_numbers:
            return True
            
        # Check if it's a positive integer
        try:
            num = int(number)
            return num >= 0
        except ValueError:
            return False
            
    def validate_column_header(self, line_num: int, line: str) -> None:
        """
        Validate the column header line.
        
        Args:
            line_num: Line number in the file
            line: The column header line
        """
        line = line.strip()
        if not line.startswith('#CHROM'):
            self.add_error(line_num, "Column header must start with #CHROM")
            return
            
        columns = line.split('\t')
        self.column_headers = columns
        
        # Required columns
        required = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 
                   'QUAL', 'FILTER', 'INFO']
        
        if len(columns) < 8:
            self.add_error(line_num,
                         f"Insufficient columns. Expected at least 8, "
                         f"got {len(columns)}")
            return
            
        for i, req in enumerate(required):
            if i >= len(columns) or columns[i] != req:
                self.add_error(line_num,
                             f"Column {i+1} should be '{req}', "
                             f"got '{columns[i] if i < len(columns) else 'missing'}'")
                             
        # If FORMAT column exists, there should be at least one sample
        if len(columns) > 8:
            if columns[8] != 'FORMAT':
                self.add_error(line_num,
                             f"Column 9 should be 'FORMAT', got '{columns[8]}'")
            else:
                self.sample_names = columns[9:]
                self.log(f"Found {len(self.sample_names)} sample(s): "
                        f"{', '.join(self.sample_names[:5])}"
                        f"{'...' if len(self.sample_names) > 5 else ''}")
                        
    def validate_data_line(self, line_num: int, line: str) -> None:
        """
        Validate a data line.
        
        Args:
            line_num: Line number in the file
            line: The data line content
        """
        line = line.strip()
        fields = line.split('\t')
        
        # Check number of fields
        expected_fields = len(self.column_headers)
        if len(fields) != expected_fields:
            self.add_error(line_num,
                         f"Expected {expected_fields} fields, "
                         f"got {len(fields)}")
            return
            
        # Validate CHROM
        chrom = fields[0]
        if self.contig_definitions and chrom not in self.contig_definitions:
            self.add_warning(line_num,
                           f"CHROM '{chrom}' not defined in ##contig headers")
            
        # Validate POS
        try:
            pos = int(fields[1])
            if pos < 1:
                self.add_error(line_num,
                             f"POS must be positive, got {pos}")
        except ValueError:
            self.add_error(line_num,
                         f"POS must be an integer, got '{fields[1]}'")
            
        # Validate REF
        ref = fields[3]
        if not re.match(r'^[ACGTNacgtn]+$', ref):
            self.add_error(line_num,
                         f"REF contains invalid characters: '{ref}'")
            
        # Validate ALT
        alt = fields[4]
        if alt != '.':
            alts = alt.split(',')
            for alt_allele in alts:
                if not re.match(r'^[ACGTNacgtn]+$|^<.+>$|^\*$', alt_allele):
                    self.add_error(line_num,
                                 f"ALT contains invalid allele: '{alt_allele}'")
                    
        # Validate QUAL
        qual = fields[5]
        if qual != '.':
            try:
                float(qual)
            except ValueError:
                self.add_error(line_num,
                             f"QUAL must be a number or '.', got '{qual}'")
                
        # Validate FILTER
        filter_val = fields[6]
        if filter_val not in ['.', 'PASS']:
            filters = filter_val.split(';')
            for filt in filters:
                if filt not in self.filter_definitions:
                    self.add_warning(line_num,
                                   f"FILTER '{filt}' not defined in headers")
                    
        # Validate INFO
        info = fields[7]
        if info != '.':
            self._validate_info_field(line_num, info)
            
        # Validate FORMAT and sample fields
        if len(fields) > 8:
            format_field = fields[8]
            sample_fields = fields[9:]
            self._validate_format_fields(line_num, format_field, sample_fields)
            
    def _validate_info_field(self, line_num: int, info: str) -> None:
        """
        Validate INFO field.
        
        Args:
            line_num: Line number in the file
            info: The INFO field content
        """
        info_items = info.split(';')
        
        for item in info_items:
            if '=' in item:
                key, value = item.split('=', 1)
            else:
                key = item
                value = None
                
            # Check if INFO field is defined
            if key not in self.info_definitions:
                self.add_warning(line_num,
                               f"INFO field '{key}' not defined in headers")
                continue
                
            definition = self.info_definitions[key]
            
            # Flag type should not have a value
            if definition.get('Type') == 'Flag':
                if value is not None:
                    self.add_error(line_num,
                                 f"INFO field '{key}' is a Flag but has a value")
            else:
                if value is None:
                    self.add_error(line_num,
                                 f"INFO field '{key}' requires a value")
                else:
                    # Validate value type
                    self._validate_value_type(line_num, f"INFO {key}",
                                            value, definition)
                                            
    def _validate_format_fields(self, line_num: int, format_field: str,
                                sample_fields: List[str]) -> None:
        """
        Validate FORMAT and sample fields.
        
        Args:
            line_num: Line number in the file
            format_field: The FORMAT field content
            sample_fields: List of sample field contents
        """
        format_keys = format_field.split(':')
        
        # Validate format keys are defined
        for key in format_keys:
            if key not in self.format_definitions:
                self.add_warning(line_num,
                               f"FORMAT field '{key}' not defined in headers")
                
        # Validate each sample field
        for i, sample in enumerate(sample_fields):
            sample_values = sample.split(':')
            
            if len(sample_values) != len(format_keys):
                self.add_error(line_num,
                             f"Sample {self.sample_names[i] if i < len(self.sample_names) else i+1} "
                             f"has {len(sample_values)} values but FORMAT has "
                             f"{len(format_keys)} keys")
                continue
                
            # Validate each value
            for key, value in zip(format_keys, sample_values):
                if key in self.format_definitions and value != '.':
                    definition = self.format_definitions[key]
                    self._validate_value_type(line_num,
                                            f"FORMAT {key} in sample "
                                            f"{self.sample_names[i] if i < len(self.sample_names) else i+1}",
                                            value, definition)
                                            
    def _validate_value_type(self, line_num: int, field_name: str,
                            value: str, definition: Dict) -> None:
        """
        Validate value against type definition.
        
        Args:
            line_num: Line number in the file
            field_name: Name of the field being validated
            value: The value to validate
            definition: The field definition containing Type and Number
        """
        value_type = definition.get('Type')
        
        # Split by comma for multiple values
        values = value.split(',')
        
        for val in values:
            if val == '.':
                continue
                
            if value_type == 'Integer':
                try:
                    int(val)
                except ValueError:
                    self.add_error(line_num,
                                 f"{field_name} should be Integer, got '{val}'")
                    
            elif value_type == 'Float':
                try:
                    float(val)
                except ValueError:
                    self.add_error(line_num,
                                 f"{field_name} should be Float, got '{val}'")
                    
            elif value_type == 'Character':
                if len(val) != 1:
                    self.add_error(line_num,
                                 f"{field_name} should be a single character, "
                                 f"got '{val}'")
                                 
    def validate(self) -> bool:
        """
        Validate the VCF file.
        
        Returns:
            True if validation passed, False otherwise
        """
        self.log(f"Starting validation of {self.filepath}")
        
        try:
            with gzip.open(self.filepath, 'rt', encoding='utf-8') if self.filepath.endswith('.gz') else open(self.filepath, 'r', encoding='utf-8') as f:
                line_num = 0
                in_header = True
                
                for line in f:
                    line_num += 1
                    
                    # Skip empty lines
                    if not line.strip():
                        continue
                        
                    # Header lines
                    if line.startswith('##'):
                        self.validate_header_line(line_num, line)
                        
                    # Column header
                    elif line.startswith('#'):
                        self.validate_column_header(line_num, line)
                        in_header = False
                        
                    # Data lines
                    else:
                        if in_header:
                            self.add_error(line_num,
                                         "Data line found before column header")
                        self.validate_data_line(line_num, line)
                        
        except FileNotFoundError:
            self.errors.append(f"File not found: {self.filepath}")
            return False
        except Exception as e:
            self.errors.append(f"Unexpected error: {str(e)}")
            return False
            
        # Check for required fileformat
        if self.vcf_version is None:
            self.errors.insert(0, "Missing required ##fileformat header")
            
        # Print results
        self._print_results()
        
        return len(self.errors) == 0
        
    def _print_results(self) -> None:
        """Print validation results."""
        print("\n" + "=" * 70)
        print(f"VCF Validation Report for: {self.filepath}")
        print("=" * 70)
        
        if self.vcf_version:
            print(f"\nVCF Version: {self.vcf_version}")
            
        print(f"\nINFO fields defined: {len(self.info_definitions)}")
        print(f"FORMAT fields defined: {len(self.format_definitions)}")
        print(f"FILTER fields defined: {len(self.filter_definitions)}")
        print(f"Samples: {len(self.sample_names)}")
        
        if self.warnings:
            print(f"\n{len(self.warnings)} Warning(s):")
            print("-" * 70)
            for warning in self.warnings[:10]:  # Show first 10 warnings
                print(warning)
            if len(self.warnings) > 10:
                print(f"... and {len(self.warnings) - 10} more warnings")
                
        if self.errors:
            print(f"\n{len(self.errors)} Error(s):")
            print("-" * 70)
            for error in self.errors[:10]:  # Show first 10 errors
                print(error)
            if len(self.errors) > 10:
                print(f"... and {len(self.errors) - 10} more errors")
                
        print("\n" + "=" * 70)
        
        if self.errors:
            print("VALIDATION FAILED ✗")
        else:
            print("VALIDATION PASSED ✓")
            
        print("=" * 70 + "\n")


def main():
    """Main function to run the VCF validator from command line."""
    parser = argparse.ArgumentParser(
        description='Validate VCF (Variant Call Format) files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s input.vcf
  %(prog)s input.vcf --verbose
  %(prog)s input.vcf -v -o report.txt
        """
    )
    
    parser.add_argument(
        'input',
        type=str,
        help='Input VCF file to validate'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        default=None,
        help='Output file for validation report (default: stdout)'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 1.0.0'
    )
    
    args = parser.parse_args()
    
    # Redirect output if specified
    if args.output:
        sys.stdout = open(args.output, 'w', encoding='utf-8')
        
    try:
        validator = VCFValidator(args.input, verbose=args.verbose)
        success = validator.validate()
        
        # Return appropriate exit code
        sys.exit(0 if success else 1)
        
    except KeyboardInterrupt:
        print("\nValidation interrupted by user")
        sys.exit(130)
    except Exception as e:
        print(f"\nFatal error: {str(e)}", file=sys.stderr)
        sys.exit(1)
    finally:
        if args.output:
            sys.stdout.close()


if __name__ == '__main__':
    main()
