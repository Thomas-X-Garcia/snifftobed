#!/usr/bin/env python3
"""
snifftobed.py - Universal VCF to BED converter with three output formats

This script converts structural variant VCF files to BED format with three output options:
1. Default: Standard BED format compatible with bedtools
2. --enhanced: Standard BED + all original SURVIVOR columns + VCF fields
3. --original: Byte-for-byte identical to SURVIVOR's vcftobed

Author: Thomas X. Garcia, PhD, HCLD
License: MIT
Version: 1.0
"""

import argparse
import sys
import os
import re
from typing import Dict, List, Tuple, Optional, NamedTuple
from dataclasses import dataclass
from datetime import datetime


class SVType:
    """Enumeration of structural variant types"""
    DEL = 0
    DUP = 1
    INV = 2
    TRA = 3
    INS = 4
    BND = 5
    CNV = 6
    UNKNOWN = -1


@dataclass
class Coordinate:
    """Represents a genomic coordinate"""
    chr: str
    pos: int


@dataclass
class StrandPair:
    """Represents strand information for breakpoints"""
    first: bool  # True = forward (+), False = reverse (-)
    second: bool


@dataclass
class SVEntry:
    """Represents a structural variant entry"""
    start: Coordinate
    stop: Coordinate
    sv_type: int
    sv_len: int
    sv_id: str
    strands: StrandPair
    cpos: Tuple[int, int]  # CIPOS
    cend: Tuple[int, int]  # CIEND
    quality: float
    genotype: str
    header: str
    # Additional fields for enhanced mode
    ref: str = ""
    alt: str = ""
    qual: str = "."
    filter: str = "."
    info: str = ""
    format: str = ""
    samples: List[str] = None


class VCFParser:
    """Parser for VCF files following SURVIVOR's implementation"""
    
    def __init__(self, min_size: int = 0):
        self.min_size = min_size
        self.entries_parsed = 0
        self.entries_kept = 0
        self.vcf_headers = []
        self.column_headers = []
        self.sample_names = []
        
    def get_type(self, type_str: str) -> int:
        """Convert SV type string to numeric code"""
        type_upper = type_str.upper()
        
        if type_upper.startswith("DEL"):
            return SVType.DEL
        elif type_upper.startswith("DUP"):
            return SVType.DUP
        elif type_upper.startswith("INV"):
            return SVType.INV
        elif type_upper.startswith("TRA"):
            return SVType.TRA
        elif (type_upper.startswith("INS") or type_upper.startswith("ALU") or 
              type_upper.startswith("LINE1") or type_upper.startswith("SVA")):
            return SVType.INS
        elif type_upper.startswith("BND"):
            return SVType.BND
        elif type_upper.startswith("CNV"):
            return SVType.CNV
        return SVType.UNKNOWN
    
    def trans_type(self, sv_type: int) -> str:
        """Convert numeric SV type to string"""
        type_map = {
            SVType.DEL: "DEL",
            SVType.DUP: "DUP",
            SVType.INV: "INV",
            SVType.TRA: "TRA",
            SVType.INS: "INS",
            SVType.BND: "BND",
            SVType.CNV: "CNV"
        }
        return type_map.get(sv_type, "NA")
    
    def parse_stop(self, info_field: str) -> Coordinate:
        """Parse END position and CHR2 from INFO field"""
        stop = Coordinate(chr="", pos=-1)
        
        # Parse END position
        end_match = re.search(r'END=(\d+)', info_field)
        if end_match:
            stop.pos = int(end_match.group(1))
        
        # Parse CHR2 for translocations
        chr2_match = re.search(r'CHR2=([^;]+)', info_field)
        if chr2_match:
            stop.chr = chr2_match.group(1)
            
        return stop
    
    def parse_strands_lumpy(self, alt_field: str) -> StrandPair:
        """Parse strand information from ALT field (Lumpy style)"""
        strands = StrandPair(first=True, second=True)
        
        bracket_count = 0
        for char in alt_field:
            if char == '[':
                if bracket_count == 0:
                    strands.first = False
                else:
                    strands.second = False
                bracket_count += 1
            elif char == ']':
                if bracket_count == 0:
                    strands.first = True
                else:
                    strands.second = True
                bracket_count += 1
            elif char == '\t':
                break
                
        return strands
    
    def parse_pos_from_alt(self, alt_field: str) -> Optional[Coordinate]:
        """Parse position from ALT field for BND entries (e.g., N[chr2:12345[)"""
        # Patterns: N[chr2:12345[ or ]chr2:12345]N or [chr2:12345[N or N]chr2:12345]
        match = re.search(r'[\[\]]([^:]+):(\d+)[\[\]]', alt_field)
        if match:
            return Coordinate(chr=match.group(1), pos=int(match.group(2)))
        return None
    
    def parse_strands(self, info_field: str) -> Optional[StrandPair]:
        """Parse strand information from INFO field"""
        # Check for CT field (Delly style)
        ct_match = re.search(r'CT=(\d)to(\d)', info_field)
        if ct_match:
            first = ct_match.group(1) != '5'
            second = ct_match.group(2) != '5'
            return StrandPair(first=first, second=second)
        
        # Check for STRANDS field
        strands_match = re.search(r'STRANDS=([+-])([+-])', info_field)
        if strands_match:
            first = strands_match.group(1) == '+'
            second = strands_match.group(2) == '+'
            return StrandPair(first=first, second=second)
        
        return None
    
    def parse_pair(self, value_str: str) -> Tuple[int, int]:
        """Parse comma-separated integer pair (e.g., CIPOS=-10,10)"""
        parts = value_str.split(',')
        if len(parts) >= 2:
            try:
                return (int(parts[0]), int(parts[1]))
            except ValueError:
                pass
        return (0, 0)
    
    def parse_vcf_entry(self, line: str) -> Optional[SVEntry]:
        """Parse a single VCF entry line"""
        if line.startswith('#') or not line.strip():
            return None
            
        fields = line.strip().split('\t')
        if len(fields) < 8:
            return None
        
        # Initialize entry
        entry = SVEntry(
            start=Coordinate(chr=fields[0], pos=int(fields[1])),
            stop=Coordinate(chr="", pos=-1),
            sv_type=SVType.UNKNOWN,
            sv_len=-1,
            sv_id=fields[2] if len(fields) > 2 else "",
            strands=StrandPair(first=True, second=True),
            cpos=(0, 0),
            cend=(0, 0),
            quality=float(fields[5]) if len(fields) > 5 and fields[5] != '.' else -1,
            genotype="./.",
            header="",
            # Additional fields for enhanced mode
            ref=fields[3] if len(fields) > 3 else "",
            alt=fields[4] if len(fields) > 4 else "",
            qual=fields[5] if len(fields) > 5 else ".",
            filter=fields[6] if len(fields) > 6 else ".",
            info=fields[7] if len(fields) > 7 else "",
            format=fields[8] if len(fields) > 8 else "",
            samples=fields[9:] if len(fields) > 9 else []
        )
        
        # Get REF and ALT
        ref = fields[3] if len(fields) > 3 else ""
        alt = fields[4] if len(fields) > 4 else ""
        
        # Parse strands from ALT field (Lumpy style)
        entry.strands = self.parse_strands_lumpy(alt)
        
        # Parse INFO field
        info = fields[7] if len(fields) > 7 else ""
        
        # Parse stop position
        entry.stop = self.parse_stop(info)
        
        # Parse SV type
        svtype_match = re.search(r'SVTYPE=([^;]+)', info)
        if svtype_match:
            entry.sv_type = self.get_type(svtype_match.group(1))
        
        # For BND entries, also check ALT field for position
        if entry.sv_type == SVType.BND and (entry.stop.pos == -1 or not entry.stop.chr):
            alt_pos = self.parse_pos_from_alt(alt)
            if alt_pos:
                entry.stop = alt_pos
        
        # Parse SV length
        svlen_match = re.search(r'SVLEN=(-?\d+)', info)
        if svlen_match:
            entry.sv_len = abs(int(svlen_match.group(1)))
        
        # Parse CIPOS and CIEND
        cipos_match = re.search(r'CIPOS=([^;]+)', info)
        if cipos_match:
            entry.cpos = self.parse_pair(cipos_match.group(1))
        else:
            # Check if entry has implicit CIPOS (some tools may not output it)
            entry.cpos = (0, 0)
            
        ciend_match = re.search(r'CIEND=([^;]+)', info)
        if ciend_match:
            entry.cend = self.parse_pair(ciend_match.group(1))
        else:
            # Check if entry has implicit CIEND
            entry.cend = (0, 0)
        
        # Parse strand information from INFO field
        info_strands = self.parse_strands(info)
        if info_strands:
            entry.strands = info_strands
        
        # Special handling for BND/TRA - SURVIVOR seems to default to ++ for BND
        if entry.sv_type == SVType.BND and info_strands is None:
            entry.strands = StrandPair(first=True, second=True)
        
        # Set default strands based on SV type if not already set
        if info_strands is None:
            if entry.sv_type == SVType.DEL or entry.sv_type == SVType.INS:
                entry.strands = StrandPair(first=True, second=False)
            elif entry.sv_type == SVType.DUP:
                entry.strands = StrandPair(first=False, second=True)
        
        # Handle missing stop position
        if entry.stop.pos == -1 and entry.sv_len != -1:
            entry.stop.pos = entry.start.pos + entry.sv_len
        
        # Handle missing stop chromosome
        if not entry.stop.chr:
            entry.stop.chr = entry.start.chr
        
        # Calculate SV length if not provided
        if entry.sv_len == -1:
            if entry.sv_type == SVType.INS or entry.sv_type == SVType.DEL:
                # Try to calculate from REF/ALT
                entry.sv_len = abs(len(ref) - len(alt))
            else:
                entry.sv_len = abs(entry.start.pos - entry.stop.pos)
        
        # Handle BND type - convert to TRA or INV based on chromosomes
        if entry.sv_type == SVType.BND:
            if entry.start.chr == entry.stop.chr:
                entry.sv_type = SVType.INV
            else:
                entry.sv_type = SVType.TRA
        
        # Build header (first 9 fields)
        entry.header = '\t'.join(fields[:9]) if len(fields) >= 9 else '\t'.join(fields)
        
        # Parse genotype if available
        if len(fields) > 9:
            gt_field = fields[9].split(':')[0]
            if gt_field:
                entry.genotype = gt_field
        
        self.entries_parsed += 1
        return entry
    
    def parse_vcf(self, filename: str, min_size: int) -> List[SVEntry]:
        """Parse VCF file and return list of SV entries"""
        entries = []
        line_count = 0
        
        print(f"[{datetime.now().strftime('%H:%M:%S')}] Parsing VCF file: {filename}")
        
        try:
            with open(filename, 'r') as f:
                for line in f:
                    line_count += 1
                    
                    if line.startswith('##'):
                        self.vcf_headers.append(line.strip())
                        continue
                    
                    if line.startswith('#'):
                        # Parse column headers
                        self.column_headers = line.strip()[1:].split('\t')
                        if len(self.column_headers) > 9:
                            self.sample_names = self.column_headers[9:]
                        continue
                    
                    entry = self.parse_vcf_entry(line)
                    if entry:
                        # Check size filter
                        size = entry.sv_len
                        if entry.sv_type == SVType.INS:
                            size = entry.sv_len
                        elif entry.sv_type not in [SVType.TRA, SVType.BND, SVType.UNKNOWN]:
                            size = abs(entry.stop.pos - entry.start.pos)
                        
                        if size >= min_size or min_size == -1:
                            entries.append(entry)
                            self.entries_kept += 1
                    
                    if line_count % 10000 == 0:
                        print(f"[{datetime.now().strftime('%H:%M:%S')}] Processed {line_count} lines, found {len(entries)} SVs")
        
        except FileNotFoundError:
            print(f"ERROR: Could not open file: {filename}", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"ERROR: Failed to parse VCF file at line {line_count}: {str(e)}", file=sys.stderr)
            sys.exit(1)
        
        print(f"[{datetime.now().strftime('%H:%M:%S')}] Parsing complete: {self.entries_parsed} entries parsed, {self.entries_kept} entries kept")
        return entries


def trans_strands(strand: bool) -> str:
    """Convert boolean strand to character representation"""
    return '+' if strand else '-'


def parse_info_field(info_string):
    """
    Parse the INFO field and extract all key-value pairs.
    Special handling for PRECISION which appears as plain text.
    """
    info_dict = {}
    
    # Check for PRECISE or IMPRECISE at the beginning
    if info_string.startswith('PRECISE'):
        info_dict['PRECISION'] = 'PRECISE'
        info_string = info_string.replace('PRECISE;', '')
    elif info_string.startswith('IMPRECISE'):
        info_dict['PRECISION'] = 'IMPRECISE'
        info_string = info_string.replace('IMPRECISE;', '')
    
    # Parse the remaining key=value pairs
    if info_string:
        info_parts = info_string.split(';')
        for part in info_parts:
            if '=' in part:
                key, value = part.split('=', 1)
                info_dict[key] = value
    
    return info_dict


def parse_format_and_sample(format_string, sample_string):
    """
    Parse FORMAT and corresponding sample data.
    Returns a dictionary mapping FORMAT fields to sample values.
    """
    format_fields = format_string.split(':')
    sample_values = sample_string.split(':')
    
    # Create dictionary mapping format fields to values
    format_dict = {}
    for i, field in enumerate(format_fields):
        if i < len(sample_values):
            format_dict[field] = sample_values[i]
        else:
            format_dict[field] = ''
    
    return format_dict


def trans_type_str(sv_type: int) -> str:
    """Convert numeric SV type to string (for output)"""
    type_map = {
        SVType.DEL: "DEL",
        SVType.DUP: "DUP", 
        SVType.INV: "INV",
        SVType.TRA: "TRA",
        SVType.INS: "INS",
        SVType.BND: "BND",
        SVType.CNV: "CNV"
    }
    return type_map.get(sv_type, "NA")


def write_bed_output(entries: List[SVEntry], output_file: str, min_size: int, max_size: int, 
                     output_format: str, parser: VCFParser):
    """Write SV entries to BED format file"""
    written = 0
    filtered = 0
    
    # Define all possible INFO columns in order (from svParser.py)
    # Note: STRAND and SVTYPE removed as they're redundant with columns 9-10 and 11
    info_columns = ['PRECISION', 'AF', 'CHR2', 'COVERAGE', 'END', 'PHASE', 
                   'RNAMES', 'STDEV_LEN', 'STDEV_POS', 'SUPPORT', 
                   'SUPPORT_LONG', 'SVLEN']
    
    # Track which FORMAT fields we've seen
    format_columns = []
    format_columns_determined = False
    
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Writing {output_format} BED output to: {output_file}")
    
    try:
        with open(output_file, 'w') as f:
            # Write header based on format
            header_written = False
            
            for entry in entries:
                # Calculate size for filtering
                if entry.sv_type == SVType.INS:
                    size = entry.sv_len
                elif entry.sv_type in [SVType.TRA, SVType.BND]:
                    size = min_size + 1  # Always include translocations
                else:
                    size = abs(entry.stop.pos - entry.start.pos)
                
                # Apply size filters
                if size > min_size and (size < max_size or max_size == -1):
                    # Parse INFO and FORMAT fields for enhanced modes
                    if output_format != 'original':
                        info_data = parse_info_field(entry.info)
                        format_data = {}
                        if entry.format and entry.samples and len(entry.samples) > 0:
                            format_data = parse_format_and_sample(entry.format, entry.samples[0])
                        
                        # Determine FORMAT columns from first entry
                        if not format_columns_determined and format_data:
                            format_columns = list(format_data.keys())
                            format_columns_determined = True
                    
                    # Write header once we know the format
                    if not header_written:
                        if output_format == 'standard':
                            # Standard bedtools-compatible format
                            header_parts = [
                                "chrom", "start", "end", "name", "score", "strand", "svtype", 
                                "chr2", "end2", "cipos1", "cipos2", "ciend1", "ciend2", "strand2"
                            ]
                            if format_columns_determined:
                                # Add VCF fields
                                header_parts.extend(["REF", "ALT", "QUAL", "FILTER"])
                                # Add INFO columns
                                header_parts.extend(info_columns)
                                # Add FORMAT columns
                                header_parts.extend(format_columns)
                            f.write('#' + '\t'.join(header_parts) + '\n')
                        
                        elif output_format == 'enhanced':
                            # Enhanced format with all columns
                            header_parts = [
                                # Standard BED columns (1-6)
                                "chrom", "start", "end", "name", "score", "strand",
                                # Additional bedtools columns (7-14)
                                "svtype", "chr2", "end2", "cipos1", "cipos2", "ciend1", "ciend2", "strand2",
                                # Original SURVIVOR columns (15-29)
                                "chr1_survivor", "start1_cipos1", "start1_cipos2", "chr2_survivor", 
                                "end2_ciend1", "end2_ciend2", "sv_id", "comma", "strand1_survivor", 
                                "strand2_survivor", "svtype_survivor", "chr1_repeat", "start1_0based", 
                                "chr2_repeat", "end2_survivor"
                            ]
                            if format_columns_determined:
                                # VCF basic columns
                                header_parts.extend(["REF", "ALT", "QUAL", "FILTER"])
                                # Add INFO columns
                                header_parts.extend(info_columns)
                                # Add FORMAT columns
                                header_parts.extend(format_columns)
                            f.write('#' + '\t'.join(header_parts) + '\n')
                        
                        # Original format doesn't need a header
                        header_written = True
                    
                    # Build output based on format
                    if output_format == 'original':
                        # Original SURVIVOR format (15 columns)
                        bed_line = [
                            entry.start.chr,                                    # 1
                            str(entry.start.pos - 2 + entry.cpos[0]),          # 2
                            str(entry.start.pos - 2 + entry.cpos[1]),          # 3
                            entry.stop.chr,                                     # 4
                            str(entry.stop.pos - 1 + entry.cend[0]),          # 5
                            str(entry.stop.pos - 1 + entry.cend[1]),          # 6
                            entry.sv_id,                                        # 7
                            ',',                                                # 8
                            trans_strands(entry.strands.first),                # 9
                            trans_strands(entry.strands.second),               # 10
                            trans_type_str(entry.sv_type),                     # 11
                            entry.start.chr,                                    # 12
                            str(entry.start.pos - 1),                          # 13
                            entry.stop.chr,                                     # 14
                            str(entry.stop.pos)                                # 15
                        ]
                    
                    elif output_format == 'standard':
                        # Standard bedtools-compatible format
                        start_0based = str(entry.start.pos - 1)
                        # For end position: use end2 for same-chr, or start+1 for translocations
                        if entry.start.chr == entry.stop.chr and entry.sv_type != SVType.TRA:
                            end_pos = str(entry.stop.pos)
                        else:
                            end_pos = str(int(start_0based) + 1)
                        
                        bed_line = [
                            entry.start.chr,                                    # 1: chrom
                            start_0based,                                       # 2: start (0-based)
                            end_pos,                                           # 3: end (1-based)
                            entry.sv_id,                                        # 4: name
                            "1000",                                            # 5: score
                            trans_strands(entry.strands.first),                # 6: strand
                            trans_type_str(entry.sv_type),                     # 7: svtype
                            entry.stop.chr,                                     # 8: chr2
                            str(entry.stop.pos),                               # 9: end2
                            str(entry.start.pos - 2 + entry.cpos[0]),          # 10: cipos1
                            str(entry.start.pos - 2 + entry.cpos[1]),          # 11: cipos2
                            str(entry.stop.pos - 1 + entry.cend[0]),          # 12: ciend1
                            str(entry.stop.pos - 1 + entry.cend[1]),          # 13: ciend2
                            trans_strands(entry.strands.second)                # 14: strand2
                        ]
                        
                        # Add VCF fields if available
                        if format_columns_determined:
                            bed_line.extend([
                                entry.ref,
                                entry.alt,
                                entry.qual,
                                entry.filter
                            ])
                            # Add INFO columns in fixed order
                            for col in info_columns:
                                bed_line.append(info_data.get(col, ''))
                            # Add FORMAT columns in determined order
                            for col in format_columns:
                                bed_line.append(format_data.get(col, ''))
                    
                    elif output_format == 'enhanced':
                        # Enhanced format: bedtools + all original columns + VCF fields
                        start_0based = str(entry.start.pos - 1)
                        if entry.start.chr == entry.stop.chr and entry.sv_type != SVType.TRA:
                            end_pos = str(entry.stop.pos)
                        else:
                            end_pos = str(int(start_0based) + 1)
                        
                        bed_line = [
                            # Standard BED columns (1-6)
                            entry.start.chr,                                    # 1: chrom
                            start_0based,                                       # 2: start (0-based)
                            end_pos,                                           # 3: end (1-based)
                            entry.sv_id,                                        # 4: name
                            "1000",                                            # 5: score
                            trans_strands(entry.strands.first),                # 6: strand
                            # Additional bedtools columns (7-14)
                            trans_type_str(entry.sv_type),                     # 7: svtype
                            entry.stop.chr,                                     # 8: chr2
                            str(entry.stop.pos),                               # 9: end2
                            str(entry.start.pos - 2 + entry.cpos[0]),          # 10: cipos1
                            str(entry.start.pos - 2 + entry.cpos[1]),          # 11: cipos2
                            str(entry.stop.pos - 1 + entry.cend[0]),          # 12: ciend1
                            str(entry.stop.pos - 1 + entry.cend[1]),          # 13: ciend2
                            trans_strands(entry.strands.second),               # 14: strand2
                            # Original SURVIVOR columns (15-29)
                            entry.start.chr,                                    # 15: chr1_survivor
                            str(entry.start.pos - 2 + entry.cpos[0]),          # 16: start1_cipos1
                            str(entry.start.pos - 2 + entry.cpos[1]),          # 17: start1_cipos2
                            entry.stop.chr,                                     # 18: chr2_survivor
                            str(entry.stop.pos - 1 + entry.cend[0]),          # 19: end2_ciend1
                            str(entry.stop.pos - 1 + entry.cend[1]),          # 20: end2_ciend2
                            entry.sv_id,                                        # 21: sv_id
                            ',',                                                # 22: comma
                            trans_strands(entry.strands.first),                # 23: strand1_survivor
                            trans_strands(entry.strands.second),               # 24: strand2_survivor
                            trans_type_str(entry.sv_type),                     # 25: svtype_survivor
                            entry.start.chr,                                    # 26: chr1_repeat
                            str(entry.start.pos - 1),                          # 27: start1_0based
                            entry.stop.chr,                                     # 28: chr2_repeat
                            str(entry.stop.pos)                                # 29: end2_survivor
                        ]
                        
                        # Add VCF fields if available
                        if format_columns_determined:
                            bed_line.extend([
                                entry.ref,
                                entry.alt,
                                entry.qual,
                                entry.filter
                            ])
                            # Add INFO columns in fixed order
                            for col in info_columns:
                                bed_line.append(info_data.get(col, ''))
                            # Add FORMAT columns in determined order
                            for col in format_columns:
                                bed_line.append(format_data.get(col, ''))
                    
                    f.write('\t'.join(bed_line) + '\n')
                    written += 1
                else:
                    filtered += 1
    
    except Exception as e:
        print(f"ERROR: Failed to write output file: {str(e)}", file=sys.stderr)
        sys.exit(1)
    
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Output complete: {written} SVs written, {filtered} filtered by size")


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='Universal VCF to BED converter with three output formats',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Output formats:

1. Default (standard BED) - Bedtools-compatible format:
   - Columns 1-6: Standard BED (chrom, start, end, name, score, strand)
   - Columns 7-14: Additional SV info (svtype, chr2, end2, confidence intervals, strand2)
   - Columns 15+: VCF fields (REF, ALT, QUAL, FILTER, expanded INFO, expanded FORMAT)
   
2. Enhanced (--enhanced) - All information preserved:
   - Columns 1-6: Standard BED format
   - Columns 7-14: Additional bedtools columns
   - Columns 15-29: Original SURVIVOR columns
   - Columns 30+: VCF fields (REF, ALT, QUAL, FILTER, expanded INFO, expanded FORMAT)
   
3. Original (--original) - SURVIVOR vcftobed format:
   - 15 columns matching SURVIVOR vcftobed exactly
   - No headers, byte-for-byte identical output

Column descriptions:

Standard BED format (default):
  1. chrom - Chromosome
  2. start - Start position (0-based)
  3. end - End position (1-based, or start+1 for translocations)
  4. name - SV ID
  5. score - Score (1000 as placeholder)
  6. strand - Start strand
  7. svtype - SV type
  8. chr2 - End chromosome
  9. end2 - Original end position
  10-13. cipos1, cipos2, ciend1, ciend2 - Confidence intervals
  14. strand2 - End strand
  15+. VCF fields (when available)

Original SURVIVOR format (--original):
  1. Start chromosome
  2. Start position - 2 + CIPOS first value
  3. Start position - 2 + CIPOS second value
  4. End chromosome
  5. End position - 1 + CIEND first value
  6. End position - 1 + CIEND second value
  7. SV ID
  8. Comma (,)
  9. Start strand (+/-)
  10. End strand (+/-)
  11. SV type (DEL/DUP/INV/INS/TRA)
  12. Start chromosome (repeated)
  13. Start position - 1 (0-based)
  14. End chromosome (repeated)
  15. End position (1-based)

Examples:
  # Default: bedtools-compatible format
  %(prog)s input.vcf
  
  # Enhanced: all columns preserved
  %(prog)s input.vcf --enhanced
  
  # Original: SURVIVOR format
  %(prog)s input.vcf --original
  
  # With filtering
  %(prog)s input.vcf --min-size 50 --max-size 100000
  
  # Use with bedtools
  %(prog)s input.vcf -o svs.bed
  bedtools intersect -a svs.bed -b genes.bed
        """
    )
    
    # Positional arguments
    parser.add_argument('vcf_file', 
                       help='Input VCF file')
    
    # Optional arguments
    parser.add_argument('--min-size', '-m',
                       type=int,
                       default=0,
                       help='Minimum SV size to include (default: 0, no minimum)')
    
    parser.add_argument('--max-size', '-M', 
                       type=int,
                       default=-1,
                       help='Maximum SV size to include (default: -1, no maximum)')
    
    parser.add_argument('--output', '-o',
                       type=str,
                       help='Output BED file (default: input filename with .bed extension)')
    
    parser.add_argument('--original',
                       action='store_true',
                       help='Output original SURVIVOR vcftobed format (15 columns)')
    
    parser.add_argument('--enhanced',
                       action='store_true',
                       help='Output enhanced format with all columns preserved')
    
    parser.add_argument('--version', '-v',
                       action='version',
                       version='%(prog)s 3.0.0')
    
    parser.add_argument('--debug', '-d',
                       action='store_true',
                       help='Enable debug output')
    
    args = parser.parse_args()
    
    # Determine output format
    if args.original and args.enhanced:
        print("ERROR: Cannot specify both --original and --enhanced", file=sys.stderr)
        sys.exit(1)
    
    if args.original:
        output_format = 'original'
    elif args.enhanced:
        output_format = 'enhanced'
    else:
        output_format = 'standard'  # Default: bedtools-compatible
    
    # Validate input file
    if not os.path.exists(args.vcf_file):
        print(f"ERROR: Input file does not exist: {args.vcf_file}", file=sys.stderr)
        sys.exit(1)
    
    # Set output filename if not specified
    if not args.output:
        base_name = os.path.splitext(args.vcf_file)[0]
        args.output = base_name + '.bed'
    
    # Validate size parameters
    if args.min_size < -1:
        print("ERROR: min_size must be >= -1", file=sys.stderr)
        sys.exit(1)
        
    if args.max_size < -1:
        print("ERROR: max_size must be >= -1", file=sys.stderr)
        sys.exit(1)
        
    if args.max_size != -1 and args.min_size > args.max_size:
        print("ERROR: min_size cannot be greater than max_size", file=sys.stderr)
        sys.exit(1)
    
    # Print run parameters
    print(f"\n{'='*60}")
    print(f"snifftobed.py - Universal VCF to BED converter")
    print(f"{'='*60}")
    print(f"Input VCF: {args.vcf_file}")
    print(f"Output BED: {args.output}")
    print(f"Output format: {output_format}")
    print(f"Min SV size: {args.min_size if args.min_size >= 0 else 'No minimum'}")
    print(f"Max SV size: {args.max_size if args.max_size >= 0 else 'No maximum'}")
    print(f"{'='*60}\n")
    
    # Parse VCF file
    parser_obj = VCFParser(min_size=args.min_size)
    entries = parser_obj.parse_vcf(args.vcf_file, args.min_size)
    
    if not entries:
        print("WARNING: No structural variants found in input file", file=sys.stderr)
    
    # Write BED output
    write_bed_output(entries, args.output, args.min_size, args.max_size, output_format, parser_obj)
    
    print(f"\n[{datetime.now().strftime('%H:%M:%S')}] Conversion complete!")
    print(f"{'='*60}\n")


if __name__ == '__main__':
    main()
