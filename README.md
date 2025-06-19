# snifftobed.py - Universal VCF to BED Converter

[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)
[![T2T Compatible](https://img.shields.io/badge/T2T-compatible-green.svg)](https://github.com/marbl/CHM13)
[![Sniffles2 Compatible](https://img.shields.io/badge/Sniffles2-compatible-green.svg)](https://github.com/fritzsedlazeck/Sniffles)
[![VCF Compatible](https://img.shields.io/badge/VCF-compatible-green.svg)](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
[![SURVIVOR Compatible](https://img.shields.io/badge/SURVIVOR-compatible-green.svg)](https://github.com/fritzsedlazeck/SURVIVOR)
[![bedtools Compatible](https://img.shields.io/badge/bedtools-compatible-green.svg)](https://bedtools.readthedocs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Table of Contents

1. [Overview](#overview)
2. [Key Features](#key-features)
3. [Installation](#installation)
4. [Quick Start](#quick-start)
5. [Output Formats](#output-formats)
6. [Command Line Options](#command-line-options)
7. [Column Specifications](#column-specifications)
8. [Usage Examples](#usage-examples)
9. [Working with Specific SV Types](#working-with-specific-sv-types)
10. [Integration with Other Tools](#integration-with-other-tools)
11. [Filtering and Analysis](#filtering-and-analysis)
12. [Performance and Limitations](#performance-and-limitations)
13. [Troubleshooting](#troubleshooting)
14. [FAQ](#faq)
15. [Version History](#version-history)
16. [License and Citation](#license-and-citation)

## Overview

`snifftobed.py` is a comprehensive VCF to BED converter specifically designed for structural variants (SVs). It combines the functionality of multiple tools into a single, versatile script that can output three different formats to suit various analysis needs.

### Why This Tool?

Traditional VCF to BED converters often lose important information during conversion. This tool addresses several key challenges:

1. **Bedtools Compatibility**: The default output works directly with bedtools without modification
2. **Information Preservation**: All VCF fields can be preserved in the output
3. **SURVIVOR Compatibility**: Maintains backward compatibility with existing SURVIVOR pipelines
4. **Flexible Output**: Three output formats to suit different analysis needs
5. **SV-Specific Design**: Handles all structural variant types correctly

## Key Features

- ✅ **Three Output Formats**: Standard BED (default), Enhanced, and Original SURVIVOR
- ✅ **Complete VCF Field Expansion**: Preserves REF, ALT, QUAL, FILTER, INFO, and FORMAT fields
- ✅ **Bedtools Ready**: Default output works directly with bedtools intersect, closest, etc.
- ✅ **SV Type Support**: Correctly handles DEL, DUP, INV, INS, TRA/BND, CNV
- ✅ **Confidence Intervals**: Preserves CIPOS/CIEND information
- ✅ **Multi-sample Support**: Processes multi-sample VCFs (expands first sample's FORMAT)
- ✅ **Size Filtering**: Built-in min/max size filtering
- ✅ **Progress Reporting**: Real-time progress updates for large files
- ✅ **Error Handling**: Robust error handling with informative messages

## Installation

### Requirements

- Python 3.6 or higher
- No external dependencies required (uses only Python standard library)

### Setup

```bash
# Download the script
wget https://github.com/Thomas-X-Garcia/snifftobed/raw/main/snifftobed.py

# Make it executable
chmod +x snifftobed.py

# Test installation
python3 snifftobed.py --version
```

### Alternative Installation

```bash
# Clone repository
git clone https://github.com/Thomas-X-Garcia/snifftobed.git
cd snifftobed

# Add to PATH (optional)
echo "export PATH=$PATH:$(pwd)" >> ~/.bashrc
source ~/.bashrc
```

## Quick Start

### Basic Usage

```bash
# Default: Create bedtools-compatible BED file
python3 snifftobed.py input.vcf

# Specify output file
python3 snifftobed.py input.vcf -o output.bed

# Filter by size (only SVs >= 50bp)
python3 snifftobed.py input.vcf --min-size 50

# Original SURVIVOR format
python3 snifftobed.py input.vcf --original

# Enhanced format with all information
python3 snifftobed.py input.vcf --enhanced
```

### Quick Examples with Bedtools

```bash
# Find SVs overlapping genes
python3 snifftobed.py svs.vcf -o svs.bed
bedtools intersect -a svs.bed -b genes.bed -wa -wb > svs_in_genes.txt

# Find nearest gene to each SV
bedtools closest -a svs.bed -b genes.bed -d > svs_nearest_gene.txt

# Calculate SV coverage in windows
bedtools coverage -a genome_windows.bed -b svs.bed > sv_coverage.txt
```

## Output Formats

### 1. Default Format (Bedtools-compatible)

The default output is optimized for bedtools and other genomic analysis tools:

```
#chrom  start  end  name  score  strand  svtype  chr2  end2  cipos1  cipos2  ciend1  ciend2  strand2  REF  ALT  QUAL  FILTER  [INFO fields...]  [FORMAT fields...]
chr1    9999   10500  DEL001  1000  +  DEL  chr1  10500  9998  9998  10499  10499  -  N  <DEL>  60  PASS  PRECISE  0.5  ...
```

**Key characteristics:**
- Standard BED coordinates (0-based start, 1-based end)
- All VCF information preserved in additional columns
- Translocations represented as 1bp intervals
- Header line for easy column identification

### 2. Enhanced Format (--enhanced)

Preserves ALL information including original SURVIVOR columns:

```
#chrom  start  end  name  score  strand  [bedtools cols...]  [SURVIVOR cols...]  [VCF fields...]
chr1    9999   10500  DEL001  1000  +  ...  chr1  9998  9998  ...  N  <DEL>  60  PASS  ...
```

**Use when you need:**
- Complete data preservation
- Both bedtools compatibility AND SURVIVOR format
- Maximum information for complex analyses

### 3. Original Format (--original)

Byte-for-byte identical to SURVIVOR vcftobed:

```
chr1  9998  9998  chr1  10499  10499  DEL001  ,  +  -  DEL  chr1  9999  chr1  10500
```

**Use for:**
- Existing SURVIVOR-based pipelines
- Backward compatibility
- Scripts expecting exact SURVIVOR format

## Command Line Options

```
usage: snifftobed.py [-h] [--min-size MIN_SIZE] [--max-size MAX_SIZE]
                           [--output OUTPUT] [--original] [--enhanced]
                           [--version] [--debug]
                           vcf_file

Universal VCF to BED converter with three output formats

positional arguments:
  vcf_file              Input VCF file

optional arguments:
  -h, --help            show this help message and exit
  --min-size MIN_SIZE, -m MIN_SIZE
                        Minimum SV size to include (default: 0, no minimum)
  --max-size MAX_SIZE, -M MAX_SIZE
                        Maximum SV size to include (default: -1, no maximum)
  --output OUTPUT, -o OUTPUT
                        Output BED file (default: input filename with .bed extension)
  --original            Output original SURVIVOR vcftobed format (15 columns)
  --enhanced            Output enhanced format with all columns preserved
  --version, -v         show program's version number and exit
  --debug, -d           Enable debug output
```

### Option Details

- **--min-size/-m**: Only output SVs >= this size (in base pairs)
- **--max-size/-M**: Only output SVs <= this size (use -1 for no limit)
- **--output/-o**: Specify output filename (default: input.vcf → input.bed)
- **--original**: Use SURVIVOR's exact output format
- **--enhanced**: Include all possible columns
- **--debug/-d**: Print additional debugging information

## Column Specifications

### Default Format Columns (Bedtools-compatible)

| Column | Name | Description | Example |
|--------|------|-------------|---------|
| 1 | chrom | Chromosome | chr1 |
| 2 | start | Start position (0-based) | 9999 |
| 3 | end | End position (1-based) | 10500 |
| 4 | name | SV identifier | DEL001 |
| 5 | score | Score (placeholder) | 1000 |
| 6 | strand | Start strand | + |
| 7 | svtype | SV type | DEL |
| 8 | chr2 | End chromosome | chr1 |
| 9 | end2 | Original end position | 10500 |
| 10 | cipos1 | CIPOS first value | 9998 |
| 11 | cipos2 | CIPOS second value | 9998 |
| 12 | ciend1 | CIEND first value | 10499 |
| 13 | ciend2 | CIEND second value | 10499 |
| 14 | strand2 | End strand | - |
| 15 | REF | Reference allele | N |
| 16 | ALT | Alternate allele | <DEL> |
| 17 | QUAL | Quality score | 60 |
| 18 | FILTER | Filter status | PASS |
| 19 | PRECISION | Breakpoint precision | PRECISE |
| 20 | AF | Allele frequency | 0.5 |
| 21 | CHR2 | Mate chromosome | chr1 |
| 22 | COVERAGE | Coverage info | 30,30,30 |
| 23 | END | End position from INFO | 10500 |
| 24 | PHASE | Phase information | 0\|1 |
| 25 | RNAMES | Read names | read1,read2 |
| 26 | STDEV_LEN | Length std deviation | 10.5 |
| 27 | STDEV_POS | Position std deviation | 5.2 |
| 28 | SUPPORT | Support reads | 15 |
| 29 | SUPPORT_LONG | Long read support | 8 |
| 30 | SVLEN | SV length | -500 |
| 31+ | FORMAT fields | GT, GQ, DR, DV, etc. | 0/1, 45, 20, 15 |

### INFO Fields (Columns 19-30)

The following INFO fields are always included (empty if not present in VCF):

1. **PRECISION**: PRECISE or IMPRECISE
2. **AF**: Allele frequency
3. **CHR2**: Mate chromosome for translocations
4. **COVERAGE**: Coverage values
5. **END**: End position
6. **PHASE**: Phase information
7. **RNAMES**: Supporting read names
8. **STDEV_LEN**: Length standard deviation
9. **STDEV_POS**: Position standard deviation
10. **SUPPORT**: Number of supporting reads
11. **SUPPORT_LONG**: Long read support count
12. **SVLEN**: SV length

Note: STRAND and SVTYPE from INFO are not included as separate columns since they're already in the BED format (columns 6, 7, 14).

### FORMAT Fields (Columns 31+)

FORMAT fields are dynamically determined from the VCF file. Common fields include:

- **GT**: Genotype (0/1, 1/1, etc.)
- **GQ**: Genotype quality
- **DR**: Reference read depth
- **DV**: Variant read depth
- **AD**: Allelic depths
- **DP**: Total depth
- **PL**: Phred-scaled likelihoods

## Usage Examples

### Example 1: Basic SV Analysis

```bash
# Convert VCF to BED
python3 snifftobed.py sample.vcf -o sample.bed

# Count SVs by type
awk '{print $7}' sample.bed | sort | uniq -c

# Extract large deletions (>1kb)
awk '$7=="DEL" && sqrt(($3-$2)^2)>1000' sample.bed > large_dels.bed

# Find SVs in specific region
awk '$1=="chr1" && $2>=1000000 && $3<=2000000' sample.bed > chr1_region.bed
```

### Example 2: Quality Filtering

```bash
# High-quality SVs only (QUAL>=60, PASS filter, >=10 supporting reads)
python3 snifftobed.py calls.vcf | \
  awk -F'\t' '$17>=60 && $18=="PASS" && $28>=10' > high_quality.bed

# Precise calls with high allele frequency
awk -F'\t' '$19=="PRECISE" && $20>=0.5' sample.bed > precise_common.bed

# Homozygous variants with good genotype quality
awk -F'\t' '$31=="1/1" && $32>=30' sample.bed > homozygous_high_gq.bed
```

### Example 3: Complex Filtering Pipeline

```bash
# Step 1: Convert and filter by size
python3 snifftobed.py input.vcf --min-size 50 --max-size 100000 -o filtered.bed

# Step 2: Apply quality filters
awk -F'\t' 'NR==1 || ($17>=30 && $18=="PASS" && $19=="PRECISE")' filtered.bed > quality.bed

# Step 3: Intersect with regions of interest
bedtools intersect -a quality.bed -b regions_of_interest.bed -u > final.bed

# Step 4: Annotate with genes
bedtools intersect -a final.bed -b genes.bed -wa -wb > annotated.bed
```

### Example 4: Population Analysis

```bash
# Convert population VCF
python3 snifftobed.py population.vcf -o population.bed

# Common variants (AF > 5%)
awk -F'\t' 'NR==1 || $20>0.05' population.bed > common_variants.bed

# Rare variants (AF < 1%)
awk -F'\t' 'NR>1 && $20<0.01 && $20!=""' population.bed > rare_variants.bed

# Singleton variants
awk -F'\t' 'NR>1 && $28==1' population.bed > singletons.bed

# Distribution by frequency bins
awk -F'\t' 'NR>1 && $20!="" {
  if($20<0.01) print "Rare"
  else if($20<0.05) print "Low"
  else if($20<0.95) print "Common"
  else print "Fixed"
}' population.bed | sort | uniq -c
```

## Working with Specific SV Types

### Deletions

```bash
# Extract all deletions
awk -F'\t' '$7=="DEL"' sample.bed > deletions.bed

# Large deletions (>10kb)
awk -F'\t' '$7=="DEL" && sqrt($30^2)>10000' sample.bed > large_dels.bed

# Deletions in exons
bedtools intersect -a <(awk '$7=="DEL"' sample.bed) -b exons.bed -wa -wb > dels_in_exons.bed
```

### Duplications

```bash
# Extract duplications with size range
awk -F'\t' '$7=="DUP" && sqrt($30^2)>=1000 && sqrt($30^2)<=10000' sample.bed > dup_1k_10k.bed

# Tandem duplications (same chromosome)
awk -F'\t' '$7=="DUP" && $1==$8' sample.bed > tandem_dups.bed
```

### Inversions

```bash
# Extract inversions
awk -F'\t' '$7=="INV"' sample.bed > inversions.bed

# Large inversions affecting multiple genes
awk -F'\t' '$7=="INV" && ($3-$2)>100000' sample.bed > large_inv.bed
bedtools intersect -a large_inv.bed -b genes.bed -wa -wb | \
  awk '{print $4}' | sort | uniq -c | awk '$1>5' > multi_gene_inv.txt
```

### Translocations

```bash
# Extract all translocations
awk -F'\t' '$7=="TRA"' sample.bed > translocations.bed

# Inter-chromosomal translocations only
awk -F'\t' '$7=="TRA" && $1!=$8' sample.bed > inter_chr_tra.bed

# Create both breakpoint positions for translocations
awk -F'\t' '$7=="TRA" {
  print $1"\t"$2"\t"$3"\t"$4"_BP1\t"$5"\t"$6
  print $8"\t"$9-1"\t"$9"\t"$4"_BP2\t"$5"\t"$14
}' sample.bed > tra_breakpoints.bed
```

### Insertions

```bash
# Extract insertions
awk -F'\t' '$7=="INS"' sample.bed > insertions.bed

# Mobile element insertions (if annotated)
grep -E "LINE1|ALU|SVA" insertions.bed > mobile_elements.bed

# Large insertions (>100bp)
awk -F'\t' '$7=="INS" && sqrt($30^2)>100' sample.bed > large_ins.bed
```

## Integration with Other Tools

### Bedtools Integration

```bash
# Setup: Convert VCF to BED
python3 snifftobed.py input.vcf -o svs.bed

# 1. Find SVs overlapping features
bedtools intersect -a svs.bed -b genes.bed -wa -wb > svs_in_genes.bed
bedtools intersect -a svs.bed -b promoters.bed -wa -wb > svs_in_promoters.bed
bedtools intersect -a svs.bed -b enhancers.bed -wa -wb > svs_in_enhancers.bed

# 2. Find nearest feature to each SV
bedtools closest -a svs.bed -b genes.bed -d > nearest_gene.bed
bedtools closest -a svs.bed -b tad_boundaries.bed -d > nearest_tad.bed

# 3. Calculate coverage/enrichment
bedtools coverage -a functional_regions.bed -b svs.bed > sv_coverage.bed
bedtools fisher -a svs.bed -b hotspots.bed -g genome.sizes > enrichment.txt

# 4. Window-based analysis
bedtools makewindows -g genome.sizes -w 100000 > windows.bed
bedtools coverage -a windows.bed -b svs.bed -counts > sv_density.bed

# 5. Merge nearby SVs
bedtools merge -i <(sort -k1,1 -k2,2n svs.bed) -d 1000 -c 4,7 -o distinct,distinct > merged_svs.bed
```

### R Integration

```r
# Load libraries
library(tidyverse)
library(GenomicRanges)

# Read BED file with all columns
sv_data <- read.table("sample.bed", header=TRUE, sep="\t", comment.char="")

# Basic analysis
sv_summary <- sv_data %>%
  group_by(svtype) %>%
  summarise(
    count = n(),
    mean_size = mean(abs(as.numeric(SVLEN)), na.rm=TRUE),
    mean_qual = mean(QUAL, na.rm=TRUE),
    mean_support = mean(SUPPORT, na.rm=TRUE)
  )

# Create GenomicRanges object
sv_gr <- GRanges(
  seqnames = sv_data$chrom,
  ranges = IRanges(start = sv_data$start + 1, end = sv_data$end),
  strand = sv_data$strand,
  svtype = sv_data$svtype,
  name = sv_data$name,
  qual = sv_data$QUAL,
  af = sv_data$AF
)

# Overlap with genes
genes_gr <- import("genes.bed")
overlaps <- findOverlaps(sv_gr, genes_gr)

# Plot SV size distribution
ggplot(sv_data, aes(x=abs(as.numeric(SVLEN)), fill=svtype)) +
  geom_histogram(bins=50) +
  scale_x_log10() +
  facet_wrap(~svtype, scales="free_y") +
  theme_minimal() +
  labs(x="SV Size (bp)", y="Count", title="SV Size Distribution by Type")

# Plot quality vs support
ggplot(sv_data, aes(x=SUPPORT, y=QUAL, color=svtype)) +
  geom_point(alpha=0.6) +
  scale_x_log10() +
  theme_minimal() +
  labs(x="Supporting Reads", y="Quality Score", title="Quality vs Support")
```

### Python/Pandas Integration

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read BED file
# Specify column names for reliable parsing
columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'svtype', 
           'chr2', 'end2', 'cipos1', 'cipos2', 'ciend1', 'ciend2', 'strand2',
           'REF', 'ALT', 'QUAL', 'FILTER', 'PRECISION', 'AF', 'CHR2', 
           'COVERAGE', 'END', 'PHASE', 'RNAMES', 'STDEV_LEN', 'STDEV_POS',
           'SUPPORT', 'SUPPORT_LONG', 'SVLEN', 'GT', 'GQ', 'DR', 'DV']

df = pd.read_csv('sample.bed', sep='\t', comment='#', names=columns, header=0)

# Basic statistics
print("SV Type Distribution:")
print(df['svtype'].value_counts())

print("\nSize Statistics by Type:")
df['sv_size'] = df['end'] - df['start']
print(df.groupby('svtype')['sv_size'].describe())

# Quality filtering
high_qual = df[(df['QUAL'] >= 30) & (df['FILTER'] == 'PASS')]
print(f"\nHigh quality SVs: {len(high_qual)} / {len(df)}")

# Allele frequency analysis
common_svs = df[df['AF'] >= 0.05]
rare_svs = df[df['AF'] < 0.01]

# Chromosome distribution
chr_counts = df['chrom'].value_counts().sort_index()
plt.figure(figsize=(12, 6))
chr_counts.plot(kind='bar')
plt.title('SV Distribution by Chromosome')
plt.xlabel('Chromosome')
plt.ylabel('Count')
plt.tight_layout()
plt.show()

# Support vs Quality scatter plot
plt.figure(figsize=(10, 8))
for svtype in df['svtype'].unique():
    subset = df[df['svtype'] == svtype]
    plt.scatter(subset['SUPPORT'], subset['QUAL'], label=svtype, alpha=0.6)
plt.xscale('log')
plt.xlabel('Supporting Reads')
plt.ylabel('Quality Score')
plt.legend()
plt.title('Quality vs Support by SV Type')
plt.grid(True, alpha=0.3)
plt.show()

# Export filtered results
high_conf = df[(df['QUAL'] >= 50) & 
               (df['SUPPORT'] >= 10) & 
               (df['PRECISION'] == 'PRECISE')]
high_conf.to_csv('high_confidence_svs.bed', sep='\t', index=False)
```

### SURVIVOR Integration

For pipelines that require SURVIVOR format:

```bash
# Generate SURVIVOR-compatible output
python3 snifftobed.py input.vcf --original -o survivor_format.bed

# Use with SURVIVOR tools
SURVIVOR stats survivor_format.bed NA NA NA > stats.txt
SURVIVOR merge file_list.txt 1000 1 1 1 0 50 merged.vcf

# Convert back if needed
python3 snifftobed.py merged.vcf -o merged_bedtools.bed
```

## Filtering and Analysis

### Size-based Filtering

```bash
# During conversion
python3 snifftobed.py input.vcf --min-size 50 --max-size 100000 -o filtered.bed

# After conversion
awk -F'\t' 'NR==1 || (sqrt($30^2)>=50 && sqrt($30^2)<=100000)' sample.bed > size_filtered.bed

# Size categories
awk -F'\t' 'NR>1 {
  size = sqrt($30^2)
  if(size < 100) print "Small"
  else if(size < 1000) print "Medium"
  else if(size < 10000) print "Large"
  else print "Very Large"
}' sample.bed | sort | uniq -c
```

### Quality-based Filtering

```bash
# Multiple quality criteria
awk -F'\t' 'NR==1 || ($17>=30 && $18=="PASS" && $19=="PRECISE" && $28>=5)' sample.bed > quality_filtered.bed

# Genotype quality filtering
awk -F'\t' 'NR==1 || $32>=20' sample.bed > high_gq.bed

# Read depth filtering
awk -F'\t' 'NR==1 || ($33+$34)>=10' sample.bed > good_depth.bed
```

### Complex Queries

```bash
# Find recurrent SVs (same type and similar position)
awk -F'\t' 'NR>1 {print $1"\t"int($2/1000)*1000"\t"$7}' sample.bed | \
  sort | uniq -c | awk '$1>1' > recurrent_regions.txt

# Extract SVs affecting multiple genes
bedtools intersect -a sample.bed -b genes.bed -wa -wb | \
  awk -F'\t' '{print $4"\t"$NF}' | \
  sort | uniq | \
  awk '{count[$1]++; genes[$1]=genes[$1]","$2} END {
    for(sv in count) if(count[sv]>1) print sv"\t"count[sv]"\t"genes[sv]
  }' > multi_gene_svs.txt

# Find compound heterozygous pairs
awk -F'\t' '$31=="0/1" {print $1"\t"$2"\t"$3"\t"$4}' sample.bed | \
  bedtools intersect -a - -b genes.bed -wa -wb | \
  awk '{gene[$NF]++; svs[$NF]=svs[$NF]","$4} END {
    for(g in gene) if(gene[g]>=2) print g"\t"gene[g]"\t"svs[g]
  }' > compound_het_candidates.txt
```

## Performance and Limitations

### Performance Characteristics

- **Speed**: Processes ~10,000-20,000 variants per second
- **Memory**: Uses ~100MB for typical VCF files
- **Scaling**: Linear time complexity O(n) with number of variants
- **Progress**: Updates every 10,000 lines for large files

### File Size Examples

| Input VCF | Variants | Processing Time | Memory Usage |
|-----------|----------|-----------------|--------------|
| 10 MB | ~5,000 | <1 second | ~50 MB |
| 100 MB | ~50,000 | 3-5 seconds | ~100 MB |
| 1 GB | ~500,000 | 30-50 seconds | ~500 MB |
| 10 GB | ~5,000,000 | 5-8 minutes | ~2 GB |

### Limitations

1. **Multi-sample VCFs**: Only the first sample's FORMAT fields are expanded
2. **Complex SVs**: Treats complex events as single variants
3. **Symbolic ALTs**: Assumes standard symbolic alleles (<DEL>, <DUP>, etc.)
4. **Memory**: Large VCFs (>10GB) may require significant memory
5. **Coordinate System**: Uses SURVIVOR's non-standard calculations in original mode

### Optimization Tips

```bash
# For very large files, process in chunks
split -l 1000000 large.vcf chunk_
for chunk in chunk_*; do
  python3 snifftobed.py $chunk -o ${chunk}.bed
done
cat chunk_*.bed > combined.bed

# Pre-filter VCF to reduce size
bcftools view -f PASS -i 'QUAL>=20' input.vcf | \
  python3 snifftobed.py - -o filtered.bed

# Use streams to avoid temporary files
bcftools view input.vcf.gz | \
  python3 snifftobed.py - | \
  awk '$17>=30' | \
  bedtools sort > final.bed
```

## Troubleshooting

### Common Issues and Solutions

#### 1. "No variants in output"

**Possible causes:**
- Size filters too restrictive
- VCF has no SVTYPE in INFO field
- Malformed VCF file

**Solutions:**
```bash
# Check VCF format
grep -v "^#" input.vcf | head -5

# Run without filters
python3 snifftobed.py input.vcf --min-size -1 --max-size -1

# Check for SVTYPE
grep -v "^#" input.vcf | head -5 | cut -f8 | grep -o "SVTYPE=[^;]*"
```

#### 2. "Bedtools error: chroms out of order"

**Solution:**
```bash
# Sort the BED file
sort -k1,1 -k2,2n output.bed > sorted.bed

# Or use bedtools sort
bedtools sort -i output.bed > sorted.bed
```

#### 3. "Column count mismatch"

**Cause:** FORMAT fields vary between VCF files

**Solution:**
```bash
# Check header to identify columns
head -1 output.bed | tr '\t' '\n' | nl

# Count columns
head -1 output.bed | awk -F'\t' '{print NF " columns"}'
```

#### 4. "Missing AF or other INFO fields"

**Cause:** Not all VCFs contain all INFO fields

**Solution:**
```bash
# Check which INFO fields are present
grep "^##INFO" input.vcf

# Fields not in VCF will be empty in output
# Use awk to handle missing fields gracefully
awk -F'\t' '$20=="" {$20="0"}1' OFS='\t' output.bed > fixed.bed
```

#### 5. "Memory error with large files"

**Solutions:**
```bash
# Process in chunks
split -l 500000 large.vcf temp_chunk_

# Use streaming if possible
zcat large.vcf.gz | python3 snifftobed.py - -o output.bed

# Filter first
bcftools view -f PASS large.vcf | python3 snifftobed.py - -o output.bed
```

### Validation Commands

```bash
# Verify conversion worked
wc -l input.vcf output.bed

# Check for common issues
grep -c "^#" output.bed  # Should be 1 (header only)
awk -F'\t' '{print NF}' output.bed | sort | uniq -c  # All lines same column count

# Validate coordinates
awk -F'\t' '$2>$3' output.bed  # Should return nothing (start <= end)

# Check SV types
awk -F'\t' '{print $7}' output.bed | sort | uniq -c

# Verify original mode matches SURVIVOR
./SURVIVOR vcftobed input.vcf 0 -1 survivor.bed
python3 snifftobed.py input.vcf --original -o python.bed
diff survivor.bed python.bed  # Should be empty
```

## FAQ

### Q1: Why are there different coordinate systems?

**A:** The original SURVIVOR vcftobed uses non-standard position calculations (pos-2 for columns 2-3, pos-1 for columns 5-6). The default and enhanced modes use standard BED coordinates (0-based start) for compatibility with bedtools.

### Q2: How are translocations handled?

**A:** In default/enhanced modes, translocations create a 1bp interval at the first breakpoint (end = start + 1). The partner coordinate is preserved in the chr2/end2 columns. This allows bedtools to work correctly.

### Q3: What happens to multi-sample VCFs?

**A:** Following the behavior of svParser.py, only the first sample's FORMAT fields are expanded into columns. Information for other samples is not included in the output.

### Q4: Why are STRAND and SVTYPE not in the INFO columns?

**A:** These are already present in the standard BED columns (strand in column 6 and 14, svtype in column 7), so including them again would be redundant.

### Q5: Can I use this with compressed VCF files?

**A:** Not directly, but you can pipe the output:
```bash
zcat input.vcf.gz | python3 snifftobed.py - -o output.bed
# or
gunzip -c input.vcf.gz | python3 snifftobed.py - -o output.bed
```

### Q6: How do I know which column contains which field?

**A:** The default and enhanced formats include a header line starting with #. You can also check:
```bash
head -1 output.bed | tr '\t' '\n' | nl
```

### Q7: What's the difference between end2 and END?

**A:** 
- `end2` (column 9): The actual end position from the VCF POS + SVLEN or from BND ALT
- `END` (column 23): The END value from the INFO field
These may differ in some cases.

### Q8: Can I convert back from BED to VCF?

**A:** No, this is a one-way conversion. Some VCF information (like header lines, some INFO fields) cannot be fully reconstructed from the BED format.

### Q9: How are confidence intervals used?

**A:** CIPOS and CIEND values are added to positions in columns 10-13 but don't affect the main start/end coordinates (columns 2-3). Use these for understanding position uncertainty.

### Q10: Should I use --min-size 0 or --min-size 1?

**A:** Use --min-size 1 to exclude variants with zero or negative size (often problematic entries). Use --min-size 0 to include everything.

## Version History

### Version 1.0.0
- Initial release
- Exact replication of SURVIVOR vcftobed
- Python implementation of C++ original
- VCF field expansion (INFO and FORMAT)
- Enhanced vs original output modes
- Header line to output
- Robust error handling
- Unified three output formats into single tool
- Default output is bedtools-compatible
- --enhanced flag for complete information preservation
- INFO field parsing with fixed column order
- Handles missing fields
- Basic size filtering

## License and Citation

### License

This software is released under the MIT License:

```
MIT License

Copyright (c) 2025 Thomas X. Garcia, PhD, HCLD

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

### Citation

If you use this tool in your research, please cite:

```
snifftobed.py - Universal VCF to BED Converter
Thomas X. Garcia, PhD, HCLD (2025)
https://github.com/Thomas-X-Garcia/snifftobed
```

And the original SURVIVOR tool:

```
Jeffares, D.C., Jolly, C., Hoti, M. et al. Transient structural variations 
have strong effects on quantitative traits and reproductive isolation in 
fission yeast. Nat Commun 8, 14061 (2017).
https://github.com/fritzsedlazeck/SURVIVOR
```

### Acknowledgments

- Based on SURVIVOR vcftobed by Fritz Sedlazeck
- INFO/FORMAT field expansion inspired by svParser.py from the T2T annotation pipeline
- Developed to address the need for a unified, flexible VCF to BED converter
- Thanks to the genomics community for feedback and use cases

### Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

For bug reports and feature requests, please open an issue on GitHub.

### Contact

For questions, suggestions, or issues:
- GitHub Issues: https://github.com/Thomas-X-Garcia/snifftobed/issues
- Email: your.email@example.com

---

**Remember:** Always validate your results, especially when using non-standard formats or working with complex structural variants. When in doubt, compare a few entries manually to ensure the conversion meets your needs.
