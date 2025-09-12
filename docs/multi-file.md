# Multi-File Processing Tutorial

Learn how to process multiple VCF files simultaneously with CheckRef for efficient whole-genome or multi-chromosome analysis.

## Learning Objectives

By the end of this tutorial, you'll be able to:
- Process multiple chromosomes simultaneously
- Use different input methods (comma-separated, glob patterns)
- Understand automatic chromosome matching
- Interpret combined results across chromosomes

## Prerequisites

- Completed the [Basic Usage](./basic-usage) tutorial
- Multiple VCF files (typically per-chromosome)
- Reference panel with matching chromosome files
- Understanding of CheckRef output structure

## Multi-File Input Methods

CheckRef supports several ways to specify multiple input files:

### Method 1: Comma-Separated List

Process specific files by listing them with commas:

```bash
nextflow run main.nf \
  --targetVcfs "chr20.vcf.gz,chr21.vcf.gz,chr22.vcf.gz" \
  --referenceDir "/path/to/reference/panels/" \
  --outputDir multi_results \
  --fixMethod correct \
  -profile singularity
```

**When to use:**
- Processing specific chromosomes
- Files in different directories
- Non-standard naming patterns

### Method 2: Glob Patterns

Use wildcards to match multiple files:

```bash
nextflow run main.nf \
  --targetVcfs "/data/vcfs/*.vcf.gz" \
  --referenceDir "/path/to/reference/panels/" \
  --outputDir multi_results \
  --fixMethod correct \
  -profile singularity
```

**Pattern examples:**
```bash
# All VCF files in directory
--targetVcfs "/data/vcfs/*.vcf.gz"

# Chromosome-specific pattern
--targetVcfs "/data/vcfs/sample_chr*.vcf.gz"  

# Specific chromosome range
--targetVcfs "/data/vcfs/chr2[0-2].vcf.gz"
```

**When to use:**
- Standard file naming conventions
- Processing entire directories
- Automated workflows

### Method 3: Mixed Paths

Combine files from different locations:

```bash
nextflow run main.nf \
  --targetVcfs "/path1/chr20.vcf.gz,/path2/chr21.vcf.gz,/path3/chr22.vcf.gz" \
  --referenceDir "/ref/panels/" \
  --outputDir mixed_results \
  -profile singularity
```

## Chromosome Detection and Matching

CheckRef automatically detects chromosome information from filenames:

### Supported Naming Patterns
```
chr22.vcf.gz          → chromosome 22
sample_chr22.vcf.gz   → chromosome 22
chr22_filtered.vcf.gz → chromosome 22
chromosome22.vcf.gz   → chromosome 22
22.vcf.gz             → chromosome 22
```

### Reference Panel Matching
CheckRef automatically matches VCF files with reference legend files:

**Example file structure:**
```
Target files:
├── sample_chr20.vcf.gz
├── sample_chr21.vcf.gz
└── sample_chr22.vcf.gz

Reference panels:
├── 1000GP_chr20.legend.gz
├── 1000GP_chr21.legend.gz  
└── 1000GP_chr22.legend.gz
```

CheckRef will automatically pair:
- `sample_chr20.vcf.gz` ↔ `1000GP_chr20.legend.gz`
- `sample_chr21.vcf.gz` ↔ `1000GP_chr21.legend.gz`
- `sample_chr22.vcf.gz` ↔ `1000GP_chr22.legend.gz`

## Complete Multi-File Example

Here's a comprehensive example processing three chromosomes:

```bash
nextflow run main.nf \
  --targetVcfs "sample_chr20.vcf.gz,sample_chr21.vcf.gz,sample_chr22.vcf.gz" \
  --referenceDir "/cbio/dbs/refpanels/h3a_reference_panels/version_7/v7hc_s/sites/" \
  --legendPattern "*chr*.legend.gz" \
  --outputDir multi_chromosome_results \
  --fixMethod correct \
  --max_cpus 8 \
  --max_memory '16.GB' \
  -profile singularity \
  -resume
```

### Command Breakdown
- **targetVcfs**: Three specific chromosome files
- **legendPattern**: Pattern to match reference files
- **max_cpus**: Parallel processing across chromosomes
- **-resume**: Resume if interrupted

## Understanding Multi-File Output

Multi-file processing creates both individual and combined results:

### Directory Structure
```
multi_chromosome_results/
├── chr20_allele_switch_results.tsv
├── chr20_allele_switch_summary.txt  
├── chr20.corrected.vcf.gz
├── chr20.corrected.vcf.gz.tbi
├── chr21_allele_switch_results.tsv
├── chr21_allele_switch_summary.txt
├── chr21.corrected.vcf.gz
├── chr21.corrected.vcf.gz.tbi
├── chr22_allele_switch_results.tsv
├── chr22_allele_switch_summary.txt
├── chr22.corrected.vcf.gz  
├── chr22.corrected.vcf.gz.tbi
├── all_chromosomes_summary.txt      # Combined summary
└── reports/                         # Execution reports
```

### Individual Chromosome Results

Each chromosome gets its own complete analysis:

**chr20_allele_switch_summary.txt:**
```
Results Summary for chr20:
Total variants in target: 45000
Common variants: 38000  
Matched variants: 34200 (90.00%)
Switched alleles: 2660 (7.00%)
Other inconsistencies: 1140 (3.00%)
```

### Combined Summary

The `all_chromosomes_summary.txt` aggregates results:

```
Multi-Chromosome Analysis Summary
=================================

Overall Statistics:
- Chromosomes processed: 3
- Total target variants: 135000
- Total common variants: 114000
- Overall match rate: 88.60%

Per-Chromosome Breakdown:
-------------------------
chr20: 38000 common variants, 90.00% matched
chr21: 41000 common variants, 87.80% matched  
chr22: 35000 common variants, 88.00% matched

Switch Type Summary:
-------------------
Matched variants: 101004 (88.60%)
Switched alleles: 7980 (7.00%)
Complementary issues: 3420 (3.00%) 
Other inconsistencies: 1596 (1.40%)

Quality Assessment: EXCELLENT
- Match rate >85%: ✓
- Error rate <5%: ✓  
- Consistent across chromosomes: ✓
```

## Parallel Processing Benefits

Multi-file processing provides significant advantages:

### Performance Scaling
```bash
# Single file processing time
chr22 alone: ~15 minutes

# Multi-file parallel processing  
chr20 + chr21 + chr22: ~18 minutes total
# (vs 45 minutes if run sequentially)
```

### Resource Optimization
```bash
# Optimize for your system
--max_cpus 16        # Use available cores
--max_memory '64.GB' # Allocate sufficient memory  
```

## Advanced Multi-File Scenarios

### Processing Whole Genome
```bash
nextflow run main.nf \
  --targetVcfs "/data/whole_genome/chr*.vcf.gz" \
  --referenceDir "/ref/panels/" \
  --outputDir whole_genome_results \
  --fixMethod correct \
  --max_cpus 20 \
  --max_memory '128.GB' \
  -profile slurm
```

### Different Reference Panels Per Chromosome
```bash
# Organize references by chromosome
reference_panels/
├── chr20/
│   └── chr20.legend.gz
├── chr21/  
│   └── chr21.legend.gz
└── chr22/
    └── chr22.legend.gz

# CheckRef will automatically find the right references
```

### Population-Specific Processing
```bash
nextflow run main.nf \
  --targetVcfs "african_samples_chr*.vcf.gz" \
  --referenceDir "/ref/panels/african_populations/" \
  --legendPattern "*AFR*chr*.legend.gz" \
  --outputDir african_analysis \
  -profile singularity
```

## Quality Control Across Chromosomes

### Consistency Checks

Compare results across chromosomes:

```bash
# Extract match rates for each chromosome
grep "Matched variants" multi_results/*_summary.txt

# Expected output:
chr20_summary.txt:Matched variants: 34200 (90.00%)
chr21_summary.txt:Matched variants: 35958 (87.80%)
chr22_summary.txt:Matched variants: 30800 (88.00%)
```

### Identifying Problematic Chromosomes

Look for outliers:

```bash
# Check for unusually low match rates
awk '/Matched variants/ {print FILENAME, $3, $4}' *_summary.txt

# Flag chromosomes with <70% match rates
awk '/Matched variants/ {
  gsub(/[()]/, "", $4)
  if ($4 < 70) print "WARNING:", FILENAME, "low match rate:", $4
}' *_summary.txt
```

## Troubleshooting Multi-File Processing

### Common Issues

#### 1. Missing Reference Files
**Error**: `No matching legend file found for chr21`

**Solution**: Check reference directory and pattern:
```bash
ls /ref/panels/ | grep chr21
# Adjust --legendPattern if needed
--legendPattern "*chr21*.legend.gz"
```

#### 2. Inconsistent Chromosome Naming
**Error**: `Cannot match chromosome naming`

**Solution**: Standardize naming or use mapping:
```bash
# Check target chromosome names
for file in *.vcf.gz; do
  echo "$file:"
  bcftools view -H "$file" | cut -f1 | sort -u | head -3
done

# Check reference names  
ls /ref/panels/ | grep -E "(chr|Chr)"
```

#### 3. Memory Issues with Large Files
**Error**: `Process exceeded memory limit`

**Solution**: Increase resources or process fewer files:
```bash
# Increase memory allocation
--max_memory '64.GB'

# Or process in smaller batches
--targetVcfs "chr20.vcf.gz,chr21.vcf.gz"  # First batch
--targetVcfs "chr22.vcf.gz,chrX.vcf.gz"   # Second batch
```

### Performance Optimization

#### Resource Allocation
```bash
# For HPC environments
nextflow run main.nf \
  --targetVcfs "chr*.vcf.gz" \
  --referenceDir "/ref/panels/" \
  --max_cpus 32 \
  --max_memory '256.GB' \
  -profile slurm
```

#### Resume Functionality
```bash
# Resume interrupted multi-file runs
nextflow run main.nf \
  --targetVcfs "chr*.vcf.gz" \
  --referenceDir "/ref/panels/" \
  -resume  # Continues from where it left off
```

## Best Practices

### File Organization
```
project/
├── input_vcfs/
│   ├── sample_chr20.vcf.gz
│   ├── sample_chr21.vcf.gz
│   └── sample_chr22.vcf.gz
├── reference_panels/
│   ├── ref_chr20.legend.gz
│   ├── ref_chr21.legend.gz  
│   └── ref_chr22.legend.gz
└── results/
    └── multi_analysis/
```

### Naming Conventions
- Use consistent chromosome identifiers (chr20, chr21, etc.)
- Include chromosome information in filenames
- Maintain parallel structure between target and reference files

### Resource Planning
- Allocate ~4GB memory per chromosome
- Use available CPU cores for parallel processing
- Plan for ~2x disk space for output files

## Next Steps

After mastering multi-file processing:

1. **[Batch Processing](./batch-processing)** - Automate large-scale analyses
2. **[Performance Optimization](./performance)** - Speed up processing
3. **[Quality Control](./quality-control)** - Validate multi-chromosome results
4. **[Custom Parameters](./custom-parameters)** - Fine-tune for your data

## Quick Reference

### Multi-File Commands
```bash
# Comma-separated files
--targetVcfs "file1.vcf.gz,file2.vcf.gz,file3.vcf.gz"

# Glob patterns
--targetVcfs "/path/to/*.vcf.gz"
--targetVcfs "/path/to/chr*.vcf.gz"  

# Mixed paths
--targetVcfs "/path1/chr20.vcf.gz,/path2/chr21.vcf.gz"

# With resource allocation
--max_cpus 16 --max_memory '64.GB'

# Resume interrupted runs
-resume
```

### Output Files
- Individual: `chr*_allele_switch_results.tsv`, `chr*.corrected.vcf.gz`
- Combined: `all_chromosomes_summary.txt`
- Reports: `reports/execution_report.html`