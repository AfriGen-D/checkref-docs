# Multi-File Processing Tutorial

Learn how to process multiple chromosomes simultaneously in a single CheckRef run.

## What You'll Learn

- Process 3 chromosomes in one command
- Use different input methods (comma-separated vs glob patterns)
- Understand combined results
- Optimize for performance

## Time Required
15-20 minutes

## Prerequisites

- Completed the [Quick Start Tutorial](./quick-start)
- Multiple VCF files (ideally 2-3 chromosomes)
- Reference panel with matching files

## Tutorial Overview

We'll process chromosomes 20, 21, and 22 together to see how CheckRef handles multiple files simultaneously.

## Step 1: Organize Your Files

Set up your files for multi-file processing:

```bash
# Check your chromosome files
ls -la chr20.vcf.gz chr21.vcf.gz chr22.vcf.gz

# Verify they're all indexed
ls -la *.vcf.gz.tbi

# Check reference panels
ls /path/to/reference/panels/ | grep -E "(chr20|chr21|chr22)"
```

## Step 2: Choose Your Input Method

### Method A: Comma-Separated List
```bash
nextflow run main.nf \
  --targetVcfs "chr20.vcf.gz,chr21.vcf.gz,chr22.vcf.gz" \
  --referenceDir /path/to/reference/panels/ \
  --outputDir multi_comma_analysis \
  --fixMethod correct \
  -profile singularity
```

### Method B: Glob Pattern  
```bash
nextflow run main.nf \
  --targetVcfs "chr*.vcf.gz" \
  --referenceDir /path/to/reference/panels/ \
  --outputDir multi_glob_analysis \
  --fixMethod correct \
  -profile singularity
```

**Choose the method that matches your file organization.**

## Step 3: Monitor Multi-File Progress

Watch the pipeline process multiple files:

```
N E X T F L O W  ~  version 23.04.0
Launching `main.nf` [clever_pasteur] DSL2

[CHECK_ALLELE_SWITCH] Submitted process > CHECK_ALLELE_SWITCH (1) # chr20
[CHECK_ALLELE_SWITCH] Submitted process > CHECK_ALLELE_SWITCH (2) # chr21  
[CHECK_ALLELE_SWITCH] Submitted process > CHECK_ALLELE_SWITCH (3) # chr22
[CHECK_ALLELE_SWITCH] Completed process > CHECK_ALLELE_SWITCH (1)
[CHECK_ALLELE_SWITCH] Completed process > CHECK_ALLELE_SWITCH (2)
[CHECK_ALLELE_SWITCH] Completed process > CHECK_ALLELE_SWITCH (3)

[CORRECT_SWITCHED_SITES] Submitted process > CORRECT_SWITCHED_SITES (1)
[CORRECT_SWITCHED_SITES] Submitted process > CORRECT_SWITCHED_SITES (2) 
[CORRECT_SWITCHED_SITES] Submitted process > CORRECT_SWITCHED_SITES (3)
```

**Key observation:** All three chromosomes process in parallel!

## Step 4: Examine Multi-File Output

Check your comprehensive results:

```bash
# List all output files
ls -la multi_comma_analysis/

# Expected structure:
# chr20_allele_switch_results.tsv
# chr20_allele_switch_summary.txt
# chr20.corrected.vcf.gz
# chr21_allele_switch_results.tsv  
# chr21_allele_switch_summary.txt
# chr21.corrected.vcf.gz
# chr22_allele_switch_results.tsv
# chr22_allele_switch_summary.txt
# chr22.corrected.vcf.gz
# all_chromosomes_summary.txt        # Combined summary!
```

## Step 5: Review Combined Results

Look at the aggregated summary:

```bash
cat multi_comma_analysis/all_chromosomes_summary.txt
```

**Example output:**
```
Multi-Chromosome Analysis Summary
=================================

Overall Statistics:
- Chromosomes processed: 3
- Total target variants: 285000
- Total common variants: 245000
- Overall match rate: 87.2%

Per-Chromosome Breakdown:
-------------------------
chr20: 82000 common variants, 89.1% matched
chr21: 79000 common variants, 86.8% matched  
chr22: 84000 common variants, 85.7% matched

Quality Assessment: EXCELLENT
```

## Step 6: Compare Individual vs Combined

**Performance benefits:**
```bash
# Time comparison
# Individual runs: ~15 min × 3 = 45 minutes
# Multi-file run: ~18 minutes total
# Savings: 60% time reduction!
```

**Resource efficiency:**
- Shared reference loading
- Parallel processing
- Single pipeline overhead

## Step 7: Validate Results

Quick validation across all chromosomes:

```bash
# Count variants in each output
echo "=== Variant Counts ==="
for file in multi_comma_analysis/chr*.corrected.vcf.gz; do
  chr=$(basename $file | cut -d'.' -f1)
  count=$(bcftools view -H $file | wc -l)
  echo "$chr: $count variants"
done

# Check corrections applied
echo "=== Corrections Applied ==="
for file in multi_comma_analysis/chr*.corrected.vcf.gz; do
  chr=$(basename $file | cut -d'.' -f1)  
  corrections=$(bcftools view $file | grep "SWITCHED=1" | wc -l)
  echo "$chr: $corrections corrections"
done
```

## Common Multi-File Scenarios

### Scenario 1: Whole Genome Processing
```bash
# Process all autosomal chromosomes
nextflow run main.nf \
  --targetVcfs "chr[1-9].vcf.gz,chr[12][0-9].vcf.gz" \
  --referenceDir /ref/panels/ \
  --outputDir whole_genome \
  --max_cpus 20 \
  --max_memory '128.GB'
```

### Scenario 2: Specific Chromosome Set
```bash
# Process just the chromosomes you need
nextflow run main.nf \
  --targetVcfs "chr6.vcf.gz,chr19.vcf.gz,chr22.vcf.gz" \
  --referenceDir /ref/panels/ \
  --outputDir hla_gwas_regions
```

### Scenario 3: Files in Different Directories
```bash
# Mixed file locations
nextflow run main.nf \
  --targetVcfs "/data1/chr20.vcf.gz,/data2/chr21.vcf.gz,/data3/chr22.vcf.gz" \
  --referenceDir /ref/panels/ \
  --outputDir mixed_locations
```

## Troubleshooting Multi-File Issues

### Issue: Some Files Not Found
```bash
# Check each file individually
for file in chr20.vcf.gz chr21.vcf.gz chr22.vcf.gz; do
  if [ -f "$file" ]; then
    echo "✓ $file exists"
  else  
    echo "✗ $file missing"
  fi
done
```

### Issue: Inconsistent Results Across Chromosomes
```bash
# Compare match rates
grep "Matched variants" multi_comma_analysis/*_summary.txt

# Look for outliers
# chr20: 89% - Good
# chr21: 45% - Investigate!
# chr22: 87% - Good
```

### Issue: Memory Problems
```bash
# Increase resources for large multi-file runs
--max_memory '64.GB' \
--max_cpus 16
```

## Best Practices

### File Organization
```bash
# Organize files consistently
data/
├── chr20.vcf.gz
├── chr21.vcf.gz  
├── chr22.vcf.gz
└── reference_panels/
    ├── chr20.legend.gz
    ├── chr21.legend.gz
    └── chr22.legend.gz
```

### Resource Allocation
```bash
# Scale resources with file count
# 2-3 files: --max_cpus 8 --max_memory '32.GB'
# 5-10 files: --max_cpus 16 --max_memory '64.GB'  
# 15+ files: --max_cpus 32 --max_memory '128.GB'
```

### Progress Monitoring
```bash
# Use monitoring for large runs
-with-timeline timeline.html \
-with-report report.html
```

## Key Takeaways

1. **Multi-file processing is faster** than individual runs
2. **Parallel processing** handles multiple chromosomes simultaneously
3. **Combined summary** gives overall quality assessment
4. **Resource scaling** improves performance for large datasets
5. **File organization** matters for efficiency

## Next Steps

Now that you can process multiple files:

**Advanced techniques:**
- [Method Selection Tutorial](./method-selection) - Choose optimal correction strategies
- [Understanding Results](/docs/understanding-results) - Interpret multi-file output

**Quality control:**
- [Quality Control Guide](/docs/quality-control) - Verify multi-file results
- [Troubleshooting Guide](/docs/troubleshooting) - Solve multi-file issues

**For comprehensive details:**
- [Multi-File Documentation](/docs/multi-file) - Complete technical reference

## Quick Reference

### Multi-File Commands
```bash
# Comma-separated
--targetVcfs "file1.vcf.gz,file2.vcf.gz,file3.vcf.gz"

# Glob pattern
--targetVcfs "chr*.vcf.gz"

# Mixed paths  
--targetVcfs "/path1/chr20.vcf.gz,/path2/chr21.vcf.gz"

# With resources
--max_cpus 16 --max_memory '64.GB'
```

### Expected Performance
- **2-3 files**: ~60% time savings vs individual runs
- **5+ files**: ~70% time savings vs individual runs
- **Memory usage**: ~4GB per chromosome + base overhead