# Single File Analysis Tutorial

Master the fundamentals of CheckRef with detailed single-file analysis.

## Learning Objectives

- Perform thorough single-chromosome analysis
- Understand detailed parameter options
- Learn quality assessment techniques
- Master result interpretation for single files

## Overview

Single file analysis is the foundation of CheckRef usage. This tutorial provides an in-depth walkthrough of analyzing one VCF file against a reference panel, with detailed explanations of each step.

## Preparing Your Single File Analysis

### Input Requirements

**Target VCF File:**
- Compressed with bgzip (`.vcf.gz`)
- Tabix indexed (`.vcf.gz.tbi`)
- Contains SNP variants
- Proper chromosome naming

**Reference Legend File:**
- Standard legend format
- Matching chromosome identifier
- Contains allele information (a0, a1)

### Example Setup

```bash
# Verify your input file
ls -la chr22.vcf.gz chr22.vcf.gz.tbi

# Check VCF contents
bcftools view -h chr22.vcf.gz | head -5
bcftools view -H chr22.vcf.gz | head -3

# Verify reference panel
ls -la /ref/panels/*chr22*.legend.gz
zcat /ref/panels/chr22.legend.gz | head -5
```

## Step-by-Step Analysis

### Step 1: Basic Analysis Command

```bash
nextflow run main.nf \
  --targetVcfs chr22.vcf.gz \
  --referenceDir /path/to/reference/panels/ \
  --outputDir chr22_analysis \
  --fixMethod remove \
  -profile singularity
```

### Step 2: Monitor Progress

CheckRef provides progress updates during execution:

```
N E X T F L O W  ~  version 23.04.0
Launching `main.nf` [peaceful_sinoussi] DSL2 - revision: abc123ef

[CHECK_ALLELE_SWITCH] Submitted process > CHECK_ALLELE_SWITCH (1)
[CHECK_ALLELE_SWITCH] Running process > CHECK_ALLELE_SWITCH (1)
[CHECK_ALLELE_SWITCH] Completed process > CHECK_ALLELE_SWITCH (1)
[REMOVE_SWITCHED_SITES] Submitted process > REMOVE_SWITCHED_SITES (1)
[REMOVE_SWITCHED_SITES] Completed process > REMOVE_SWITCHED_SITES (1)

Pipeline completed successfully!
```

### Step 3: Examine Output Structure

```bash
# Check output directory
ls -la chr22_analysis/

# Expected files:
chr22_analysis/
├── chr22_allele_switch_results.tsv    # Detailed results
├── chr22_allele_switch_summary.txt    # Summary statistics
├── chr22.noswitch.vcf.gz              # Cleaned VCF
├── chr22.noswitch.vcf.gz.tbi          # VCF index
└── reports/                           # Execution reports
    ├── execution_report.html
    ├── timeline_report.html
    └── dag_report.html
```

## Detailed Result Analysis

### Examining the Switch Results File

```bash
# Look at the detailed results
head -20 chr22_analysis/chr22_allele_switch_results.tsv

# Count different switch types
cut -f7 chr22_analysis/chr22_allele_switch_results.tsv | sort | uniq -c

# Example output:
#  4250 MATCH
#   350 SWITCH  
#   200 COMPLEMENT
#    75 COMPLEMENT_SWITCH
#    25 OTHER
```

### Understanding Each Switch Type

**MATCH entries:**
```
CHROM  POS      TARGET_REF TARGET_ALT REF_REF REF_ALT STATUS
22     16050075 A          G          A       G       MATCH
```
- Perfect alignment between target and reference
- No action needed
- High confidence variants

**SWITCH entries:**
```  
CHROM  POS      TARGET_REF TARGET_ALT REF_REF REF_ALT STATUS
22     16050115 G          A          A       G       SWITCH
```
- REF and ALT are flipped
- Easily correctable
- Common occurrence

**COMPLEMENT entries:**
```
CHROM  POS      TARGET_REF TARGET_ALT REF_REF REF_ALT STATUS  
22     16050213 C          T          G       A       COMPLEMENT
```
- Opposite strand representation
- Requires strand flip correction
- More complex but manageable

### Summary Statistics Interpretation

```bash
# Review the summary file
cat chr22_analysis/chr22_allele_switch_summary.txt

# Example content:
Results Summary:
Total variants in target: 125,000
Total variants in reference: 1,741,597
Common variants: 4,900
Matched variants: 4,250 (86.73%)
Switched alleles: 350 (7.14%)
Complementary strand issues: 200 (4.08%)
Complement + switch issues: 75 (1.53%)
Other inconsistencies: 25 (0.51%)
Target overlap: 3.92%
Reference overlap: 0.28%
```

**Key Metrics Explained:**

1. **Target overlap (3.92%)**: 3.92% of your variants are in the reference panel
2. **Match rate (86.73%)**: High quality - most overlapping variants align well
3. **Switch rate (7.14%)**: Normal range - indicates minor orientation issues
4. **Error rate (0.51%)**: Very low - suggests good data quality

## Quality Assessment

### Assessing Results Quality

**Excellent Quality Indicators:**
- Match rate >85%
- Switch rate 5-15%
- Error rate <5%
- Target overlap >50% (for well-covered chromosomes)

**Good Quality Indicators:**
- Match rate 70-85%
- Switch rate 10-25%
- Error rate 5-10%
- Target overlap 30-50%

**Concerning Indicators:**
- Match rate <70%
- Error rate >10%
- Very low target overlap <10%

### Detailed Quality Checks

```bash
# Check variant distribution across positions
bcftools view -H chr22.vcf.gz | cut -f2 | \
  awk '{print int($1/1000000)"M-"int($1/1000000+1)"M"}' | \
  sort | uniq -c

# Check allele frequency distribution
bcftools +fill-tags chr22.vcf.gz -- -t AF | \
  bcftools query -f '%AF\n' | \
  awk 'BEGIN{rare=0; common=0} $1<0.05{rare++} $1>=0.05{common++} END{print "Rare:", rare, "Common:", common}'

# Examine problematic variants
grep "OTHER" chr22_analysis/chr22_allele_switch_results.tsv | head -10
```

## Advanced Single File Analysis

### Custom Parameter Usage

```bash
# High-memory analysis for large files
nextflow run main.nf \
  --targetVcfs large_chr22.vcf.gz \
  --referenceDir /ref/panels/ \
  --outputDir chr22_advanced \
  --max_memory '32.GB' \
  --max_cpus 8 \
  --fixMethod correct \
  -profile singularity

# Specific reference pattern matching  
nextflow run main.nf \
  --targetVcfs chr22.vcf.gz \
  --referenceDir /ref/panels/ \
  --legendPattern "*chr22*v2*.legend.gz" \
  --outputDir chr22_specific \
  -profile docker
```

### Quality Control Steps

```bash
# Pre-analysis QC
bcftools stats chr22.vcf.gz > pre_analysis_stats.txt

# Check for duplicate variants
bcftools view -H chr22.vcf.gz | cut -f1,2 | sort | uniq -d

# Verify chromosome consistency
bcftools view -H chr22.vcf.gz | cut -f1 | sort -u
# Should show only "22" or "chr22"

# Post-analysis validation
bcftools view -H chr22_analysis/chr22.noswitch.vcf.gz | wc -l
# Compare with original count
```

### Comparing Remove vs Correct Methods

Run both methods on the same file:

```bash
# Remove method
nextflow run main.nf \
  --targetVcfs chr22.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod remove \
  --outputDir chr22_remove

# Correct method
nextflow run main.nf \
  --targetVcfs chr22.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod correct \
  --outputDir chr22_correct

# Compare results
echo "Original variants:"
bcftools view -H chr22.vcf.gz | wc -l

echo "Remove method retained:"
bcftools view -H chr22_remove/chr22.noswitch.vcf.gz | wc -l

echo "Correct method retained:"  
bcftools view -H chr22_correct/chr22.corrected.vcf.gz | wc -l

echo "Corrections applied:"
bcftools view chr22_correct/chr22.corrected.vcf.gz | grep "SWITCHED=1" | wc -l
```

## Common Single File Scenarios

### Scenario 1: High-Quality Array Data

Typical results for well-processed SNP array data:

```bash
# Expected summary for array data:
# Match rate: 85-95%
# Switch rate: 3-8% 
# Error rate: <2%
# Target overlap: 60-80%

# Recommended approach:
--fixMethod correct  # Retain maximum variants
```

### Scenario 2: Whole Genome Sequencing Data

Results for WGS data may vary more:

```bash
# Expected summary for WGS:
# Match rate: 70-90%
# Switch rate: 5-15%
# Error rate: 2-8%
# Target overlap: 40-70%

# Recommended approach:
--fixMethod remove   # Conservative for novel variants
```

### Scenario 3: Targeted Sequencing

Results for exome or targeted sequencing:

```bash
# Expected summary:
# Match rate: 75-90%  
# Switch rate: 5-12%
# Error rate: 3-10%
# Target overlap: 30-60%

# Approach depends on downstream analysis needs
```

## Troubleshooting Single File Issues

### Low Match Rates

```bash
# Investigate low match rates
grep -v "MATCH" chr22_analysis/chr22_allele_switch_results.tsv | \
  cut -f7 | sort | uniq -c

# Check if reference panel is appropriate
ls /ref/panels/ | grep -i population
# Use population-matched reference if available
```

### High Error Rates

```bash
# Examine problematic variants
grep "OTHER" chr22_analysis/chr22_allele_switch_results.tsv | \
  head -20

# Check if variants are complex
bcftools view -H chr22.vcf.gz | \
  awk 'length($4)>1 || length($5)>1 {print $1":"$2, $4">"$5}' | \
  head -10

# Filter to SNPs only if needed
bcftools view -v snps chr22.vcf.gz > chr22_snps_only.vcf.gz
```

### Unexpected Results

```bash
# Validate against known databases
# Check a few variants manually
bcftools query -f '%CHROM:%POS %REF>%ALT\n' chr22.vcf.gz | head -5

# Look up in dbSNP or other databases to verify allele orientations
```

## Best Practices for Single File Analysis

### File Preparation
1. Always validate VCF format before analysis
2. Ensure proper indexing with tabix
3. Use consistent chromosome naming
4. Apply basic quality filters if needed

### Parameter Selection
1. Start with default parameters
2. Adjust memory based on file size
3. Choose fix method based on downstream needs
4. Use appropriate container profile

### Result Validation
1. Always review summary statistics
2. Check a sample of problematic variants manually
3. Validate corrections when using correct method
4. Compare with expected patterns for your data type

## Next Steps

After mastering single file analysis:

1. **[Multi-File Processing](./multi-file)** - Scale to multiple chromosomes
2. **[Understanding Results](./understanding-results)** - Deep dive into interpretation  
3. **[Correction Methods](./correction-methods)** - Choose optimal fix strategies
4. **[Quality Control](./quality-control)** - Validate your results

## Quick Reference

### Single File Command Template
```bash
nextflow run main.nf \
  --targetVcfs YOURFILE.vcf.gz \
  --referenceDir /path/to/reference/panels/ \
  --outputDir OUTPUTNAME \
  --fixMethod [remove|correct] \
  [--max_memory MEMORY] \
  [--max_cpus CPUS] \
  -profile [docker|singularity|standard]
```

### Quality Thresholds
- **Match rate**: >85% excellent, 70-85% good, <70% investigate
- **Switch rate**: 5-15% normal, >25% investigate  
- **Error rate**: <5% excellent, 5-10% acceptable, >10% investigate