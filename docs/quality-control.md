# Quality Control and Validation

Comprehensive guide to CheckRef's quality control procedures, validation workflows, and automated quality assessment features.

## Overview

CheckRef includes multiple layers of quality control to ensure reliable allele switch detection and correction:

1. **Input Validation** - Pre-processing file integrity checks
2. **Analysis Quality Control** - Real-time monitoring during processing
3. **Correction Validation** - Post-processing verification of results
4. **Automated Reporting** - Comprehensive quality metrics and recommendations

## Input Validation Workflow

### Automatic VCF Validation

CheckRef automatically validates all input VCF files before processing:

```bash
# Validation runs automatically for all input files
nextflow run main.nf \
  --targetVcfs "chr20.vcf.gz,chr21.vcf.gz,chr22.vcf.gz" \
  --referenceDir /ref/panels/ \
  --outputDir results/
```

### Validation Checks Performed

| Check | Purpose | Action on Failure |
|-------|---------|-------------------|
| **File Existence** | Verify files are accessible | Skip file, continue with others |
| **File Size** | Detect empty/corrupted files | Skip file, log warning |
| **Compression Integrity** | Validate gzip compression | Skip file, suggest recompression |
| **VCF Format** | Verify bcftools can read file | Skip file, suggest format check |
| **Data Content** | Ensure variant records exist | Skip file, log warning |

### Validation Reports

Validation results are saved to `results/validation/`:

```
results/validation/
├── chr20_validation_report.txt    # Detailed validation results
├── chr21_validation_report.txt    # Per-chromosome reports
├── chr22_validation_report.txt    # Include pass/fail status
└── validation_summary.txt         # Overall validation summary
```

### Example Validation Report

```
====================================
VCF VALIDATION REPORT FOR CHR 22
====================================
File: chr22.vcf.gz
Validation Date: 2024-09-14 15:30:45

File size: 45,234,567 bytes
Data lines found: 1,234,567

✅ VALIDATION PASSED: File appears to be valid
File format: Valid VCF
Compression: gzip compressed data
Status: Ready for processing
```

## Analysis Quality Control

### Real-Time Monitoring

During analysis, CheckRef monitors several quality metrics:

#### Overlap Assessment
```bash
# CheckRef automatically calculates:
# - Target-reference overlap rates
# - Common variant counts
# - Match/switch ratios
```

#### Quality Thresholds

| Metric | Excellent | Good | Investigate |
|--------|-----------|------|-------------|
| **Target Overlap** | >70% | 50-70% | <50% |
| **Match Rate** | >85% | 70-85% | <70% |
| **Switch Rate** | <15% | 15-30% | >30% |

### Automatic Quality Flags

CheckRef flags potential quality issues:

```bash
# Common quality warnings:
⚠️  Low overlap rate (45%) - check population match
⚠️  High switch rate (35%) - verify reference panel
⚠️  No common variants found - check chromosome naming
```

## Correction Validation (Correct Method Only)

### Automated Validation Features

When using `--fixMethod correct`, CheckRef includes comprehensive validation:

```bash
nextflow run main.nf \
  --targetVcfs chr22.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod correct \
  --outputDir results/
```

### Validation Processes

1. **Switch Accuracy Verification**
   - Confirms expected switches were applied
   - Validates REF/ALT allele swaps
   - Checks SWITCHED=1 flags in output VCF

2. **Allele Frequency Validation**
   - Compares frequencies before/after correction
   - Flags unexpected frequency changes
   - Identifies potential correction errors

3. **Consistency Checks**
   - Verifies corrected variants match reference
   - Validates genotype consistency
   - Checks for introduced artifacts

### Validation Parameters

Control validation sensitivity:

```bash
nextflow run main.nf \
  --targetVcfs chr22.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod correct \
  --validation_af_threshold 0.05 \    # Flag AF changes >5%
  --validation_fold_threshold 2.0 \   # Flag 2-fold changes
  --outputDir results/
```

### Validation Output Files

```
results/verification/
├── chr22_verification_results.txt     # Pass/fail status
├── af_comparison.txt                  # Frequency analysis
├── switch_validation.txt              # Switch accuracy
└── validation_summary.html           # Interactive report
```

### Example Validation Results

**Successful Validation:**
```
✅ VERIFICATION PASSED for chr22
Switch accuracy: 98.5% (1,234 of 1,253 switches verified)
Allele frequency changes: 0.2% of variants flagged
Mean AF difference: 0.003
Maximum AF difference: 0.045
Status: High-quality corrections
```

**Quality Concerns:**
```
⚠️  VERIFICATION WARNING for chr22
Switch accuracy: 87.2% (1,095 of 1,253 switches verified)
Allele frequency changes: 3.4% of variants flagged
Mean AF difference: 0.025
Maximum AF difference: 0.156
Status: Review flagged variants manually
```

## Quality Control Best Practices

### Pre-Analysis QC

1. **File Preparation**
```bash
# Verify VCF integrity
bcftools view -h your_file.vcf.gz | head -5
bcftools stats your_file.vcf.gz > file_stats.txt

# Check chromosome naming consistency
bcftools view -H your_file.vcf.gz | cut -f1 | sort | uniq
```

2. **Reference Panel Verification**
```bash
# Verify reference panel structure
ls /ref/panels/ | grep chr22
zcat /ref/panels/chr22.legend.gz | head -5

# Check population match
# Use population-specific panels when available
```

### During Analysis QC

1. **Monitor Progress**
```bash
# Check validation reports as they're generated
tail -f results/validation/*_validation_report.txt

# Monitor Nextflow progress
nextflow log last -f name,status,exit
```

2. **Resource Monitoring**
```bash
# Check resource usage
nextflow log last -f name,cpus,memory,time
```

### Post-Analysis QC

1. **Review Validation Reports**
```bash
# Check overall success
grep "VALIDATION PASSED\|VERIFICATION PASSED" results/validation/*.txt

# Review any failures
grep "FAILED\|WARNING" results/validation/*.txt
```

2. **Quality Metrics Assessment**
```bash
# Check summary statistics
cat results/*_summary.txt

# Review correction statistics (correct method)
cat results/correction_stats.txt
```

## Troubleshooting Quality Issues

### Low Overlap Rates

**Symptoms:** <50% target-reference overlap

**Causes & Solutions:**
```bash
# Population mismatch
--referenceDir /population-specific/panels/

# Chromosome naming differences
# Check: chr22 vs 22 vs chromosome22
bcftools view -H input.vcf.gz | cut -f1 | head

# Build mismatch (hg19 vs hg38)
# Use matching reference panel build
```

### High Switch Rates

**Symptoms:** >30% variants flagged as switches

**Causes & Solutions:**
```bash
# Strand issues
# Use population-matched reference panels

# Quality filters needed
bcftools view -i 'QUAL>=20 && INFO/DP>=10' input.vcf.gz

# Wrong reference panel
# Verify panel matches your population
```

### Validation Failures

**Symptoms:** Correction validation fails

**Causes & Solutions:**
```bash
# Adjust validation thresholds
--validation_af_threshold 0.10
--validation_fold_threshold 3.0

# Review specific failures
grep "FAILED" results/verification/*.txt

# Consider remove method instead
--fixMethod remove
```

## Quality Control Workflows

### Standard QC Workflow

```bash
# 1. Run with validation enabled (default)
nextflow run main.nf \
  --targetVcfs chr22.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod correct \
  --outputDir qc_results/

# 2. Review validation reports
cat qc_results/validation/*_validation_report.txt

# 3. Check verification results
cat qc_results/verification/*_verification_results.txt

# 4. Assess overall quality
grep "✅\|⚠️\|❌" qc_results/validation/*.txt
```

### Strict QC Workflow

```bash
# Use stricter validation thresholds
nextflow run main.nf \
  --targetVcfs chr22.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod correct \
  --validation_af_threshold 0.01 \    # Flag 1% AF changes
  --validation_fold_threshold 1.5 \   # Flag 1.5-fold changes
  --outputDir strict_qc_results/
```

### Permissive QC Workflow

```bash
# Use more permissive thresholds
nextflow run main.nf \
  --targetVcfs chr22.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod correct \
  --validation_af_threshold 0.10 \    # Flag 10% AF changes
  --validation_fold_threshold 3.0 \   # Flag 3-fold changes
  --outputDir permissive_qc_results/
```

## Advanced Quality Control

### Custom Validation Scripts

Create custom validation workflows:

```bash
# Custom post-processing validation
#!/bin/bash
RESULTS_DIR=$1

# Check all chromosomes passed validation
for chr in 20 21 22; do
    if grep -q "VALIDATION PASSED" ${RESULTS_DIR}/validation/chr${chr}_validation_report.txt; then
        echo "✅ Chr${chr}: Validation passed"
    else
        echo "❌ Chr${chr}: Validation failed"
        exit 1
    fi
done

# Check correction success rates
if [ -f "${RESULTS_DIR}/correction_stats.txt" ]; then
    total_corrected=$(grep "Corrected=" ${RESULTS_DIR}/correction_stats.txt | \
                     sed 's/.*Corrected=\([0-9]*\).*/\1/' | \
                     awk '{sum+=$1} END {print sum}')
    echo "Total corrections: ${total_corrected}"
fi
```

### Integration with External QC Tools

```bash
# Integrate with bcftools stats
bcftools stats results/chr22.corrected.vcf.gz > chr22_stats.txt

# Compare with original
bcftools stats original/chr22.vcf.gz > chr22_original_stats.txt
diff chr22_original_stats.txt chr22_stats.txt

# Use VCFtools for additional QC
vcftools --gzvcf results/chr22.corrected.vcf.gz \
         --freq --out chr22_corrected_freq

# Hardy-Weinberg equilibrium testing
vcftools --gzvcf results/chr22.corrected.vcf.gz \
         --hardy --out chr22_corrected_hwe
```

### Automated QC Reporting

```bash
# Generate comprehensive QC report
#!/bin/bash
RESULTS_DIR=$1
REPORT_FILE="${RESULTS_DIR}/qc_summary_report.html"

cat > ${REPORT_FILE} << EOF
<!DOCTYPE html>
<html>
<head><title>CheckRef QC Report</title></head>
<body>
<h1>CheckRef Quality Control Report</h1>
<h2>Validation Summary</h2>
EOF

# Add validation results
for report in ${RESULTS_DIR}/validation/*_validation_report.txt; do
    chr=$(basename $report | sed 's/_validation_report.txt//')
    if grep -q "VALIDATION PASSED" $report; then
        echo "<p>✅ ${chr}: Passed</p>" >> ${REPORT_FILE}
    else
        echo "<p>❌ ${chr}: Failed</p>" >> ${REPORT_FILE}
    fi
done

echo "</body></html>" >> ${REPORT_FILE}
```

## Quality Control Parameters Reference

### Validation Control Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--skip_validation` | boolean | `false` | Skip all validation steps |
| `--validation_af_threshold` | float | `0.05` | AF change threshold for flagging |
| `--validation_fold_threshold` | float | `2.0` | Fold change threshold for flagging |

### Quality Thresholds

| Metric | Parameter | Recommended Range |
|--------|-----------|-------------------|
| **AF Change Threshold** | `--validation_af_threshold` | 0.01-0.10 |
| **Fold Change Threshold** | `--validation_fold_threshold` | 1.5-3.0 |

## Next Steps

After implementing quality control:

1. **[Understanding Results](./understanding-results)** - Interpret QC metrics
2. **[Troubleshooting](./troubleshooting)** - Solve quality issues
3. **[Multi-File Processing](./multi-file)** - Scale QC to multiple chromosomes
4. **[Correction Methods](./correction-methods)** - Choose optimal correction approach

## Quick Reference

### Essential QC Commands

```bash
# Check validation status
grep "VALIDATION\|VERIFICATION" results/validation/*.txt

# Review quality metrics
cat results/*_summary.txt

# Check correction success (correct method)
cat results/correction_stats.txt

# Identify quality issues
grep "WARNING\|FAILED" results/validation/*.txt
```

### Quality Thresholds Summary

- **Excellent Quality**: >95% switch accuracy, <1% AF changes
- **Good Quality**: 90-95% switch accuracy, 1-5% AF changes
- **Investigate**: <90% switch accuracy, >5% AF changes
