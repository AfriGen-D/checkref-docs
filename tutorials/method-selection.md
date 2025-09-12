# Method Selection Tutorial

Learn how to choose between CheckRef's "remove" and "correct" methods for your specific analysis needs.

## What You'll Learn

- When to use remove vs correct methods
- How to test both methods on your data
- How to evaluate which method works better
- How to make informed decisions for different analysis types

## Time Required
20 minutes

## Prerequisites

- Completed [Quick Start Tutorial](./quick-start)
- Understanding of your downstream analysis goals
- Sample VCF file for testing

## The Two Methods

CheckRef offers two approaches for handling allele switches:

| Method | Action | Best For |
|--------|--------|----------|
| **remove** | Delete switched variants | Quality-first, imputation prep |
| **correct** | Fix by swapping REF/ALT | Power-first, association studies |

## Tutorial: Testing Both Methods

We'll run both methods on the same file and compare results.

### Step 1: Run Remove Method

```bash
nextflow run main.nf \
  --targetVcfs test_chr22.vcf.gz \
  --referenceDir /path/to/reference/panels/ \
  --fixMethod remove \
  --outputDir method_test_remove \
  -profile singularity
```

### Step 2: Run Correct Method

```bash
nextflow run main.nf \
  --targetVcfs test_chr22.vcf.gz \
  --referenceDir /path/to/reference/panels/ \
  --fixMethod correct \
  --outputDir method_test_correct \
  -profile singularity
```

### Step 3: Compare the Results

Create a comparison script:

```bash
cat > compare_methods.sh << 'EOF'
#!/bin/bash

echo "=== Method Comparison ==="
echo "Analysis: $(basename $(pwd))"
echo "Date: $(date)"
echo "========================"

# Original file stats
original_variants=$(bcftools view -H test_chr22.vcf.gz | wc -l)
echo "Original variants: $original_variants"

# Remove method results
if [ -f method_test_remove/chr22.noswitch.vcf.gz ]; then
  remove_variants=$(bcftools view -H method_test_remove/chr22.noswitch.vcf.gz | wc -l)
  remove_retention=$(echo "scale=1; ($remove_variants * 100) / $original_variants" | bc)
  echo "Remove method: $remove_variants variants (${remove_retention}% retention)"
fi

# Correct method results  
if [ -f method_test_correct/chr22.corrected.vcf.gz ]; then
  correct_variants=$(bcftools view -H method_test_correct/chr22.corrected.vcf.gz | wc -l)
  correct_retention=$(echo "scale=1; ($correct_variants * 100) / $original_variants" | bc)
  corrections=$(bcftools view method_test_correct/chr22.corrected.vcf.gz | grep "SWITCHED=1" | wc -l)
  echo "Correct method: $correct_variants variants (${correct_retention}% retention)"
  echo "Corrections applied: $corrections"
fi

# Quality comparison
echo ""
echo "=== Quality Metrics ==="
for method in remove correct; do
  summary_file="method_test_${method}/chr22_allele_switch_summary.txt"
  if [ -f "$summary_file" ]; then
    echo "${method^} method:"
    match_rate=$(grep "Matched variants" $summary_file | awk '{print $3}' | tr -d '()')
    error_rate=$(grep "Other inconsistencies" $summary_file | awk '{print $3}' | tr -d '()')
    echo "  Match rate: $match_rate"
    echo "  Error rate: $error_rate"
  fi
done
EOF

chmod +x compare_methods.sh
./compare_methods.sh
```

### Step 4: Analyze the Comparison

**Example output:**
```
=== Method Comparison ===
Original variants: 125,000

Remove method: 108,000 variants (86.4% retention)
Correct method: 123,500 variants (98.8% retention)
Corrections applied: 16,500

=== Quality Metrics ===
Remove method:
  Match rate: 100.0%  (only perfect matches remain)
  Error rate: 0.0%    (all problems removed)

Correct method:  
  Match rate: 87.5%   (includes corrected variants)
  Error rate: 1.2%    (only uncorrectable problems remain)
```

## Decision Framework

Use this framework to choose the best method:

### Choose REMOVE when:

**Preparing for imputation:**
```bash
# Imputation servers require high-quality variants
nextflow run main.nf \
  --targetVcfs pre_imputation.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod remove
```

**Quality is more important than quantity:**
```bash
# Fine-mapping studies, functional analyses
nextflow run main.nf \
  --targetVcfs gwas_hits.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod remove
```

**Conservative approach preferred:**
- When you can't validate corrections
- For novel or rare variants
- When downstream tools are sensitive to orientation issues

### Choose CORRECT when:

**Association studies needing power:**
```bash
# GWAS where variant count matters
nextflow run main.nf \
  --targetVcfs gwas_cohort.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod correct
```

**Rare variant analysis:**
```bash
# Every variant counts
nextflow run main.nf \
  --targetVcfs exome_rare_variants.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod correct
```

**Maximum data retention needed:**
- Large-scale population studies
- When you can validate corrections
- For well-characterized variant sets

## Real-World Scenarios

### Scenario 1: GWAS Preparation

**Your goal:** Genome-wide association study
**Key need:** Statistical power (many variants)
**Recommendation:** CORRECT method

```bash
nextflow run main.nf \
  --targetVcfs "gwas_chr*.vcf.gz" \
  --referenceDir /ref/panels/ \
  --fixMethod correct \
  --outputDir gwas_corrected
```

**Why:** GWAS benefits from maximum variant count. The `SWITCHED=1` flag lets you track corrected variants if needed.

### Scenario 2: Imputation Server Submission

**Your goal:** Submit to Michigan Imputation Server
**Key need:** Data quality (server requirements)
**Recommendation:** REMOVE method

```bash
nextflow run main.nf \
  --targetVcfs "pre_imputation_chr*.vcf.gz" \
  --referenceDir /ref/panels/ \
  --fixMethod remove \
  --outputDir imputation_ready
```

**Why:** Imputation servers are strict about data quality. Removing problematic variants prevents server rejection.

### Scenario 3: Rare Disease Study

**Your goal:** Exome analysis for rare disease
**Key need:** Retain functional variants
**Recommendation:** CORRECT method + validation

```bash
nextflow run main.nf \
  --targetVcfs exome_family.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod correct \
  --outputDir rare_disease_corrected

# Then validate corrections for coding variants
bcftools view rare_disease_corrected/exome_family.corrected.vcf.gz | \
  grep "SWITCHED=1" | \
  bcftools query -f '%CHROM:%POS %REF>%ALT\n'
```

### Scenario 4: Population Genetics

**Your goal:** Population structure analysis
**Key need:** Allele frequency accuracy
**Recommendation:** REMOVE method

```bash
nextflow run main.nf \
  --targetVcfs "population_chr*.vcf.gz" \
  --referenceDir /ref/panels/ \
  --fixMethod remove \
  --outputDir population_clean
```

**Why:** Population genetics analyses are sensitive to allele orientation. Clean data is more important than variant count.

## Advanced Decision Making

### Hybrid Approach

You can use both methods strategically:

```bash
# Step 1: Get comprehensive results with correct method
nextflow run main.nf \
  --targetVcfs input.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod correct \
  --outputDir full_correct

# Step 2: Create high-confidence subset for critical analyses
bcftools view -e 'INFO/SWITCHED=1' \
  full_correct/input.corrected.vcf.gz > high_confidence.vcf.gz
```

### Method-Specific Quality Filters

**For remove method - be less strict on input:**
```bash
# Keep more variants initially, let CheckRef do the filtering
bcftools view -i 'QUAL>=10' input.vcf.gz | \
nextflow run main.nf \
  --targetVcfs /dev/stdin \
  --referenceDir /ref/panels/ \
  --fixMethod remove
```

**For correct method - be more strict on input:**
```bash
# Apply strict filters, then correct remaining issues
bcftools view -i 'QUAL>=30 && INFO/DP>=20' input.vcf.gz > filtered.vcf.gz
nextflow run main.nf \
  --targetVcfs filtered.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod correct
```

## Validation Strategies

### For Remove Method
```bash
# Check retention rates by variant type
bcftools view -H method_test_remove/chr22.noswitch.vcf.gz | \
  awk '{print length($4), length($5)}' | \
  awk '$1==1 && $2==1 {snp++} $1>1 || $2>1 {indel++} 
       END {print "SNPs:", snp, "INDELs:", indel}'
```

### For Correct Method  
```bash
# Validate a sample of corrections
bcftools view method_test_correct/chr22.corrected.vcf.gz | \
  grep "SWITCHED=1" | head -5 | \
  cut -f1,2,4,5

# Check these positions in reference panel
# Manually verify corrections make sense
```

## Common Decision Points

### High Switch Rate (>20%)
- **Consider:** Population mismatch with reference panel
- **Try:** Different reference panel first
- **Then:** Use remove method for safety

### Low Match Rate (<70%)
- **Consider:** Data quality issues
- **Try:** More strict input filtering
- **Then:** Use remove method initially

### Critical Analysis
- **Use:** Remove method for initial analysis
- **Then:** Compare with correct method results
- **Validate:** Key findings with both methods

## Method Performance

### Typical Results

**Remove method:**
- Retention: 80-90%
- Quality: 100% match rate
- Speed: Slightly faster
- Downstream: Universal compatibility

**Correct method:**
- Retention: 95-99%
- Quality: 85-95% match rate  
- Speed: Slightly slower
- Downstream: May need `SWITCHED=1` handling

## Next Steps

After choosing your method:

**For further learning:**
- [Results Validation Tutorial](./results-validation) - Verify your chosen method worked
- [QC Setup Tutorial](./qc-setup) - Monitor method performance over time

**For comprehensive details:**
- [Correction Methods Documentation](/docs/correction-methods) - Complete technical reference

## Quick Reference

### Decision Tree
```
Need maximum quality? → REMOVE
Need maximum variants? → CORRECT
Preparing for imputation? → REMOVE  
Running association study? → CORRECT
Can't validate corrections? → REMOVE
Have validation resources? → CORRECT
```

### Method Commands
```bash
# Remove method
--fixMethod remove

# Correct method
--fixMethod correct

# Check corrections
bcftools view output.vcf.gz | grep "SWITCHED=1" | wc -l
```

### Validation Commands
```bash
# Compare retention
echo "Original: $(bcftools view -H input.vcf.gz | wc -l)"
echo "Final: $(bcftools view -H output.vcf.gz | wc -l)"

# Check correction rate
corrections=$(bcftools view output.vcf.gz | grep "SWITCHED=1" | wc -l)
switches=$(grep -c "SWITCH" results.tsv)
echo "Correction success: $corrections/$switches"
```