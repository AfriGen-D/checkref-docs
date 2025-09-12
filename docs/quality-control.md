# Quality Control Tutorial

Learn comprehensive quality control strategies for validating CheckRef results and ensuring data integrity.

## Learning Objectives

By the end of this tutorial, you'll be able to:
- Implement systematic QC checks for CheckRef results
- Identify and investigate quality issues
- Validate corrections and removals
- Establish QC pipelines for routine use

## QC Overview

Quality control for CheckRef involves multiple validation layers:

1. **Input validation** - Verify data quality before processing
2. **Process monitoring** - Track pipeline execution 
3. **Result validation** - Confirm output quality
4. **Correction verification** - Validate fixes were applied correctly

## Pre-Analysis QC

### Input File Validation

Before running CheckRef, validate your input files:

```bash
# Check VCF file integrity
bcftools view -h input.vcf.gz | head -10
bcftools stats input.vcf.gz > input_stats.txt

# Verify file format and compression
file input.vcf.gz
tabix -l input.vcf.gz

# Check for common issues
bcftools view -H input.vcf.gz | head -5 | cut -f1-5
```

**Key checks:**
- File is properly bgzip compressed
- Tabix index exists and is valid
- VCF header contains required fields
- Variants are properly formatted
- Chromosome naming is consistent

### Reference Panel Validation

```bash
# Verify reference files exist and are readable
ls -la /ref/panels/*.legend.gz

# Check reference format
zcat /ref/panels/chr22.legend.gz | head -5

# Validate chromosome naming consistency
zcat /ref/panels/chr22.legend.gz | cut -f1 | head -1000 | \
  grep -o "chr[0-9XY]*\|^[0-9XY]*" | sort -u
```

### Expected Input Quality Metrics

```bash
# Generate comprehensive input statistics
bcftools stats input.vcf.gz

# Key metrics to check:
# - Number of SNPs vs INDELs (should be mostly SNPs)
# - Ti/Tv ratio (should be ~2.0-2.1 for human data)
# - Missing data rate (should be <10%)
# - Allele frequency distribution
```

## Process Monitoring QC

### Real-Time Monitoring

Track pipeline progress and resource usage:

```bash
# Run with monitoring enabled
nextflow run main.nf \
  --targetVcfs input.vcf.gz \
  --referenceDir /ref/panels/ \
  --outputDir qc_analysis \
  -with-trace trace.txt \
  -with-timeline timeline.html \
  -with-report report.html

# Monitor progress in another terminal
tail -f .nextflow.log

# Check resource usage
nextflow log last -f name,status,exit,realtime,rss
```

### Process-Level QC

```bash
# Check individual process success
grep "CHECK_ALLELE_SWITCH" trace.txt
grep "REMOVE_SWITCHED_SITES\|CORRECT_SWITCHED_SITES" trace.txt

# Monitor memory and CPU usage
awk -F'\t' 'NR>1 {print $2, $10, $11}' trace.txt | \
  sort -k2 -nr | head -10
```

## Result Validation QC

### Summary Statistics QC

Establish quality thresholds for summary statistics:

```bash
# Extract key metrics from summary
grep -E "(Total variants|Matched variants|Other inconsistencies)" \
  qc_analysis/*_summary.txt

# Calculate quality scores
awk '
/Matched variants/ { match_rate = $3; gsub(/[()]/, "", match_rate) }
/Other inconsistencies/ { error_rate = $3; gsub(/[()]/, "", error_rate) }
END { 
  print "Match rate:", match_rate
  print "Error rate:", error_rate
  if (match_rate > 85) print "PASS: High match rate"
  else if (match_rate > 70) print "WARN: Moderate match rate"  
  else print "FAIL: Low match rate"
  
  if (error_rate < 5) print "PASS: Low error rate"
  else if (error_rate < 10) print "WARN: Moderate error rate"
  else print "FAIL: High error rate"
}' qc_analysis/*_summary.txt
```

### Detailed Results QC

```bash
# Validate switch classification distribution
cut -f7 qc_analysis/*_allele_switch_results.tsv | sort | uniq -c | \
  awk '{
    total += $1; 
    type[$2] = $1
  } 
  END {
    for (t in type) {
      pct = (type[t]/total)*100
      printf "%-20s: %6d (%5.1f%%)\n", t, type[t], pct
    }
    
    # Quality assessment
    match_pct = (type["MATCH"]/total)*100
    other_pct = (type["OTHER"]/total)*100
    
    if (match_pct > 80) print "✓ PASS: Match rate > 80%"
    else print "✗ WARN: Low match rate"
    
    if (other_pct < 10) print "✓ PASS: Error rate < 10%" 
    else print "✗ WARN: High error rate"
  }'
```

### Cross-Validation Checks

Compare results with expected patterns:

```bash
# Check allele frequency consistency
bcftools +fill-tags input.vcf.gz -- -t AF,AC | \
  bcftools query -f '%CHROM:%POS\t%AF\n' > input_af.txt

bcftools +fill-tags qc_analysis/*.corrected.vcf.gz -- -t AF,AC | \
  bcftools query -f '%CHROM:%POS\t%AF\n' > output_af.txt

# Compare allele frequencies for corrected sites
join -t$'\t' input_af.txt output_af.txt | \
  awk -F'\t' '$2 != $3 {print $1, $2, $3}' | head -10
```

## Correction Validation QC

### Verify Correct Method Results

When using the correct method, validate that corrections were applied properly:

```bash
# Count expected vs actual corrections
expected=$(grep -E "SWITCH|COMPLEMENT" qc_analysis/*_results.tsv | wc -l)
actual=$(bcftools view qc_analysis/*.corrected.vcf.gz | grep "SWITCHED=1" | wc -l)

echo "Expected corrections: $expected"
echo "Actual corrections: $actual"

if [ "$expected" -eq "$actual" ]; then
  echo "✓ PASS: All expected corrections applied"
else
  echo "✗ WARN: Correction count mismatch"
fi
```

### Manual Verification of Corrections

```bash
# Select a few corrected variants for manual verification
bcftools view qc_analysis/*.corrected.vcf.gz | \
  grep "SWITCHED=1" | head -5 | \
  cut -f1,2,4,5

# Check these positions in the original file
bcftools view input.vcf.gz 22:16050075 | grep -v "^#"

# Check against reference panel
zcat /ref/panels/chr22.legend.gz | grep "16050075"
```

### Validate Remove Method Results

For the remove method, confirm appropriate variants were removed:

```bash
# Count removed variants
original=$(bcftools view -H input.vcf.gz | wc -l)
remaining=$(bcftools view -H qc_analysis/*.noswitch.vcf.gz | wc -l)
removed=$((original - remaining))

echo "Original variants: $original"
echo "Remaining variants: $remaining" 
echo "Removed variants: $removed"

# Check removal rate
removal_rate=$(echo "scale=2; ($removed * 100) / $original" | bc)
echo "Removal rate: ${removal_rate}%"

# Validate removed positions match switch results
grep -v "MATCH" qc_analysis/*_results.tsv | wc -l
```

## Multi-File QC

When processing multiple files, implement comparative QC:

### Cross-Chromosome Consistency

```bash
# Compare match rates across chromosomes
for file in results/*_summary.txt; do
  chr=$(basename $file | cut -d'_' -f1)
  match_rate=$(grep "Matched variants" $file | awk '{print $3}' | tr -d '()')
  echo "$chr: $match_rate"
done | sort -k2 -nr

# Flag outlier chromosomes
awk '{
  chr[NR] = $1; rate[NR] = $2; sum += $2
} 
END {
  avg = sum/NR
  for (i=1; i<=NR; i++) {
    diff = (rate[i] - avg)
    if (diff < 0) diff = -diff
    if (diff > 15) print "OUTLIER:", chr[i], rate[i]"% (avg:", avg"%)"
  }
}' chromosome_rates.txt
```

### Cohort-Level QC

```bash
# Generate cohort summary
cat > cohort_qc.sh << 'EOF'
#!/bin/bash

total_variants=0
total_common=0
total_matched=0

for summary in results/*_summary.txt; do
  variants=$(grep "Total variants in target" $summary | awk '{print $5}')
  common=$(grep "Common variants" $summary | awk '{print $3}') 
  matched=$(grep "Matched variants" $summary | awk '{print $3}')
  
  total_variants=$((total_variants + variants))
  total_common=$((total_common + common))
  total_matched=$((total_matched + matched))
done

echo "Cohort Summary:"
echo "Total variants processed: $total_variants"
echo "Total common variants: $total_common"
echo "Total matched variants: $total_matched"
echo "Overall match rate: $(echo "scale=2; ($total_matched * 100) / $total_common" | bc)%"
EOF

chmod +x cohort_qc.sh
./cohort_qc.sh
```

## Automated QC Pipeline

### QC Script Template

Create an automated QC script:

```bash
cat > checkref_qc.sh << 'EOF'
#!/bin/bash

# CheckRef Quality Control Script
set -e

RESULTS_DIR=$1
ORIGINAL_VCF=$2

if [ $# -ne 2 ]; then
  echo "Usage: $0 <results_dir> <original_vcf>"
  exit 1
fi

echo "=== CheckRef QC Report ==="
echo "Analysis: $(basename $RESULTS_DIR)"  
echo "Date: $(date)"
echo "=========================="

# Check files exist
if [ ! -d "$RESULTS_DIR" ]; then
  echo "❌ Results directory not found: $RESULTS_DIR"
  exit 1
fi

if [ ! -f "$ORIGINAL_VCF" ]; then
  echo "❌ Original VCF not found: $ORIGINAL_VCF"
  exit 1
fi

# Parse summary statistics
summary_file=$(find $RESULTS_DIR -name "*_summary.txt" | head -1)

if [ -f "$summary_file" ]; then
  echo "📊 Summary Statistics:"
  match_rate=$(grep "Matched variants" $summary_file | awk '{print $3}' | tr -d '()')
  error_rate=$(grep "Other inconsistencies" $summary_file | awk '{print $3}' | tr -d '()')
  
  echo "  Match rate: $match_rate"
  echo "  Error rate: $error_rate"
  
  # Quality assessment
  match_num=$(echo $match_rate | tr -d '%')
  error_num=$(echo $error_rate | tr -d '%')
  
  if (( $(echo "$match_num > 80" | bc -l) )); then
    echo "  ✅ Match rate: EXCELLENT"
  elif (( $(echo "$match_num > 70" | bc -l) )); then
    echo "  ⚠️  Match rate: GOOD"
  else
    echo "  ❌ Match rate: POOR"
  fi
  
  if (( $(echo "$error_num < 5" | bc -l) )); then
    echo "  ✅ Error rate: EXCELLENT"
  elif (( $(echo "$error_num < 10" | bc -l) )); then
    echo "  ⚠️  Error rate: ACCEPTABLE"
  else
    echo "  ❌ Error rate: HIGH"
  fi
fi

# Check output files
echo ""
echo "📁 Output Files:"
for ext in "noswitch.vcf.gz" "corrected.vcf.gz"; do
  output_file=$(find $RESULTS_DIR -name "*.$ext" | head -1)
  if [ -f "$output_file" ]; then
    size=$(stat -f%z "$output_file" 2>/dev/null || stat -c%s "$output_file")
    echo "  ✅ $(basename $output_file): ${size} bytes"
    
    # Count variants
    count=$(bcftools view -H "$output_file" | wc -l)
    echo "     Variants: $count"
  fi
done

# Correction validation (if applicable)
corrected_file=$(find $RESULTS_DIR -name "*.corrected.vcf.gz" | head -1)
if [ -f "$corrected_file" ]; then
  echo ""
  echo "🔧 Correction Validation:"
  corrections=$(bcftools view "$corrected_file" | grep "SWITCHED=1" | wc -l)
  echo "  Corrections applied: $corrections"
  
  if [ $corrections -gt 0 ]; then
    echo "  ✅ Corrections detected in output"
  fi
fi

echo ""
echo "=========================="
echo "QC Report Complete"
EOF

chmod +x checkref_qc.sh
```

### Usage Example

```bash
# Run QC on your analysis
./checkref_qc.sh analysis_results/ input.vcf.gz

# Example output:
# === CheckRef QC Report ===
# Analysis: analysis_results
# Date: Thu Sep 12 10:30:00 2024
# ==========================
# 📊 Summary Statistics:
#   Match rate: 87.5%
#   Error rate: 2.1%
#   ✅ Match rate: EXCELLENT  
#   ✅ Error rate: EXCELLENT
# 📁 Output Files:
#   ✅ chr22.corrected.vcf.gz: 8547392 bytes
#      Variants: 123456
# 🔧 Correction Validation:
#   Corrections applied: 8234
#   ✅ Corrections detected in output
# ==========================
# QC Report Complete
```

## QC Best Practices

### Establish Baseline Metrics

For your specific data types and populations, establish expected ranges:

```bash
# Document expected metrics for your dataset
cat > qc_thresholds.txt << EOF
# QC Thresholds for [Your Dataset]
# Population: [e.g., African]
# Platform: [e.g., Illumina GWAS Array]

Match Rate Thresholds:
- Excellent: >85%
- Good: 70-85%  
- Investigate: <70%

Error Rate Thresholds:
- Excellent: <3%
- Acceptable: 3-8%
- Investigate: >8%

Switch Rate Thresholds:
- Normal: 5-15%
- Investigate: >20%

Target Overlap Thresholds:
- Good: >60%
- Acceptable: 30-60%
- Investigate: <30%
EOF
```

### Regular QC Monitoring

```bash
# Create QC monitoring script for batch processing
cat > monitor_qc.sh << 'EOF'
#!/bin/bash

# Monitor multiple CheckRef runs
for results_dir in results_*/; do
  if [ -d "$results_dir" ]; then
    echo "Checking $results_dir..."
    ./checkref_qc.sh "$results_dir" "input/$(basename $results_dir | sed 's/results_//')"
    echo ""
  fi
done
EOF

chmod +x monitor_qc.sh
```

## Troubleshooting QC Failures

### Low Match Rates

```bash
# Investigate low match rates
results_file="results/*_allele_switch_results.tsv"

# Check switch type distribution
echo "Switch type distribution:"
cut -f7 $results_file | sort | uniq -c

# Examine problematic variants
echo "Sample problematic variants:"
grep -v "MATCH" $results_file | head -10

# Check if issue is population-specific
echo "Consider using population-matched reference panel"
```

### High Error Rates

```bash
# Investigate high error rates  
echo "High error variants:"
grep "OTHER" $results_file | head -20

# Check for systematic issues
grep "OTHER" $results_file | cut -f1,2 | \
  awk '{print $1}' | sort | uniq -c | \
  sort -nr | head -5
```

## Next Steps

After implementing QC:

1. **[Performance](./performance)** - Optimize QC for large datasets
2. **[Batch Processing](./batch-processing)** - Scale QC to multiple samples
3. **[Custom Parameters](./custom-parameters)** - Adjust thresholds based on QC results

## Quick Reference

### Essential QC Commands
```bash
# Pre-analysis validation
bcftools stats input.vcf.gz
file input.vcf.gz && tabix -l input.vcf.gz

# Process monitoring
-with-trace -with-report -with-timeline

# Result validation  
./checkref_qc.sh results/ input.vcf.gz

# Correction verification
bcftools view results/*.corrected.vcf.gz | grep "SWITCHED=1" | wc -l
```

### Quality Thresholds
- **Match rate**: >85% excellent, 70-85% good
- **Error rate**: <5% excellent, <10% acceptable  
- **Target overlap**: >50% good, 30-50% acceptable