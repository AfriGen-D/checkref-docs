# Understanding CheckRef Results

Learn how to interpret and analyze CheckRef output files to make informed decisions about your data.

## Learning Objectives

By the end of this tutorial, you'll understand:
- The structure of CheckRef output files
- How to interpret allele switch statistics
- What different switch types mean
- How to validate your results

## CheckRef Output Structure

When CheckRef completes, it creates several important files:

```
results/
├── chr22_allele_switch_results.tsv    # Detailed variant-by-variant analysis
├── chr22_allele_switch_summary.txt    # Statistical summary
├── chr22.noswitch.vcf.gz             # Fixed VCF (remove method)
├── chr22.corrected.vcf.gz            # Fixed VCF (correct method)
└── reports/                          # Nextflow execution reports
```

## Detailed Results File

### File Format
The `*_allele_switch_results.tsv` file contains detailed information for each variant:

```
CHROM   POS      TARGET_REF  TARGET_ALT  REF_REF  REF_ALT  STATUS
22      16050075 A           G           A        G        MATCH
22      16050115 G           A           A        G        SWITCH
22      16050213 C           T           G        A        COMPLEMENT
22      16050298 T           C           A        G        COMPLEMENT_SWITCH
22      16050350 A           C           T        G        OTHER
```

### Column Descriptions

| Column | Description |
|--------|-------------|
| `CHROM` | Chromosome identifier |
| `POS` | Genomic position (1-based) |
| `TARGET_REF` | Reference allele in your VCF |
| `TARGET_ALT` | Alternate allele in your VCF |
| `REF_REF` | Reference allele in the reference panel |
| `REF_ALT` | Alternate allele in the reference panel |
| `STATUS` | Classification of the allele relationship |

### Status Classifications

#### MATCH
- **Meaning**: Alleles match perfectly between target and reference
- **Example**: Target A/G matches Reference A/G
- **Action**: No correction needed
- **Interpretation**: High-quality, well-oriented variant

#### SWITCH
- **Meaning**: REF and ALT alleles are flipped
- **Example**: Target G/A vs Reference A/G
- **Action**: Can be corrected by swapping REF/ALT
- **Interpretation**: Common issue, easily fixable

#### COMPLEMENT
- **Meaning**: Alleles are strand complements (A↔T, C↔G)
- **Example**: Target C/T vs Reference G/A
- **Action**: Requires strand flip correction
- **Interpretation**: Different strand representation

#### COMPLEMENT_SWITCH
- **Meaning**: Both complementary and switched
- **Example**: Target T/C vs Reference A/G
- **Action**: Requires both strand flip and allele swap
- **Interpretation**: Complex but correctable

#### OTHER
- **Meaning**: Alleles don't match any expected pattern
- **Example**: Target A/C vs Reference T/G
- **Action**: Usually removed from analysis
- **Interpretation**: Potential genotyping error or complex variant

## Summary Statistics File

### File Format
The `*_allele_switch_summary.txt` provides overview statistics:

```
Results Summary:
Total variants in target: 1321
Total variants in reference: 1741597
Common variants: 24
Matched variants: 4 (16.67%)
Switched alleles: 4 (16.67%)
Complementary strand issues: 8 (33.33%)
Complement + switch issues: 6 (25.00%)
Other inconsistencies: 2 (8.33%)
Target overlap: 1.82%
Reference overlap: 0.00%
```

### Key Metrics

#### Overlap Statistics
- **Total variants in target**: Total SNPs in your VCF file
- **Total variants in reference**: Total SNPs in reference panel
- **Common variants**: Variants present in both datasets
- **Target overlap**: Percentage of your variants found in reference
- **Reference overlap**: Percentage of reference variants in your data

#### Quality Indicators
- **Matched variants**: Well-oriented, high-quality variants
- **Switch rates**: Proportion of variants with allele orientation issues

### Interpreting Summary Statistics

#### Good Quality Indicators
```
Matched variants: 4500 (85.00%)    # High match rate
Switched alleles: 400 (7.50%)      # Moderate switch rate  
Other inconsistencies: 100 (2.00%) # Low error rate
Target overlap: 78.50%             # Good reference coverage
```

#### Concerning Patterns
```
Matched variants: 1200 (45.00%)    # Low match rate
Other inconsistencies: 800 (30.00%) # High error rate
Target overlap: 15.20%             # Poor reference coverage
```

## Analyzing Your Results

### Step 1: Check Overall Quality

Look for these patterns in your summary:

**Excellent Quality (>80% matches)**
- Proceed with confidence
- Most variants are well-oriented
- Minimal data quality issues

**Good Quality (60-80% matches)**
- Acceptable for most analyses
- Consider the switch correction method
- Monitor downstream results

**Poor Quality (<60% matches)**
- Investigate data preparation issues
- Check chromosome naming conventions
- Verify reference panel compatibility

### Step 2: Examine Switch Patterns

**Normal Switch Rates (5-15%)**
- Expected variation between datasets
- Proceed with chosen correction method

**High Switch Rates (>25%)**
- May indicate systematic issues
- Check strand consistency
- Verify reference panel version

**Very High "OTHER" Category (>10%)**
- Possible data quality issues
- May need additional QC steps
- Consider different reference panel

### Step 3: Validate Corrections

#### For Remove Method
Check how many variants remain:
```bash
# Count original variants
bcftools view -H original.vcf.gz | wc -l

# Count remaining variants
bcftools view -H results/chr22.noswitch.vcf.gz | wc -l

# Calculate retention rate
echo "Retention rate: $((remaining * 100 / original))%"
```

#### For Correct Method
Verify corrections were applied:
```bash
# Check corrected variants
bcftools view results/chr22.corrected.vcf.gz | grep "SWITCHED=1" | wc -l

# Compare with expected corrections
grep "SWITCH\|COMPLEMENT" results/chr22_allele_switch_results.tsv | wc -l
```

## Real-World Examples

### Example 1: High-Quality Dataset
```
Results Summary:
Total variants in target: 125000
Common variants: 95000
Matched variants: 85500 (90.00%)
Switched alleles: 6650 (7.00%)  
Other inconsistencies: 2850 (3.00%)
Target overlap: 76.00%
```

**Interpretation**: Excellent quality data with high match rates and good reference coverage. Safe to proceed with either correction method.

### Example 2: Reference Panel Mismatch
```
Results Summary:
Total variants in target: 89000
Common variants: 12000
Matched variants: 5400 (45.00%)
Other inconsistencies: 4800 (40.00%)
Target overlap: 13.48%
```

**Interpretation**: Poor reference panel match. Consider using a different reference panel or checking population ancestry.

### Example 3: Strand Issues
```
Results Summary:
Common variants: 78000
Matched variants: 23400 (30.00%)
Complementary strand issues: 39000 (50.00%)
Complement + switch issues: 11700 (15.00%)
```

**Interpretation**: Systematic strand differences. The correct method should handle these well.

## Troubleshooting Results

### Low Match Rates

**Possible Causes:**
1. Wrong reference panel population
2. Different genome builds (hg19 vs hg38)
3. Chromosome naming differences (chr22 vs 22)
4. Poor quality target data

**Solutions:**
```bash
# Check chromosome naming
bcftools view -H target.vcf.gz | cut -f1 | sort -u

# Check reference panel chromosomes
ls reference_panels/ | grep chr22

# Verify genome build
bcftools view -h target.vcf.gz | grep "##reference"
```

### High Error Rates

**Possible Causes:**
1. Genotyping errors in target data
2. Reference panel quality issues
3. Population structure differences

**Solutions:**
- Apply stricter quality filters to target VCF
- Use population-matched reference panels
- Consider multiple reference panels

### Unexpected Switch Patterns

**Investigation Steps:**
```bash
# Check specific problematic variants
grep "OTHER" results/chr22_allele_switch_results.tsv | head -10

# Look at genomic context
bcftools view target.vcf.gz 22:16050075-16050080

# Compare with reference
zcat reference_panels/chr22.legend.gz | grep "16050075"
```

## Next Steps

Once you understand your results:

1. **Good Results**: Proceed to [Multi-File Processing](./multi-file)
2. **Quality Issues**: Review [Data Preparation](./data-preparation)
3. **Performance Concerns**: Check [Performance Optimization](./performance)
4. **Need Corrections**: Learn about [Correction Methods](./correction-methods)

## Quick Reference

### Key Quality Thresholds
- **Match Rate**: >70% good, >85% excellent
- **Error Rate**: <10% good, <5% excellent  
- **Overlap Rate**: >50% adequate, >70% good

### Essential Commands
```bash
# View summary quickly
cat results/*_summary.txt

# Count each status type
cut -f7 results/*_results.tsv | sort | uniq -c

# Check correction success
bcftools view results/*.corrected.vcf.gz | grep "SWITCHED=1" | wc -l
```