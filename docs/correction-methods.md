# Correction Methods Tutorial

Master the two different approaches CheckRef offers for handling allele switches: removal and correction.

## Learning Objectives

By the end of this tutorial, you'll understand:
- When to use each correction method
- How removal and correction methods work
- Trade-offs between data quantity and quality
- How to validate correction success

## The Two Methods

CheckRef offers two approaches to handle allele switches:

| Method | Action | Output File | Best For |
|--------|--------|-------------|----------|
| **Remove** | Delete switched variants | `*.noswitch.vcf.gz` | Quality-first approaches |
| **Correct** | Fix alleles by swapping | `*.corrected.vcf.gz` | Retain maximum variants |

## Method 1: Remove Switched Sites

### How It Works

The remove method identifies variants with allele switches and excludes them from the output VCF:

```bash
nextflow run main.nf \
  --targetVcfs sample.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod remove \
  --outputDir remove_results
```

### Process Steps

1. **Identify switches**: Find variants with STATUS = SWITCH, COMPLEMENT, COMPLEMENT_SWITCH
2. **Create exclusion list**: Generate BED file of positions to remove
3. **Filter VCF**: Use `bcftools view -T ^exclude.bed` to remove sites
4. **Index output**: Create tabix index for the cleaned VCF

### What Gets Removed

```
Original VCF variants: 100,000
├── MATCH: 85,000        → Kept
├── SWITCH: 8,000        → Removed
├── COMPLEMENT: 5,000    → Removed  
├── COMPLEMENT_SWITCH: 1,500 → Removed
└── OTHER: 500           → Removed

Final VCF variants: 85,000
```

### Advantages
- **Highest quality**: Only well-oriented variants remain
- **No correction errors**: No risk of introducing mistakes
- **Downstream compatibility**: Works with all analysis tools
- **Conservative approach**: Minimizes false positives

### Disadvantages
- **Reduced variant count**: May lose 10-20% of variants
- **Power loss**: Fewer variants for association studies
- **Ascertainment bias**: May favor certain variant types

### When to Use Remove Method

**Imputation preparation:**
```bash
# Preparing for imputation server submission
nextflow run main.nf \
  --targetVcfs pre_imputation.vcf.gz \
  --referenceDir /imputation/ref/panels/ \
  --fixMethod remove
```

**Quality-critical analyses:**
```bash
# Fine-mapping or functional studies
nextflow run main.nf \
  --targetVcfs gwas_hits.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod remove
```

## Method 2: Correct Switched Sites

### How It Works

The correct method identifies switched variants and fixes them by swapping REF and ALT alleles:

```bash
nextflow run main.nf \
  --targetVcfs sample.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod correct \
  --outputDir correct_results
```

### Process Steps

1. **Identify switches**: Find correctable switch types
2. **Generate correction script**: Create Python script for allele swapping
3. **Apply corrections**: Swap REF/ALT for switched variants
4. **Add flags**: Mark corrected variants with `SWITCHED=1` in INFO field
5. **Sort and index**: Ensure proper VCF format

### What Gets Corrected

```
Original switches: 15,000
├── SWITCH: 8,000           → Corrected (REF↔ALT swap)
├── COMPLEMENT: 5,000       → Corrected (strand flip + swap)
├── COMPLEMENT_SWITCH: 1,500 → Corrected (complex correction)
└── OTHER: 500              → Removed (uncorrectable)

Final result: 99,500 variants (14,500 corrected)
```

### Correction Examples

**Simple switch (SWITCH):**
```
Before: REF=G ALT=A (vs reference REF=A ALT=G)  
After:  REF=A ALT=G + INFO=SWITCHED=1
```

**Complement switch (COMPLEMENT):**
```  
Before: REF=C ALT=T (vs reference REF=G ALT=A)
After:  REF=G ALT=A + INFO=SWITCHED=1
```

**Complex switch (COMPLEMENT_SWITCH):**
```
Before: REF=T ALT=C (vs reference REF=A ALT=G)
After:  REF=A ALT=G + INFO=SWITCHED=1  
```

### Advantages
- **Maximum variant retention**: Keeps ~95% of variants
- **Higher power**: More variants for association studies  
- **Transparent process**: Corrections are flagged and trackable
- **Reversible**: Can identify corrected variants later

### Disadvantages
- **Complexity**: More complex process with potential errors
- **Validation needed**: Requires verification of corrections
- **Tool compatibility**: Some tools may not handle SWITCHED flag
- **Population genetics**: May affect allele frequency estimates

### When to Use Correct Method

**Association studies:**
```bash
# GWAS where variant count matters
nextflow run main.nf \
  --targetVcfs gwas_cohort.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod correct
```

**Rare variant analysis:**
```bash
# Exome sequencing where every variant counts
nextflow run main.nf \
  --targetVcfs exome_rare_variants.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod correct
```

## Comparing Methods Side-by-Side

Let's compare both methods on the same dataset:

### Remove Method Results
```bash
nextflow run main.nf \
  --targetVcfs chr22.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod remove \
  --outputDir comparison_remove

# Results:
# Original variants: 125,000
# Final variants: 108,000 (86.4% retention)
# Removed: 17,000 variants
```

### Correct Method Results  
```bash
nextflow run main.nf \
  --targetVcfs chr22.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod correct \
  --outputDir comparison_correct

# Results:
# Original variants: 125,000
# Final variants: 123,500 (98.8% retention)  
# Corrected: 16,500 variants
# Removed: 1,500 uncorrectable variants
```

### Quality Comparison

```bash
# Check correction success rate
bcftools view comparison_correct/chr22.corrected.vcf.gz | \
  grep "SWITCHED=1" | wc -l
# Output: 16,500 corrected variants

# Verify against expected corrections
grep -E "SWITCH|COMPLEMENT" comparison_remove/chr22_allele_switch_results.tsv | \
  wc -l
# Output: 16,500 switch events detected

# Success rate: 100% (all detectable switches were corrected)
```

## Validation Strategies

### Method 1: Cross-Reference Validation

Compare your results with known datasets:

```bash
# Check corrected alleles against dbSNP
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' corrected.vcf.gz | \
  head -10 > corrected_variants.txt

# Validate against reference database
# (implementation depends on available resources)
```

### Method 2: Frequency Consistency

Check if corrected allele frequencies make sense:

```bash
# Calculate allele frequencies before and after
bcftools +fill-tags original.vcf.gz -- -t AF | \
  bcftools query -f '%POS\t%AF\n' > original_af.txt

bcftools +fill-tags corrected.vcf.gz -- -t AF | \
  bcftools query -f '%POS\t%AF\n' > corrected_af.txt

# Compare frequencies for corrected sites
```

### Method 3: Downstream Analysis Check

Test both methods in your analysis pipeline:

```bash
# Run association test with both methods
plink --vcf remove_results/chr22.noswitch.vcf.gz --assoc
plink --vcf correct_results/chr22.corrected.vcf.gz --assoc

# Compare results for consistency
```

## Decision Framework

Use this framework to choose the best method for your analysis:

### Choose REMOVE when:
- Data quality is paramount
- Preparing for imputation
- Limited computational resources
- Conservative approach preferred
- Downstream tools are sensitive to orientation

### Choose CORRECT when:
- Maximum variant count needed
- Power is critical (association studies)
- Rare variant analysis
- You can validate corrections
- Analysis tools handle flagged variants

### Hybrid Approach

You can also use both methods strategically:

```bash
# Step 1: Run correct method for maximum information
nextflow run main.nf \
  --targetVcfs sample.vcf.gz \
  --referenceDir /ref/panels/ \
  --fixMethod correct

# Step 2: Filter out high-risk corrections for critical analyses
bcftools view -e 'INFO/SWITCHED=1 && QUAL<30' \
  corrected.vcf.gz > high_confidence.vcf.gz
```

## Advanced Correction Scenarios

### Population-Specific Corrections

Different populations may have different switch patterns:

```bash
# European population (typically fewer strand issues)
nextflow run main.nf \
  --targetVcfs european_samples.vcf.gz \
  --referenceDir /ref/1000GP_EUR/ \
  --fixMethod correct

# African population (may have more complex patterns)
nextflow run main.nf \
  --targetVcfs african_samples.vcf.gz \
  --referenceDir /ref/H3Africa/ \
  --fixMethod remove  # More conservative
```

### Platform-Specific Approaches

Different genotyping platforms may need different strategies:

```bash
# SNP array data (well-characterized sites)
--fixMethod correct

# Whole genome sequencing (novel variants)
--fixMethod remove

# Targeted sequencing (functional sites)  
--fixMethod remove
```

## Monitoring Correction Quality

### Real-Time Validation

CheckRef provides built-in validation:

```bash
# Check the correction verification output
cat results/verification/chr22_verification_results.txt

# Expected output:
# ✓ Corrections verified for chr22
# Sites corrected: 16,500
# Sites failed: 0
# Verification success rate: 100%
```

### Manual Quality Checks

```bash
# Check a few corrected sites manually
bcftools view corrected.vcf.gz 22:16050075-16050080 | \
  grep "SWITCHED=1"

# Compare with original
bcftools view original.vcf.gz 22:16050075-16050080

# Check reference panel
zcat /ref/panels/chr22.legend.gz | grep "16050075"
```

## Next Steps

After mastering correction methods:

1. **[Quality Control](./quality-control)** - Validate your corrections
2. **[Performance](./performance)** - Optimize correction speed
3. **[Custom Parameters](./custom-parameters)** - Fine-tune correction behavior
4. **[Batch Processing](./batch-processing)** - Apply methods at scale

## Quick Reference

### Method Comparison
| Aspect | Remove | Correct |
|--------|--------|---------|
| **Variant retention** | ~85% | ~98% |
| **Quality** | Highest | Good |
| **Complexity** | Simple | Moderate |
| **Validation** | Not needed | Recommended |
| **Downstream compatibility** | Universal | Good |

### Key Commands
```bash
# Remove method
--fixMethod remove

# Correct method  
--fixMethod correct

# Check corrected variants
bcftools view file.vcf.gz | grep "SWITCHED=1"

# Count corrections
bcftools view -H file.vcf.gz | grep "SWITCHED=1" | wc -l
```