# Quick Start Tutorial

Get up and running with CheckRef in 10 minutes. This tutorial will walk you through your first successful analysis.

## What You'll Accomplish

- Run CheckRef on a sample VCF file
- Understand the basic output
- Choose between remove and correct methods
- Know what to do next

## Prerequisites

- Nextflow installed
- Docker or Singularity available
- A VCF file to analyze
- 10 minutes

## Step 1: Prepare Your Files (2 minutes)

Check that you have the required files:

```bash
# Check your VCF file
ls -la your_file.vcf.gz your_file.vcf.gz.tbi

# Verify it's readable
bcftools view -h your_file.vcf.gz | head -3
```

Your file should be:
- Compressed with bgzip (`.vcf.gz`)
- Indexed (`.vcf.gz.tbi` file present)
- Contain SNP variants

## Step 2: Run CheckRef (5 minutes)

Run your first CheckRef analysis:

```bash
nextflow run main.nf \
  --targetVcfs your_file.vcf.gz \
  --referenceDir /path/to/reference/panels/ \
  --outputDir my_first_analysis \
  -profile singularity
```

**What happens:**
- CheckRef compares your VCF against reference panels
- Identifies variants with allele switches
- Creates cleaned output files
- Usually completes in 2-5 minutes for typical chromosome files

## Step 3: Check Your Results (2 minutes)

Look at what CheckRef created:

```bash
# Check output directory
ls -la my_first_analysis/

# Quick summary
cat my_first_analysis/*_summary.txt
```

**Key numbers to look for:**
```
Matched variants: 4250 (85.00%)     # Good: >80%
Switched alleles: 350 (7.00%)       # Normal: 5-15% 
Other inconsistencies: 25 (0.50%)   # Good: <5%
```

## Step 4: Understand Your Method (1 minute)

CheckRef used the **remove method** by default:
- **What it did**: Removed variants with allele switches
- **Output file**: `*.noswitch.vcf.gz`
- **Result**: Higher quality, fewer variants

**Alternative - correct method:**
```bash
# To fix switches instead of removing them
nextflow run main.nf \
  --targetVcfs your_file.vcf.gz \
  --referenceDir /path/to/reference/panels/ \
  --fixMethod correct \
  --outputDir corrected_analysis
```

## Success Indicators

Your analysis succeeded if you see:
- **Match rate >70%**: Good alignment with reference
- **Error rate <10%**: Low problematic variants  
- **Output VCF created**: Fixed file generated
- **No error messages**: Pipeline completed successfully

## What Each Method Gives You

| Method | Output File | Variants Kept | Quality |
|--------|-------------|---------------|---------|
| **remove** | `*.noswitch.vcf.gz` | ~85-90% | Highest |
| **correct** | `*.corrected.vcf.gz` | ~95-98% | High |

## Quick Validation

Verify your results worked:

```bash
# Count variants before and after
echo "Original: $(bcftools view -H your_file.vcf.gz | wc -l)"
echo "Final: $(bcftools view -H my_first_analysis/*.vcf.gz | wc -l)"

# Should show reasonable retention (>80% for remove, >95% for correct)
```

## Common First-Run Issues

**"No legend files found":**
```bash
# Check reference directory
ls /path/to/reference/panels/
# Use --legendPattern if files have different names
```

**"Out of memory":**
```bash
# Add more memory
--max_memory '16.GB'
```

**"Container not found":**
```bash
# Try different profile
-profile docker
# Or use standard (no containers)
-profile standard
```

## Next Steps

Congratulations! You've completed your first CheckRef analysis. Here's what to explore next:

**If your results look good:**
- [Multi-File Tutorial](./multi-file-tutorial) - Process multiple chromosomes
- [Method Selection Tutorial](./method-selection) - Choose remove vs correct

**If you had issues:**
- [Understanding Results](/docs/understanding-results) - Interpret your output
- [Troubleshooting Guide](/docs/troubleshooting) - Solve common problems

**For more details:**
- [Documentation](/docs/) - Comprehensive reference material
- [Examples](/examples/) - More complex scenarios

## Quick Reference

### Essential Command
```bash
nextflow run main.nf \
  --targetVcfs YOUR_FILE.vcf.gz \
  --referenceDir /path/to/reference/ \
  --outputDir results \
  [--fixMethod correct] \
  -profile [docker|singularity|standard]
```

### Key Output Files
- `*_summary.txt` - Quality metrics
- `*.noswitch.vcf.gz` - Cleaned VCF (remove method)
- `*.corrected.vcf.gz` - Fixed VCF (correct method)

### Success Thresholds
- Match rate: >70% acceptable, >85% excellent
- Error rate: <10% acceptable, <5% excellent