# Troubleshooting Common Issues

Solve the most common problems encountered when using CheckRef.

## Learning Objectives

By the end of this tutorial, you'll be able to:
- Diagnose common CheckRef errors
- Fix input data issues
- Resolve resource and environment problems
- Optimize problematic runs

## Quick Diagnostic Checklist

Before diving into specific issues, run through this checklist:

```bash
# 1. Check file accessibility
ls -la /path/to/target.vcf.gz
ls -la /path/to/reference/panels/

# 2. Verify file formats
file /path/to/target.vcf.gz
bcftools view -h /path/to/target.vcf.gz | head -5

# 3. Check Nextflow status
nextflow info
nextflow version

# 4. Verify container availability (if using)
singularity --version
docker --version
```

## Input File Issues

### Issue 1: File Not Found

**Error Message:**
```
ERROR ~ Input file does not exist: /path/to/file.vcf.gz
```

**Diagnosis:**
```bash
# Check if file exists
ls -la /path/to/file.vcf.gz

# Check permissions
ls -la $(dirname /path/to/file.vcf.gz)

# Verify working directory
pwd
```

**Solutions:**
```bash
# Use absolute paths
--targetVcfs /full/path/to/file.vcf.gz

# Check current directory
--targetVcfs ./relative/path/to/file.vcf.gz

# Fix permissions if needed
chmod 644 /path/to/file.vcf.gz
```

### Issue 2: VCF Format Problems

**Error Message:**
```
ERROR ~ Invalid VCF format or corrupted file
```

**Diagnosis:**
```bash
# Check file integrity
bcftools view -h target.vcf.gz | head -10

# Validate VCF format
bcftools view -H target.vcf.gz | head -5

# Check compression
file target.vcf.gz
```

**Solutions:**
```bash
# Recompress if needed
bcftools view target.vcf | bgzip > target.vcf.gz
tabix -p vcf target.vcf.gz

# Fix VCF header issues
bcftools view target.vcf.gz | bcftools reheader -s samples.txt | \
  bgzip > fixed.vcf.gz

# Validate and fix
bcftools norm -f reference.fa -c w input.vcf.gz | \
  bgzip > normalized.vcf.gz
```

### Issue 3: Empty VCF Files

**Error Message:**
```
ERROR ~ No variants found in target VCF
```

**Diagnosis:**
```bash
# Count variants
bcftools view -H target.vcf.gz | wc -l

# Check for specific chromosome
bcftools view -H target.vcf.gz | cut -f1 | sort -u

# Look for filtering issues
bcftools view -f PASS target.vcf.gz | bcftools view -H | wc -l
```

**Solutions:**
```bash
# Remove strict filters
bcftools view -i 'FILTER="."' target.vcf.gz > unfiltered.vcf.gz

# Check SNPs only
bcftools view -v snps target.vcf.gz > snps_only.vcf.gz

# Verify chromosome naming
bcftools annotate --rename-chrs chr_name_mapping.txt target.vcf.gz
```

## Reference Panel Issues

### Issue 4: No Reference Files Found

**Error Message:**
```
ERROR ~ No legend files found matching pattern '*.legend.gz' in directory
```

**Diagnosis:**
```bash
# List reference directory contents
ls -la /path/to/reference/dir/

# Check for different file extensions
ls -la /path/to/reference/dir/ | grep -E "\.(legend|leg)"

# Test pattern matching
ls /path/to/reference/dir/*.legend.gz
```

**Solutions:**
```bash
# Adjust pattern to match your files
--legendPattern "*.legend"           # Uncompressed
--legendPattern "*chr22*.legend.gz"  # Specific chromosome
--legendPattern "*.leg.gz"           # Different extension

# Use specific directory structure
--referenceDir /path/to/chr22/
--legendPattern "reference.legend.gz"
```

### Issue 5: Chromosome Mismatch

**Error Message:**
```
ERROR ~ Cannot match chromosome 'chr22' with any reference file
```

**Diagnosis:**
```bash
# Check target chromosome naming
bcftools view -H target.vcf.gz | cut -f1 | sort -u

# Check reference file names
ls /ref/panels/ | grep -E "(22|chr22)"

# Look inside reference files
zcat /ref/panels/reference.legend.gz | head -5
```

**Solutions:**
```bash
# Option 1: Rename chromosomes in VCF
bcftools annotate --rename-chrs chr_rename.txt target.vcf.gz

# chr_rename.txt content:
# chr1 1
# chr22 22
# chrX X

# Option 2: Use pattern that matches reference files
--legendPattern "*22*.legend.gz"     # If files are named with numbers
--legendPattern "*chr22*.legend.gz"  # If files include 'chr'

# Option 3: Organize reference files by chromosome
mkdir -p /ref/panels/chr22/
cp chr22.legend.gz /ref/panels/chr22/
--referenceDir /ref/panels/chr22/
```

## Memory and Resource Issues

### Issue 6: Out of Memory

**Error Message:**
```
ERROR ~ Process exceeded memory limit (4 GB)
```

**Diagnosis:**
```bash
# Check system memory
free -h

# Check file sizes
ls -lh *.vcf.gz

# Monitor memory usage during run
top -p $(pgrep -f nextflow)
```

**Solutions:**
```bash
# Increase memory allocation
--max_memory '16.GB'

# Process-specific memory
nextflow run main.nf \
  --targetVcfs large_file.vcf.gz \
  --referenceDir /ref/panels/ \
  --max_memory '32.GB'

# For very large files, consider splitting:
bcftools view -r 22:1-30000000 input.vcf.gz > part1.vcf.gz
bcftools view -r 22:30000001-50000000 input.vcf.gz > part2.vcf.gz
```

### Issue 7: CPU Resource Problems

**Error Message:**
```
ERROR ~ Process failed due to timeout
```

**Diagnosis:**
```bash
# Check available CPUs
nproc

# Monitor CPU usage
htop

# Check Nextflow resource allocation
nextflow log last -f status,cpus,memory,time
```

**Solutions:**
```bash
# Increase CPU allocation
--max_cpus 8

# Increase time limits
--max_time '48.h'

# Use HPC profile for cluster environments
nextflow run main.nf \
  --targetVcfs input.vcf.gz \
  --referenceDir /ref/panels/ \
  -profile slurm \
  --max_cpus 16 \
  --max_memory '64.GB'
```

## Container and Environment Issues

### Issue 8: Container Not Found

**Error Message:**
```
ERROR ~ Container 'mamana/vcf-processing:latest' not found
```

**Diagnosis:**
```bash
# Check container availability
docker images | grep vcf-processing
singularity cache list

# Test container manually
docker run mamana/vcf-processing:latest bcftools --version
```

**Solutions:**
```bash
# Pull container manually
docker pull mamana/vcf-processing:latest

# For Singularity
singularity pull docker://mamana/vcf-processing:latest

# Use different profile if containers unavailable
-profile standard  # Uses local tools instead of containers
```

### Issue 9: Tool Not Found

**Error Message:**
```
ERROR ~ bcftools: command not found
```

**Diagnosis:**
```bash
# Check if tools are installed
which bcftools
which python3
which bgzip

# Check PATH
echo $PATH
```

**Solutions:**
```bash
# Install missing tools
# Ubuntu/Debian:
sudo apt-get install bcftools tabix python3

# CentOS/RHEL:
sudo yum install bcftools python3

# Using conda:
conda install -c bioconda bcftools
conda install -c conda-forge python=3.9

# Use container profile
-profile docker    # or -profile singularity
```

## Data Quality Issues

### Issue 10: Poor Match Rates

**Symptoms:**
```
Results Summary:
Matched variants: 450 (45.00%)  # Low match rate
Other inconsistencies: 300 (30.00%)  # High error rate
```

**Diagnosis:**
```bash
# Check data quality metrics
bcftools stats target.vcf.gz > target_stats.txt
cat target_stats.txt

# Examine problematic variants
grep "OTHER" results/*_allele_switch_results.tsv | head -10

# Check reference panel version
ls -la /ref/panels/
cat /ref/panels/README  # If available
```

**Solutions:**
```bash
# Apply quality filters
bcftools view -i 'QUAL>=20 && INFO/DP>=10' target.vcf.gz > filtered.vcf.gz

# Try different reference panel
--referenceDir /alternative/ref/panels/

# Check population match
# Use population-specific reference panels for better results

# Filter for common variants only
bcftools view -i 'INFO/AF>0.01' target.vcf.gz > common_variants.vcf.gz
```

### Issue 11: Inconsistent Results Across Chromosomes

**Symptoms:**
```
chr20: 90% match rate
chr21: 45% match rate  # Inconsistent
chr22: 88% match rate
```

**Diagnosis:**
```bash
# Check individual chromosome files
bcftools stats chr21.vcf.gz
bcftools view -H chr21.vcf.gz | wc -l

# Compare chromosome naming
bcftools view -H chr21.vcf.gz | cut -f1 | sort -u
ls /ref/panels/ | grep chr21
```

**Solutions:**
```bash
# Reprocess problematic chromosome
nextflow run main.nf \
  --targetVcfs chr21.vcf.gz \
  --referenceDir /ref/panels/ \
  --legendPattern "*chr21*.legend.gz" \
  --outputDir chr21_rerun

# Check chromosome-specific issues
bcftools view chr21.vcf.gz | grep "^##" | grep -i reference
```

## Pipeline Execution Issues

### Issue 12: Pipeline Hangs

**Symptoms:**
- Pipeline stops progressing
- No error messages
- Processes stuck in RUNNING state

**Diagnosis:**
```bash
# Check pipeline status
nextflow log last -f name,status,exit,duration

# Monitor system resources
htop
df -h  # Check disk space

# Check work directory
ls -la work/
du -sh work/
```

**Solutions:**
```bash
# Resume with debugging
nextflow run main.nf \
  --targetVcfs input.vcf.gz \
  --referenceDir /ref/panels/ \
  -resume \
  -with-trace \
  -with-report \
  -with-timeline

# Clean work directory if disk is full
nextflow clean -f
rm -rf work/

# Increase resource limits
--max_time '72.h'
--max_memory '128.GB'
```

### Issue 13: Resume Fails

**Error Message:**
```
ERROR ~ Cannot resume pipeline - work directory corrupted
```

**Diagnosis:**
```bash
# Check work directory integrity
ls -la work/
find work/ -name "*.exitcode" | head -10

# Check Nextflow cache
ls -la .nextflow/
```

**Solutions:**
```bash
# Clean and restart
nextflow clean -f
rm -rf .nextflow/
rm -rf work/

# Start fresh run
nextflow run main.nf \
  --targetVcfs input.vcf.gz \
  --referenceDir /ref/panels/ \
  --outputDir fresh_run

# Use different work directory
nextflow run main.nf \
  --targetVcfs input.vcf.gz \
  --referenceDir /ref/panels/ \
  -work-dir /tmp/nextflow_work
```

## Performance Issues

### Issue 14: Slow Processing

**Symptoms:**
- Long processing times
- Single-threaded execution
- Low resource utilization

**Diagnosis:**
```bash
# Monitor resource usage
htop
iotop  # Check I/O usage

# Check parallelization
nextflow log last -f name,status,cpus,memory

# Profile bottlenecks
time nextflow run main.nf --targetVcfs small_test.vcf.gz --referenceDir /ref/
```

**Solutions:**
```bash
# Increase parallelization
--max_cpus $(nproc)

# Use faster storage for work directory
-work-dir /fast/ssd/work

# Optimize I/O
export NXF_OPTS='-Xms2g -Xmx8g'

# Use local reference copies
cp -r /network/ref/panels/ /local/tmp/ref/
--referenceDir /local/tmp/ref/
```

## Advanced Troubleshooting

### Debug Mode

Enable detailed logging:

```bash
nextflow run main.nf \
  --targetVcfs debug.vcf.gz \
  --referenceDir /ref/panels/ \
  -with-trace trace.txt \
  -with-report report.html \
  -with-timeline timeline.html \
  -with-dag flowchart.html
```

### Manual Process Testing

Test individual steps manually:

```bash
# Test allele switch detection manually
python3 bin/check_allele_switch.py \
  input.vcf.gz \
  reference.legend.gz \
  output_prefix

# Test VCF correction manually
python3 bin/correct_vcf.py \
  input.vcf.gz \
  switch_results.tsv \
  corrected_output.vcf.gz
```

### Log Analysis

Examine detailed logs:

```bash
# Find failed processes
grep -r "ERROR" work/

# Check specific process logs
cat work/*/script.log

# Analyze resource usage
awk -F'\t' '$4=="COMPLETED" {print $1, $8, $9}' trace.txt | sort -k3 -n
```

## Getting Help

### When to Seek Additional Support

Contact support when:
- Multiple troubleshooting steps fail
- Error messages are unclear
- Performance issues persist despite optimization
- Results seem inconsistent with expectations

### Information to Provide

When reporting issues, include:

1. **Command used:**
   ```bash
   nextflow run main.nf --targetVcfs ... --referenceDir ...
   ```

2. **Error messages:**
   ```
   Full error output from terminal
   ```

3. **Environment details:**
   ```bash
   nextflow version
   uname -a
   free -h
   ```

4. **File information:**
   ```bash
   ls -lh *.vcf.gz
   bcftools view -h input.vcf.gz | head -5
   ```

5. **Log files:**
   - `.nextflow.log`
   - `work/*/script.log` (for failed processes)

## Quick Reference

### Emergency Commands
```bash
# Clean everything and restart
nextflow clean -f && rm -rf .nextflow/ work/

# Quick test with small dataset
--targetVcfs small_test.vcf.gz --max_memory 4.GB --max_cpus 2

# Resume with debugging
-resume -with-trace -with-report

# Force container refresh  
docker pull mamana/vcf-processing:latest
```

### Common Parameter Fixes
```bash
# Memory issues
--max_memory '32.GB'

# Time issues
--max_time '48.h'

# File access issues
--targetVcfs $(realpath input.vcf.gz)

# Reference panel issues
--legendPattern "*.legend" --referenceDir $(realpath /ref/panels/)
```