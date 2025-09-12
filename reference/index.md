# CheckRef Reference

Complete reference documentation for CheckRef command-line parameters, profiles, and configuration options.

## Reference Sections

- [**Parameters**](./parameters) - Complete parameter reference
- [**Profiles**](./profiles) - Execution profiles (docker, singularity, etc.)
- [**Configuration**](./configuration) - Advanced configuration options

## Quick Reference

### Essential Parameters
```bash
--targetVcfs     # Path(s) to target VCF file(s) (required)
--referenceDir   # Directory with reference legend files (required)
--outputDir      # Output directory (default: 'results')
--fixMethod      # 'remove' or 'correct' (default: 'remove')
```

### Common Usage Patterns
```bash
# Single file
nextflow run main.nf --targetVcfs file.vcf.gz --referenceDir /ref/panels/

# Multiple files
nextflow run main.nf --targetVcfs "chr*.vcf.gz" --referenceDir /ref/panels/

# Correct method
nextflow run main.nf --targetVcfs file.vcf.gz --referenceDir /ref/panels/ --fixMethod correct
```