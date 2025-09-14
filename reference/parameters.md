# Parameters Reference

Complete reference for all pipeline parameters in CheckRef.

## Core Parameters

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `--targetVcfs` | `string` | Path(s) to target VCF file(s). Comma-separated for multiple files |
| `--referenceDir` | `string` | Directory containing reference legend files |

### Essential Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--outputDir` | `string` | `'results'` | Output directory for results |
| `--fixMethod` | `string` | `'remove'` | Correction method: 'remove' or 'correct' |

## Input/Output Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--legendPattern` | `string` | `"*.legend.gz"` | Pattern to match legend files |
| `--email` | `string` | `null` | Email for completion notification |
| `--email_on_fail` | `string` | `null` | Email on pipeline failure |

## Validation Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--skip_validation` | `boolean` | `false` | Skip validation workflow (correct method only) |
| `--validation_af_threshold` | `float` | `0.05` | Threshold for flagging AF changes |
| `--validation_fold_threshold` | `float` | `2.0` | Threshold for flagging fold changes |

## Resource Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--max_cpus` | `integer` | `16` | Maximum number of CPUs |
| `--max_memory` | `string` | `'128.GB'` | Maximum memory allocation |
| `--max_time` | `string` | `'240.h'` | Maximum time per job |

## Example Usage

### Basic Usage
```bash
nextflow run main.nf \
  --targetVcfs chr22.vcf.gz \
  --referenceDir /path/to/reference/panels/ \
  --outputDir results/
```

### Advanced Usage
```bash
nextflow run main.nf \
  --targetVcfs "chr20.vcf.gz,chr21.vcf.gz,chr22.vcf.gz" \
  --referenceDir /path/to/reference/panels/ \
  --outputDir results/ \
  --fixMethod correct \
  --max_cpus 32 \
  --max_memory '128.GB' \
  --validation_af_threshold 0.10 \
  --email user@example.com
```


