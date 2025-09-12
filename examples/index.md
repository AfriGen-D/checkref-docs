# Examples

This section provides practical examples for running {{ PROJECT_NAME }} in different scenarios.

## Basic Examples

### Example 1: Minimal Configuration

The simplest way to run the pipeline:

```bash
nextflow run main.nf \
  --input samples.csv \
  --outdir results/
```

**Sample sheet format:**
```csv
sample,fastq_1,fastq_2
sample1,data/sample1_R1.fastq.gz,data/sample1_R2.fastq.gz
sample2,data/sample2_R1.fastq.gz,data/sample2_R2.fastq.gz
```

### Example 2: Custom Reference Genome

```bash
nextflow run main.nf \
  --input samples.csv \
  --outdir results/ \
  --fasta /path/to/custom_genome.fa \
  --bwa_index /path/to/custom_genome_bwa/
```

## Analysis-Specific Examples

### {{ EXAMPLE_ANALYSIS_1 }}

```bash
nextflow run main.nf \
  --input samples.csv \
  --outdir results/ \
  --genome GRCh38 \
  --{{ PARAM1 }} {{ VALUE1 }} \
  --{{ PARAM2 }} {{ VALUE2 }}
```

**Parameters file (params_{{ EXAMPLE_1 }}.yml):**
```yaml
input: 'samples.csv'
outdir: 'results_{{ EXAMPLE_1 }}/'
genome: 'GRCh38'
{{ PARAM1 }}: {{ VALUE1 }}
{{ PARAM2 }}: {{ VALUE2 }}
skip_multiqc: false
```

### {{ EXAMPLE_ANALYSIS_2 }}

```bash
nextflow run main.nf \
  --input samples.csv \
  --outdir results/ \
  --genome GRCh38 \
  --{{ PARAM3 }} {{ VALUE3 }} \
  --{{ PARAM4 }} {{ VALUE4 }}
```

## Environment-Specific Examples

### Local Execution with Docker

```bash
nextflow run main.nf \
  --input samples.csv \
  --outdir results/ \
  -profile docker
```

### HPC Cluster (SLURM)

```bash
nextflow run main.nf \
  --input samples.csv \
  --outdir results/ \
  -profile slurm \
  --max_cpus 128 \
  --max_memory '512.GB' \
  --max_time '48.h'
```

### Cloud Execution (AWS)

```bash
nextflow run main.nf \
  --input s3://bucket/samples.csv \
  --outdir s3://bucket/results/ \
  -profile awsbatch \
  -work-dir s3://bucket/work/
```

## Advanced Examples

### Large-Scale Analysis

For processing hundreds of samples:

```bash
nextflow run main.nf \
  --input large_cohort_samples.csv \
  --outdir results_large/ \
  --genome GRCh38 \
  --max_cpus 64 \
  --max_memory '256.GB' \
  -profile slurm \
  -resume
```

**Configuration file (large_cohort.config):**
```groovy
params {
    max_cpus = 64
    max_memory = '256.GB'
    max_time = '168.h'
}

process {
    withLabel: 'process_high_memory' {
        memory = '128.GB'
        time = '24.h'
    }
    
    withName: '{{ HIGH_MEMORY_PROCESS }}' {
        cpus = 16
        memory = '64.GB'
    }
}
```

### Custom Quality Filters

```bash
nextflow run main.nf \
  --input samples.csv \
  --outdir results_filtered/ \
  --genome GRCh38 \
  --min_mapping_quality 30 \
  --min_base_quality 25 \
  --remove_duplicates true
```

### Skip Certain Steps

```bash
nextflow run main.nf \
  --input samples.csv \
  --outdir results_custom/ \
  --genome GRCh38 \
  --skip_fastqc \
  --skip_multiqc
```

## Sample Sheets

### Paired-End Samples

```csv
sample,fastq_1,fastq_2
WGS001,/data/WGS001_R1.fastq.gz,/data/WGS001_R2.fastq.gz
WGS002,/data/WGS002_R1.fastq.gz,/data/WGS002_R2.fastq.gz
WGS003,/data/WGS003_R1.fastq.gz,/data/WGS003_R2.fastq.gz
```

### Single-End Samples

```csv
sample,fastq_1,fastq_2
RNA001,/data/RNA001.fastq.gz,
RNA002,/data/RNA002.fastq.gz,
RNA003,/data/RNA003.fastq.gz,
```

### Mixed Sample Types

```csv
sample,fastq_1,fastq_2,sample_type
WGS001,/data/WGS001_R1.fastq.gz,/data/WGS001_R2.fastq.gz,WGS
WES001,/data/WES001_R1.fastq.gz,/data/WES001_R2.fastq.gz,WES
RNA001,/data/RNA001_R1.fastq.gz,/data/RNA001_R2.fastq.gz,RNA
```

## Complete Workflows

### Workflow 1: {{ WORKFLOW_NAME_1 }}

**Description:** {{ WORKFLOW_DESC_1 }}

```bash
# Step 1: Prepare input
cat > samples.csv << EOF
sample,fastq_1,fastq_2
{{ SAMPLE_1 }},{{ PATH_1_R1 }},{{ PATH_1_R2 }}
{{ SAMPLE_2 }},{{ PATH_2_R1 }},{{ PATH_2_R2 }}
EOF

# Step 2: Run pipeline
nextflow run main.nf \
  --input samples.csv \
  --outdir {{ WORKFLOW_1_OUTPUT }}/ \
  --genome {{ WORKFLOW_1_GENOME }} \
  --{{ WORKFLOW_1_PARAM1 }} {{ WORKFLOW_1_VALUE1 }} \
  --{{ WORKFLOW_1_PARAM2 }} {{ WORKFLOW_1_VALUE2 }}

# Step 3: Check results
ls {{ WORKFLOW_1_OUTPUT }}/
```

### Workflow 2: {{ WORKFLOW_NAME_2 }}

**Description:** {{ WORKFLOW_DESC_2 }}

```yaml
# params_{{ WORKFLOW_2 }}.yml
input: 'samples.csv'
outdir: '{{ WORKFLOW_2_OUTPUT }}/'
genome: '{{ WORKFLOW_2_GENOME }}'
{{ WORKFLOW_2_PARAM1 }}: {{ WORKFLOW_2_VALUE1 }}
{{ WORKFLOW_2_PARAM2 }}: {{ WORKFLOW_2_VALUE2 }}
max_cpus: 32
max_memory: '128.GB'
```

```bash
nextflow run main.nf -params-file params_{{ WORKFLOW_2 }}.yml
```

## Testing Examples

### Quick Test

```bash
nextflow run main.nf -profile test
```

### Test with Custom Data

```bash
nextflow run main.nf \
  --input test_samples.csv \
  --outdir test_results/ \
  -profile test,docker
```

## Troubleshooting Examples

### Resume Failed Run

```bash
nextflow run main.nf \
  --input samples.csv \
  --outdir results/ \
  -resume
```

### Debug Mode

```bash
nextflow run main.nf \
  --input samples.csv \
  --outdir results/ \
  -with-trace \
  -with-report \
  -with-timeline
```

For more troubleshooting tips, see the [troubleshooting guide](/guide/troubleshooting).