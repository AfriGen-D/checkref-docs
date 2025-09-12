# Parameters Reference

Complete reference for all pipeline parameters in {{ PROJECT_NAME }}.

## Input/Output Parameters

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `--input` | `string` | Path to sample sheet CSV file |
| `--outdir` | `string` | Output directory for results |

### Optional I/O Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--email` | `string` | `null` | Email for completion notification |
| `--email_on_fail` | `string` | `null` | Email on pipeline failure |
| `--plaintext_email` | `boolean` | `false` | Send plain-text emails |
| `--max_multiqc_email_size` | `string` | `25.MB` | Maximum MultiQC email attachment size |

## Reference Genome Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--genome` | `string` | `{{ DEFAULT_GENOME }}` | Reference genome assembly |
| `--fasta` | `string` | `null` | Path to reference FASTA file |
| `--fasta_fai` | `string` | `null` | Path to reference FASTA index |
| `--bwa_index` | `string` | `null` | Path to BWA index |
| `--save_reference` | `boolean` | `false` | Save generated reference files |

## {{ ANALYSIS_TYPE }} Parameters

### Core Analysis Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--{{ PARAM1 }}` | `{{ TYPE1 }}` | `{{ DEFAULT1 }}` | {{ DESCRIPTION1 }} |
| `--{{ PARAM2 }}` | `{{ TYPE2 }}` | `{{ DEFAULT2 }}` | {{ DESCRIPTION2 }} |
| `--{{ PARAM3 }}` | `{{ TYPE3 }}` | `{{ DEFAULT3 }}` | {{ DESCRIPTION3 }} |

### Quality Control Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--skip_fastqc` | `boolean` | `false` | Skip FastQC quality control |
| `--skip_multiqc` | `boolean` | `false` | Skip MultiQC report generation |
| `--fastqc_args` | `string` | `''` | Additional FastQC arguments |

### Filtering Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--min_mapping_quality` | `integer` | `20` | Minimum mapping quality |
| `--remove_duplicates` | `boolean` | `true` | Remove duplicate reads |
| `--min_base_quality` | `integer` | `20` | Minimum base quality |

## Resource Parameters

### Computational Resources

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--max_cpus` | `integer` | `16` | Maximum number of CPUs |
| `--max_memory` | `string` | `128.GB` | Maximum memory allocation |
| `--max_time` | `string` | `240.h` | Maximum execution time |

### Process-Specific Resources

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--{{ PROCESS1 }}_cpus` | `integer` | `4` | CPUs for {{ PROCESS1 }} |
| `--{{ PROCESS1 }}_memory` | `string` | `16.GB` | Memory for {{ PROCESS1 }} |
| `--{{ PROCESS2 }}_cpus` | `integer` | `8` | CPUs for {{ PROCESS2 }} |
| `--{{ PROCESS2 }}_memory` | `string` | `32.GB` | Memory for {{ PROCESS2 }} |

## Advanced Parameters

### Container Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--singularity_pull_docker_container` | `boolean` | `false` | Pull Singularity from Docker Hub |
| `--docker_registry` | `string` | `'quay.io'` | Docker registry to use |

### Reporting Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--multiqc_config` | `string` | `null` | Custom MultiQC config file |
| `--multiqc_title` | `string` | `null` | Custom MultiQC report title |
| `--custom_config_version` | `string` | `'master'` | nf-core/configs version |

### Execution Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--custom_config_base` | `string` | `null` | Custom config base directory |
| `--hostnames` | `string` | `null` | Institutional config hostname |
| `--config_profile_name` | `string` | `null` | Institutional config name |
| `--config_profile_description` | `string` | `null` | Institutional config description |

## Parameter Files

### YAML Format

```yaml
# params.yml
input: 'samples.csv'
outdir: 'results/'
genome: 'GRCh38'
{{ PARAM1 }}: {{ VALUE1 }}
{{ PARAM2 }}: {{ VALUE2 }}
max_cpus: 16
max_memory: '64.GB'
```

### JSON Format

```json
{
  "input": "samples.csv",
  "outdir": "results/",
  "genome": "GRCh38",
  "{{ PARAM1 }}": "{{ VALUE1 }}",
  "{{ PARAM2 }}": "{{ VALUE2 }}",
  "max_cpus": 16,
  "max_memory": "64.GB"
}
```

## Parameter Validation

The pipeline validates parameters to ensure:

- Required parameters are provided
- File paths exist and are accessible
- Numeric values are within valid ranges
- Boolean values are properly formatted
- Memory/time specifications use valid units

### Common Validation Errors

**Missing required parameter:**
```
ERROR ~ Parameter '--input' is required but was not provided
```

**Invalid file path:**
```
ERROR ~ Input file does not exist: /path/to/missing/file.csv
```

**Invalid memory format:**
```
ERROR ~ Invalid memory specification: '64GB' (should be '64.GB')
```

## Examples

### Basic Parameter Set

```bash
nextflow run main.nf \
  --input samples.csv \
  --outdir results/ \
  --genome GRCh38
```

### Advanced Parameter Set

```bash
nextflow run main.nf \
  --input samples.csv \
  --outdir results/ \
  --genome GRCh38 \
  --{{ PARAM1 }} {{ VALUE1 }} \
  --{{ PARAM2 }} {{ VALUE2 }} \
  --max_cpus 32 \
  --max_memory '128.GB' \
  --email user@example.com
```

For complete examples, see the [Examples section](/examples/).