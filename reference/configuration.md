# Advanced Configuration

Comprehensive guide to configuring CheckRef for different environments and use cases.

## Configuration Files

### Primary Configuration

CheckRef uses `nextflow.config` for all configuration:

```groovy
// nextflow.config
params {
    // Core parameters
    targetVcfs = null
    referenceDir = null
    outputDir = "results"
    fixMethod = "remove"
    
    // Resource limits
    maxCpus = 4
    maxMemory = 8.GB
    maxTime = 24.h
}
```

### Custom Configuration Files

Create environment-specific configurations:

```bash
# Create custom config
cat > my_config.config << EOF
params {
    maxCpus = 16
    maxMemory = 64.GB
    outputDir = "/scratch/checkref_results"
}

process {
    withName: CHECK_ALLELE_SWITCH {
        cpus = 8
        memory = 32.GB
        time = 12.h
    }
}
EOF

# Use custom config
nextflow run main.nf -c my_config.config [options]
```

## Resource Configuration

### Global Resource Limits

```groovy
params {
    maxCpus = 16        // Maximum CPUs per process
    maxMemory = 64.GB   // Maximum memory per process  
    maxTime = 48.h      // Maximum time per process
}
```

### Process-Specific Resources

```groovy
process {
    // Default for all processes
    cpus = 1
    memory = 4.GB
    time = 1.h
    
    // Specific process configuration
    withName: CHECK_ALLELE_SWITCH {
        cpus = 4
        memory = 16.GB
        time = 8.h
    }
    
    withName: CORRECT_SWITCHED_SITES {
        cpus = 2
        memory = 8.GB
        time = 4.h
    }
    
    withName: VALIDATE_VCF_FILES {
        cpus = 1
        memory = 2.GB
        time = 30.min
    }
}
```

### Dynamic Resource Allocation

```groovy
process {
    withName: CHECK_ALLELE_SWITCH {
        cpus = { task.attempt == 1 ? 4 : 8 }
        memory = { task.attempt == 1 ? 16.GB : 32.GB }
        time = { task.attempt == 1 ? 4.h : 8.h }
        errorStrategy = 'retry'
        maxRetries = 2
    }
}
```

## Container Configuration

### Docker Configuration

```groovy
docker {
    enabled = true
    runOptions = '-u $(id -u):$(id -g)'
    registry = 'docker.io'
    temp = 'auto'
}

process {
    container = 'mamana/vcf-processing:latest'
    
    withName: CHECK_ALLELE_SWITCH {
        container = 'mamana/vcf-processing:v2.0'
    }
}
```

### Singularity Configuration

```groovy
singularity {
    enabled = true
    autoMounts = true
    cacheDir = '/tmp/singularity-cache'
    runOptions = '--cleanenv'
}

process {
    container = 'docker://mamana/vcf-processing:latest'
}
```

### Custom Container Registry

```groovy
docker {
    enabled = true
    registry = 'your-registry.com'
}

process {
    container = 'your-registry.com/checkref:latest'
}
```

## Executor Configuration

### Local Executor

```groovy
process {
    executor = 'local'
    cpus = 4
    memory = 16.GB
}
```

### SLURM Configuration

```groovy
process {
    executor = 'slurm'
    queue = 'normal'
    clusterOptions = '--account=project123 --partition=compute'
    
    withName: CHECK_ALLELE_SWITCH {
        queue = 'highmem'
        clusterOptions = '--account=project123 --partition=highmem --exclusive'
    }
}
```

### PBS/Torque Configuration

```groovy
process {
    executor = 'pbs'
    queue = 'workq'
    clusterOptions = '-A project123'
}
```

### AWS Batch Configuration

```groovy
aws {
    region = 'us-east-1'
    batch {
        cliPath = '/home/ec2-user/miniconda/bin/aws'
        jobRole = 'arn:aws:iam::123456789:role/BatchExecutionRole'
        jobQueue = 'checkref-queue'
    }
}

process {
    executor = 'awsbatch'
    container = 'mamana/vcf-processing:latest'
}
```

## Storage Configuration

### Work Directory

```groovy
workDir = '/fast/scratch/checkref_work'
```

### Output Publishing

```groovy
process {
    publishDir = [
        path: "${params.outputDir}",
        mode: 'copy',
        overwrite: true
    ]
    
    withName: CHECK_ALLELE_SWITCH {
        publishDir = [
            [path: "${params.outputDir}/analysis", mode: 'copy', pattern: "*.tsv"],
            [path: "${params.outputDir}/summaries", mode: 'copy', pattern: "*.txt"]
        ]
    }
}
```

### Temporary Directory

```groovy
env {
    TMPDIR = '/fast/tmp'
    TEMP = '/fast/tmp'
}
```

## Validation Configuration

### Validation Parameters

```groovy
params {
    // Validation control
    skip_validation = false
    validation_af_threshold = 0.05
    validation_fold_threshold = 2.0
    
    // Quality thresholds
    min_overlap_rate = 0.5
    max_switch_rate = 0.3
}
```

### Custom Validation Thresholds

```groovy
// Strict validation
params {
    validation_af_threshold = 0.01
    validation_fold_threshold = 1.5
    min_overlap_rate = 0.7
}

// Permissive validation  
params {
    validation_af_threshold = 0.10
    validation_fold_threshold = 3.0
    min_overlap_rate = 0.3
}
```

## Reporting Configuration

### Execution Reports

```groovy
report {
    enabled = true
    file = "${params.outputDir}/reports/execution_report.html"
    overwrite = true
}

timeline {
    enabled = true
    file = "${params.outputDir}/reports/timeline_report.html"
    overwrite = true
}

dag {
    enabled = true
    file = "${params.outputDir}/reports/dag_report.html"
    overwrite = true
}

trace {
    enabled = true
    file = "${params.outputDir}/reports/trace.txt"
    fields = 'task_id,hash,native_id,process,tag,name,status,exit,module,container,cpus,time,disk,memory,attempt,submit,start,complete,duration,realtime,queue,%cpu,%mem,rss,vmem,peak_rss,peak_vmem,rchar,wchar,syscr,syscw,read_bytes,write_bytes'
}
```

### Custom Report Location

```groovy
report.file = "/shared/reports/checkref_${workflow.runName}.html"
timeline.file = "/shared/reports/timeline_${workflow.runName}.html"
```

## Environment Configuration

### Environment Variables

```groovy
env {
    PYTHONPATH = '/opt/conda/lib/python3.9/site-packages'
    BCFTOOLS_PLUGINS = '/usr/local/libexec/bcftools'
    TMPDIR = '/fast/tmp'
}
```

### Module Loading (HPC)

```groovy
process {
    module = ['bcftools/1.17', 'python/3.9']
    
    withName: CHECK_ALLELE_SWITCH {
        module = ['bcftools/1.17', 'python/3.9', 'R/4.2']
    }
}
```

## Performance Tuning

### Memory Optimization

```groovy
process {
    withName: CHECK_ALLELE_SWITCH {
        memory = { vcf_size < 1.GB ? 8.GB : 
                  vcf_size < 5.GB ? 16.GB : 32.GB }
    }
}
```

### CPU Optimization

```groovy
process {
    cpus = { Math.min(params.maxCpus, 
                     task.attempt * 2) }
}
```

### I/O Optimization

```groovy
process {
    scratch = true  // Use local scratch space
    stageInMode = 'symlink'
    stageOutMode = 'rsync'
}
```

## Error Handling

### Retry Configuration

```groovy
process {
    errorStrategy = 'retry'
    maxRetries = 3
    
    withName: CHECK_ALLELE_SWITCH {
        errorStrategy = { task.exitStatus in [130,143,137,104,134,139] ? 
                         'retry' : 'finish' }
        maxRetries = 2
    }
}
```

### Ignore Errors

```groovy
process {
    withName: VALIDATE_VCF_FILES {
        errorStrategy = 'ignore'
    }
}
```

## Configuration Examples

### High-Performance Setup

```groovy
// hpc_config.config
params {
    maxCpus = 32
    maxMemory = 128.GB
    maxTime = 48.h
}

process {
    executor = 'slurm'
    queue = 'highmem'
    clusterOptions = '--account=project123 --exclusive'
    
    withName: CHECK_ALLELE_SWITCH {
        cpus = 16
        memory = 64.GB
        time = 24.h
    }
}

singularity {
    enabled = true
    autoMounts = true
    cacheDir = '/shared/singularity-cache'
}
```

### Cloud Setup

```groovy
// cloud_config.config
aws {
    region = 'us-east-1'
    batch.jobQueue = 'checkref-queue'
}

process {
    executor = 'awsbatch'
    container = 'mamana/vcf-processing:latest'
    
    withName: CHECK_ALLELE_SWITCH {
        cpus = 8
        memory = 32.GB
    }
}

docker.enabled = true
```

## Configuration Validation

### Test Configuration

```bash
# Dry run to test configuration
nextflow run main.nf -profile test -c my_config.config

# Check configuration parsing
nextflow config -profile hpc
```

### Debug Configuration

```groovy
// Enable debug logging
trace.enabled = true
trace.file = 'trace.txt'

// Verbose logging
process.echo = true
```

## Next Steps

- **[Parameters](./parameters)** - Complete parameter reference
- **[Profiles](./profiles)** - Execution profile options
- **[Troubleshooting](/docs/troubleshooting)** - Configuration troubleshooting
