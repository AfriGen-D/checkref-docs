# Resource Usage and Optimization

Comprehensive guide to understanding and optimizing resource usage in CheckRef workflows.

## Resource Requirements Overview

### Minimum System Requirements

| Component | Minimum | Recommended | High-Performance |
|-----------|---------|-------------|------------------|
| **CPU** | 2 cores | 8 cores | 32+ cores |
| **RAM** | 4 GB | 16 GB | 64+ GB |
| **Storage** | 10 GB | 50 GB | 500+ GB |
| **Network** | 1 Mbps | 100 Mbps | 1+ Gbps |

### Process-Specific Resource Usage

| Process | CPU | Memory | Disk I/O | Network | Duration |
|---------|-----|--------|----------|---------|----------|
| **VALIDATE_VCF_FILES** | Low | Low | Medium | None | 1-5 min |
| **CHECK_ALLELE_SWITCH** | High | Medium | High | None | 10-60 min |
| **REMOVE_SWITCHED_SITES** | Medium | Low | High | None | 5-15 min |
| **CORRECT_SWITCHED_SITES** | Medium | Medium | High | None | 10-30 min |
| **VERIFY_CORRECTIONS** | High | Medium | High | None | 10-60 min |
| **CREATE_SUMMARY** | Low | Low | Low | None | 1-2 min |

## Resource Scaling Guidelines

### Data Size Impact

#### Small Datasets (< 1 GB VCF)
```groovy
process {
    cpus = 2
    memory = 4.GB
    time = 2.h
    
    withName: CHECK_ALLELE_SWITCH {
        cpus = 4
        memory = 8.GB
        time = 4.h
    }
}
```

#### Medium Datasets (1-10 GB VCF)
```groovy
process {
    cpus = 4
    memory = 16.GB
    time = 8.h
    
    withName: CHECK_ALLELE_SWITCH {
        cpus = 8
        memory = 32.GB
        time = 12.h
    }
}
```

#### Large Datasets (> 10 GB VCF)
```groovy
process {
    cpus = 8
    memory = 32.GB
    time = 24.h
    
    withName: CHECK_ALLELE_SWITCH {
        cpus = 16
        memory = 64.GB
        time = 48.h
    }
}
```

### Dynamic Resource Allocation

```groovy
process {
    // Scale resources based on file size
    cpus = { 
        def vcf_size = file(params.targetVcfs).size()
        vcf_size < 1.GB ? 2 : 
        vcf_size < 10.GB ? 4 : 8
    }
    
    memory = { 
        def vcf_size = file(params.targetVcfs).size()
        vcf_size < 1.GB ? 4.GB : 
        vcf_size < 10.GB ? 16.GB : 32.GB
    }
    
    // Retry with more resources
    errorStrategy = 'retry'
    maxRetries = 2
    
    cpus = { task.attempt == 1 ? 4 : 8 }
    memory = { task.attempt == 1 ? 16.GB : 32.GB }
}
```

## Memory Optimization

### Memory Usage Patterns

#### VALIDATE_VCF_FILES
- **Peak Usage:** 100-500 MB
- **Pattern:** Constant low usage
- **Optimization:** Minimal memory needed

#### CHECK_ALLELE_SWITCH
- **Peak Usage:** 2-32 GB (depends on VCF size)
- **Pattern:** Gradual increase during variant loading
- **Optimization:** Stream processing for large files

```python
# Memory-efficient variant processing
def process_variants_streaming(vcf_file, legend_file):
    """Process variants in chunks to minimize memory usage"""
    chunk_size = 10000  # Process 10k variants at a time
    
    with pysam.VariantFile(vcf_file) as vcf:
        variant_chunk = []
        
        for variant in vcf:
            variant_chunk.append(variant)
            
            if len(variant_chunk) >= chunk_size:
                process_chunk(variant_chunk, legend_file)
                variant_chunk = []  # Clear memory
        
        # Process remaining variants
        if variant_chunk:
            process_chunk(variant_chunk, legend_file)
```

#### CORRECT_SWITCHED_SITES
- **Peak Usage:** 1-8 GB
- **Pattern:** Steady usage during VCF processing
- **Optimization:** Line-by-line processing

### Memory Monitoring

```bash
# Monitor memory usage during execution
#!/bin/bash
PROCESS_NAME="CHECK_ALLELE_SWITCH"
LOG_FILE="memory_usage.log"

while true; do
    MEMORY=$(ps aux | grep $PROCESS_NAME | grep -v grep | awk '{sum+=$6} END {print sum/1024}')
    echo "$(date): ${MEMORY} MB" >> $LOG_FILE
    sleep 30
done
```

## CPU Optimization

### CPU Usage Patterns

| Process | CPU Pattern | Parallelization | Optimization Strategy |
|---------|-------------|-----------------|----------------------|
| **VALIDATE_VCF_FILES** | Single-threaded | File-level | Process multiple files |
| **CHECK_ALLELE_SWITCH** | Multi-threaded | Algorithm-level | Parallel comparison |
| **CORRECT_SWITCHED_SITES** | Single-threaded | None | I/O optimization |

### Parallel Processing Strategies

#### Chromosome-Level Parallelization
```groovy
// Process multiple chromosomes simultaneously
workflow {
    target_vcfs_ch = Channel.fromPath(params.targetVcfs)
        .map { vcf -> tuple(extractChromosome(vcf.name), vcf) }
    
    // Each chromosome processed in parallel
    CHECK_ALLELE_SWITCH(target_vcfs_ch)
}
```

#### Process-Level Parallelization
```groovy
// Parallel execution within processes
process CHECK_ALLELE_SWITCH {
    cpus 8
    
    script:
    """
    # Use multiple threads for comparison
    check_allele_switch.py \\
        --vcf ${vcf_file} \\
        --legend ${legend_file} \\
        --threads ${task.cpus} \\
        --output results.tsv
    """
}
```

### CPU Affinity Optimization

```groovy
// Bind processes to specific CPU cores
process {
    clusterOptions = '--cpus-per-task=8 --cpu-bind=cores'
    
    withName: CHECK_ALLELE_SWITCH {
        clusterOptions = '--cpus-per-task=16 --cpu-bind=cores --exclusive'
    }
}
```

## Storage Optimization

### I/O Patterns

#### Sequential I/O (Optimal)
- VCF file reading
- Legend file parsing
- Output file writing

#### Random I/O (Suboptimal)
- Indexed VCF access
- Reference lookups

### Storage Configuration

#### Fast Local Storage
```groovy
// Use fast local storage for work directory
workDir = '/fast/local/storage/checkref_work'

// Use network storage for final outputs
process {
    publishDir = [
        path: '/shared/network/storage/results',
        mode: 'copy'
    ]
}
```

#### Temporary File Management
```groovy
process {
    // Use local scratch space
    scratch = true
    
    // Clean up temporary files
    afterScript = 'rm -rf $TMPDIR/checkref_*'
}
```

### Disk Space Monitoring

```bash
# Monitor disk usage during workflow
monitor_disk_usage() {
    local work_dir=$1
    local threshold=80  # Alert at 80% usage
    
    while true; do
        usage=$(df $work_dir | tail -1 | awk '{print $5}' | sed 's/%//')
        
        if [ $usage -gt $threshold ]; then
            echo "WARNING: Disk usage at ${usage}% in $work_dir"
            # Clean up old work files
            find $work_dir -name "work_*" -mtime +1 -delete
        fi
        
        sleep 300  # Check every 5 minutes
    done
}
```

## Network Optimization

### Data Transfer Patterns

#### Input Data Access
- Reference panel downloads
- VCF file transfers
- Container image pulls

#### Output Data Transfer
- Result file uploads
- Report generation
- Backup operations

### Network Configuration

```groovy
// Optimize for network file systems
process {
    // Stage files locally before processing
    stageInMode = 'copy'
    stageOutMode = 'rsync'
    
    // Use compression for transfers
    publishDir = [
        path: '/network/storage/results',
        mode: 'copy',
        compress: true
    ]
}
```

## Resource Monitoring and Profiling

### Built-in Monitoring

```groovy
// Enable comprehensive resource tracking
trace {
    enabled = true
    file = 'trace.txt'
    fields = 'task_id,hash,native_id,process,tag,name,status,exit,module,container,cpus,time,disk,memory,attempt,submit,start,complete,duration,realtime,queue,%cpu,%mem,rss,vmem,peak_rss,peak_vmem'
}

report {
    enabled = true
    file = 'report.html'
}

timeline {
    enabled = true
    file = 'timeline.html'
}
```

### Custom Resource Monitoring

```bash
# Resource usage profiler
#!/bin/bash
profile_resources() {
    local process_name=$1
    local output_file=$2
    
    echo "timestamp,cpu_percent,memory_mb,disk_read_mb,disk_write_mb" > $output_file
    
    while true; do
        if pgrep -f $process_name > /dev/null; then
            stats=$(ps -o pid,pcpu,rss,comm -C $process_name | tail -n +2)
            
            while read -r pid cpu mem comm; do
                # Get I/O stats
                if [ -f /proc/$pid/io ]; then
                    read_bytes=$(grep "read_bytes" /proc/$pid/io | awk '{print $2}')
                    write_bytes=$(grep "write_bytes" /proc/$pid/io | awk '{print $2}')
                    
                    echo "$(date +%s),$cpu,$((mem/1024)),$((read_bytes/1024/1024)),$((write_bytes/1024/1024))" >> $output_file
                fi
            done <<< "$stats"
        fi
        
        sleep 10
    done
}
```

## Performance Tuning Recommendations

### General Optimization

1. **Use appropriate execution profiles**
   ```bash
   # For HPC environments
   nextflow run main.nf -profile hpc,singularity
   
   # For cloud environments
   nextflow run main.nf -profile cloud,docker
   ```

2. **Optimize work directory location**
   ```bash
   # Use fast local storage
   nextflow run main.nf -w /fast/local/storage/work
   ```

3. **Enable resume for development**
   ```bash
   # Resume from previous run
   nextflow run main.nf -resume
   ```

### Process-Specific Tuning

#### For Large VCF Files
```groovy
process {
    withName: CHECK_ALLELE_SWITCH {
        // Increase memory and time
        memory = 64.GB
        time = 48.h
        
        // Use more CPUs for parallel processing
        cpus = 16
        
        // Enable retries with more resources
        errorStrategy = 'retry'
        maxRetries = 2
        memory = { task.attempt == 1 ? 64.GB : 128.GB }
    }
}
```

#### For Many Small Files
```groovy
process {
    // Reduce per-process overhead
    memory = 2.GB
    time = 1.h
    
    // Increase parallelism
    maxForks = 20
}
```

### Environment-Specific Optimization

#### HPC Clusters
```groovy
process {
    executor = 'slurm'
    queue = 'normal'
    
    // Use exclusive nodes for large jobs
    withName: CHECK_ALLELE_SWITCH {
        queue = 'highmem'
        clusterOptions = '--exclusive --mem=0'
    }
}
```

#### Cloud Environments
```groovy
process {
    // Use spot instances for cost optimization
    machineType = 'n1-standard-8'
    preemptible = true
    
    // Scale resources dynamically
    cpus = { Math.min(8, task.attempt * 2) }
    memory = { Math.min(32.GB, task.attempt * 8.GB) }
}
```

## Troubleshooting Resource Issues

### Common Resource Problems

1. **Out of Memory Errors**
   ```bash
   # Symptoms: Process killed with exit code 137
   # Solution: Increase memory allocation
   process.memory = '32.GB'
   ```

2. **Timeout Errors**
   ```bash
   # Symptoms: Process killed with exit code 143
   # Solution: Increase time limit
   process.time = '24.h'
   ```

3. **Disk Space Issues**
   ```bash
   # Symptoms: "No space left on device"
   # Solution: Clean work directory or use larger storage
   ```

### Resource Debugging

```bash
# Debug resource allocation
nextflow run main.nf -with-trace trace.txt

# Analyze resource usage
awk -F'\t' '{print $4, $10, $11, $12}' trace.txt | \
    sort -k2 -nr | head -10  # Top 10 memory users
```

## Related Documentation

- [Workflow Overview](./index) - High-level architecture
- [Process Flow](./process-flow) - Detailed process descriptions
- [Configuration](/reference/configuration) - Advanced configuration options
