# Execution Profiles

CheckRef includes several predefined execution profiles to run the pipeline in different computing environments.

## Available Profiles

### Standard Profile (Default)

**Usage:**
```bash
nextflow run main.nf -profile standard [options]
# or simply (standard is default)
nextflow run main.nf [options]
```

**Configuration:**
- **Executor:** Local execution
- **Best for:** Small datasets, testing, local development
- **Resource limits:** Uses system defaults

**Example:**
```bash
nextflow run main.nf \
  --targetVcfs chr22.vcf.gz \
  --referenceDir /ref/panels/ \
  -profile standard
```

### HPC Profile

**Usage:**
```bash
nextflow run main.nf -profile hpc [options]
```

**Configuration:**
- **Executor:** SLURM scheduler
- **Queue:** normal
- **Best for:** High-performance computing clusters
- **Account:** Requires cluster account configuration

**Example:**
```bash
nextflow run main.nf \
  --targetVcfs "chr*.vcf.gz" \
  --referenceDir /ref/panels/ \
  -profile hpc
```

**Customization:**
```bash
# Modify cluster options in nextflow.config
process.clusterOptions = '--account=your_project --partition=compute'
```

### Docker Profile

**Usage:**
```bash
nextflow run main.nf -profile docker [options]
```

**Configuration:**
- **Container engine:** Docker
- **Auto-pull:** Enabled
- **Best for:** Consistent environments, cloud computing

**Prerequisites:**
```bash
# Install Docker
sudo apt-get install docker.io
sudo usermod -aG docker $USER

# Verify installation
docker --version
docker run hello-world
```

**Example:**
```bash
nextflow run main.nf \
  --targetVcfs chr22.vcf.gz \
  --referenceDir /ref/panels/ \
  -profile docker
```

### Singularity Profile

**Usage:**
```bash
nextflow run main.nf -profile singularity [options]
```

**Configuration:**
- **Container engine:** Singularity/Apptainer
- **Auto-mounts:** Enabled
- **Best for:** HPC environments, shared systems

**Prerequisites:**
```bash
# Install Singularity (varies by system)
# On Ubuntu/Debian:
sudo apt-get install singularity-container

# Verify installation
singularity --version
```

**Example:**
```bash
nextflow run main.nf \
  --targetVcfs "chr{20,21,22}.vcf.gz" \
  --referenceDir /ref/panels/ \
  -profile singularity
```

## Profile Combinations

### Multiple Profiles

Combine profiles for specific environments:

```bash
# Docker + HPC
nextflow run main.nf -profile hpc,docker [options]

# Singularity + HPC  
nextflow run main.nf -profile hpc,singularity [options]
```

### Custom Profile Creation

Create custom profiles in `nextflow.config`:

```groovy
profiles {
    mylab {
        process.executor = 'slurm'
        process.queue = 'gpu'
        process.clusterOptions = '--account=lab123 --gres=gpu:1'
        singularity.enabled = true
        singularity.autoMounts = true
    }
    
    cloud {
        process.executor = 'awsbatch'
        aws.region = 'us-east-1'
        aws.batch.cliPath = '/home/ec2-user/miniconda/bin/aws'
        docker.enabled = true
    }
}
```

## Resource Configuration by Profile

### Standard Profile Resources

```groovy
process {
    cpus = 1
    memory = 4.GB
    time = 1.h
}
```

### HPC Profile Resources

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

## Container Images

### Default Container

All profiles use the same container image:
- **Image:** `mamana/vcf-processing:latest`
- **Contents:** bcftools, Python 3, required dependencies
- **Registry:** Docker Hub

### Custom Containers

Specify custom containers:

```groovy
process {
    withName: CHECK_ALLELE_SWITCH {
        container = 'your-registry/custom-image:v1.0'
    }
}
```

## Profile Selection Guide

| Use Case | Recommended Profile | Notes |
|----------|-------------------|-------|
| **Testing/Development** | `standard` | Quick local testing |
| **Small datasets (<1GB)** | `standard` or `docker` | Local execution sufficient |
| **Large datasets (>1GB)** | `hpc,singularity` | Use cluster resources |
| **Cloud computing** | `docker` | Consistent environments |
| **Shared HPC systems** | `hpc,singularity` | No root access needed |
| **Production pipelines** | `hpc,singularity` | Scalable and reproducible |

## Troubleshooting Profiles

### Docker Issues

```bash
# Permission denied
sudo usermod -aG docker $USER
newgrp docker

# Container pull failures
docker pull mamana/vcf-processing:latest

# Check Docker daemon
sudo systemctl status docker
```

### Singularity Issues

```bash
# Cache directory permissions
export SINGULARITY_CACHEDIR=/tmp/singularity-cache
mkdir -p $SINGULARITY_CACHEDIR

# Manual container pull
singularity pull docker://mamana/vcf-processing:latest
```

### HPC Issues

```bash
# Check SLURM availability
sinfo
squeue -u $USER

# Verify account access
sacctmgr show user $USER

# Test job submission
sbatch --wrap="echo 'test job'" --account=your_project
```

## Performance Optimization

### Resource Tuning

Adjust resources based on data size:

```bash
# Small datasets
--max_cpus 4 --max_memory 8.GB

# Large datasets  
--max_cpus 16 --max_memory 64.GB

# Very large datasets
--max_cpus 32 --max_memory 128.GB
```

### Profile-Specific Optimization

```groovy
profiles {
    performance {
        process.cpus = 8
        process.memory = 32.GB
        process.time = 24.h
        
        withName: CHECK_ALLELE_SWITCH {
            cpus = 16
            memory = 64.GB
        }
    }
}
```

## Next Steps

- **[Configuration](./configuration)** - Advanced configuration options
- **[Parameters](./parameters)** - Complete parameter reference
- **[Troubleshooting](/docs/troubleshooting)** - Solve execution issues
