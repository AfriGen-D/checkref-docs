---
layout: home

hero:
  name: "{{ PROJECT_NAME }}"
  text: "{{ HERO_TEXT }}"
  tagline: "{{ PROJECT_TAGLINE }}"
  actions:
    - theme: brand
      text: Get Started
      link: /guide/getting-started
    - theme: alt
      text: View on GitHub
      link: https://github.com/AfriGen-D/{{ REPO_NAME }}

features:
  - icon: 🔬
    title: Nextflow Workflow
    details: Scalable and reproducible bioinformatics pipeline built with Nextflow.
  - icon: 🧬
    title: African Genomics Focus
    details: Optimized for population-specific analysis of African genomic datasets.
  - icon: ☁️
    title: Cloud-Ready
    details: Run on local systems, HPC clusters, or cloud platforms with ease.
  - icon: 📊
    title: Rich Outputs
    details: Generate comprehensive reports, QC metrics, and visualizations.
  - icon: 🔧
    title: Configurable
    details: Flexible parameters and profiles for different use cases.
  - icon: 📝
    title: Well Documented
    details: Complete documentation with examples and best practices.
---

## Quick Start

Run the pipeline with minimal configuration:

```bash
# Clone the repository
git clone https://github.com/AfriGen-D/{{ REPO_NAME }}.git
cd {{ REPO_NAME }}

# Run with test data
nextflow run main.nf -profile test

# Run with your data
nextflow run main.nf --input samples.csv --outdir results/
```

## Pipeline Overview

This Nextflow pipeline provides:

- **Quality Control**: Comprehensive QC of input data
- **Processing**: {{ PROCESSING_DESCRIPTION }}
- **Analysis**: {{ ANALYSIS_DESCRIPTION }}
- **Reporting**: Automated report generation with MultiQC

## Documentation

- [Getting Started](/guide/getting-started) - Installation and basic usage
- [Parameters](/api/parameters) - Complete parameter reference
- [Examples](/examples/) - Example configurations and use cases
- [Workflow Details](/workflow/) - Technical pipeline documentation

## Requirements

- Nextflow ≥ {{ MIN_NEXTFLOW_VERSION }}
- Docker, Singularity, or Conda
- {{ ADDITIONAL_REQUIREMENTS }}

## Support

- [GitHub Issues](https://github.com/AfriGen-D/{{ REPO_NAME }}/issues)
- [Helpdesk](https://helpdesk.afrigen-d.org)
- [AfriGen-D](https://afrigen-d.org)