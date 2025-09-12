---
layout: home

hero:
  name: "CheckRef"
  text: "Allele Switch Detection & Correction"
  tagline: "A Nextflow pipeline for identifying and fixing allele switches between target VCF files and reference panels"
  actions:
    - theme: brand
      text: Get Started
      link: /tutorials/quick-start
    - theme: alt
      text: View on GitHub
      link: https://github.com/AfriGen-D/checkref

features:
  - icon: 🔍
    title: Allele Switch Detection
    details: Automatically identifies variants with flipped REF/ALT alleles, strand issues, and other inconsistencies.
  - icon: 🔧
    title: Two Correction Methods
    details: Choose to remove problematic variants or correct them by swapping alleles.
  - icon: 📊
    title: Multi-File Processing
    details: Process multiple chromosomes simultaneously with automatic matching to reference panels.
  - icon: 🧬
    title: Population Genomics Ready
    details: Works with any reference panel and supports diverse population datasets.
  - icon: ⚡
    title: Fast & Scalable
    details: Parallel processing with Nextflow for efficient analysis of large cohorts.
  - icon: 📋
    title: Comprehensive QC
    details: Detailed reports and quality metrics to validate your corrections.
---

## Quick Start

Get started with CheckRef in just a few commands:

```bash
# Clone the repository
git clone https://github.com/AfriGen-D/checkref.git
cd checkref

# Run with your data
nextflow run main.nf \
  --targetVcfs your_file.vcf.gz \
  --referenceDir /path/to/reference/panels/ \
  --outputDir results \
  -profile singularity
```

## What CheckRef Does

This Nextflow pipeline provides:

- **Allele Switch Detection**: Compares your VCF against reference panels to identify orientation issues
- **Smart Correction**: Choose to remove problematic variants or fix them by swapping REF/ALT alleles
- **Quality Assessment**: Comprehensive statistics and validation reports
- **Multi-File Support**: Process multiple chromosomes simultaneously with automatic matching

## Documentation

- [**Tutorials**](/tutorials/) - Step-by-step learning exercises
- [**Documentation**](/docs/) - Comprehensive reference material
- [**Parameters**](/api/parameters) - Complete parameter documentation
- [**Examples**](/examples/) - Ready-to-use configurations

## Requirements

- Nextflow ≥ 21.04.0
- Docker, Singularity, or Conda
- Target VCF files (bgzipped and indexed)
- Reference panel legend files

## Support

- [GitHub Issues](https://github.com/AfriGen-D/checkref/issues)
- [Helpdesk](https://helpdesk.afrigen-d.org)
- [AfriGen-D](https://afrigen-d.org)