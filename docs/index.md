# CheckRef Documentation

Comprehensive documentation for the CheckRef pipeline - your complete reference for understanding, configuring, and troubleshooting CheckRef analyses.

## Overview

CheckRef is a Nextflow pipeline that identifies and corrects allele switches between target VCF files and reference panel legend files. This documentation provides detailed information about every aspect of CheckRef usage.

## Documentation Sections

### Core Concepts
- [**Understanding Results**](./understanding-results) - Complete guide to interpreting CheckRef output and validation reports
- [**Correction Methods**](./correction-methods) - Detailed comparison of remove vs correct approaches with automated validation
- [**Single File Analysis**](./single-file) - In-depth single-chromosome processing

### Processing Workflows  
- [**Multi-File Processing**](./multi-file) - Comprehensive multi-chromosome analysis
- [**Quality Control**](./quality-control) - Complete QC procedures, validation workflows, and automated quality assessment

### Reference Materials
- [**Troubleshooting**](./troubleshooting) - Complete problem-solving guide including validation issues

## New Features

### Automated Validation
CheckRef now includes built-in validation for correction results:
- **Automatic quality assessment** when using `--fixMethod correct`
- **Allele frequency validation** to detect unexpected changes
- **Switch accuracy verification** to confirm corrections were applied
- **Interactive HTML reports** for comprehensive result review

### Enhanced Quality Control
- Real-time validation during pipeline execution
- Comprehensive validation reports with actionable recommendations
- Configurable validation thresholds for different data types
- Integration with existing QC workflows

## When to Use Documentation vs Tutorials

**Use Documentation when you need to:**
- Understand how CheckRef processes work internally
- Get comprehensive parameter information
- Implement quality control procedures
- Troubleshoot complex issues
- Reference detailed command options

**Use Tutorials when you want to:**
- Get started quickly with CheckRef
- Follow step-by-step workflows
- Learn specific techniques
- Practice with guided examples

## Documentation Organization

Each documentation page provides:
- **Comprehensive coverage** of the topic
- **Detailed examples** with real commands
- **Advanced configuration** options
- **Best practices** and recommendations
- **Troubleshooting** for common issues
- **Cross-references** to related topics

## Getting Started with Documentation

If you're new to CheckRef documentation:

1. **Start with [Understanding Results](./understanding-results)** - Learn how to interpret CheckRef output
2. **Review [Correction Methods](./correction-methods)** - Understand your correction options  
3. **Explore [Quality Control](./quality-control)** - Learn validation procedures
4. **Reference [Troubleshooting](./troubleshooting)** - When you need to solve problems

## Documentation vs Tutorial Content

| Documentation | Tutorials |
|---------------|-----------|
| Comprehensive reference | Step-by-step learning |
| Detailed explanations | Focused exercises |
| All options covered | Essential steps only |
| Problem-solving focus | Achievement-oriented |
| Reference material | Learning material |

## Contributing to Documentation

Found an issue or want to improve the documentation? 

- Visit our [GitHub repository](https://github.com/AfriGen-D/checkref)
- Submit issues or suggestions
- Contact [helpdesk.afrigen-d.org](https://helpdesk.afrigen-d.org)

## Related Resources

- [**Tutorials**](/tutorials/) - Step-by-step learning exercises
- [**API Reference**](/api/parameters) - Complete parameter documentation
- [**Examples**](/examples/) - Ready-to-use examples
- [**Guide**](/guide/getting-started) - Getting started information
