# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a VitePress documentation template for AfriGen-D (African Genomics Database) Nextflow bioinformatics pipelines. It's designed to be a reusable template that can be customized for different genomics projects and pipelines.

### Template Architecture

The repository serves as a documentation template with extensive use of template variables:
- **Template Variables**: All content uses `{{ VARIABLE_NAME }}` placeholder syntax for customization
- **Purpose**: Creates standardized, branded documentation sites for AfriGen-D genomics tools
- **Target**: Nextflow pipelines, bioinformatics tools, and genomics data resources
- **Deployment**: GitHub Pages with automated workflows

## Development Commands

### Documentation Development
```bash
# Start development server with hot reload
npm run docs:dev

# Build static documentation site
npm run docs:build  

# Preview built documentation locally
npm run docs:preview
```

### Dependencies
```bash
# Install dependencies (VitePress only)
npm install

# Clean install
npm ci
```

## Architecture & Structure

### Core Configuration
- **VitePress Config**: `.vitepress/config.ts` - Main configuration with navigation, sidebar, and theming
- **Package Configuration**: `package.json` - Template with AfriGen-D branding and GitHub Pages setup
- **GitHub Actions**: `.github/workflows/deploy.yml` - Automated GitHub Pages deployment

### Documentation Structure
```
/
├── index.md              # Homepage with hero section and features
├── guide/                # User guides and tutorials
│   ├── getting-started.md
│   ├── installation.md
│   ├── quick-start.md
│   ├── configuration.md
│   └── input-files.md
├── api/                  # Technical reference
│   └── parameters.md     # Complete parameter documentation
├── examples/             # Usage examples and workflows
│   └── index.md
├── workflow/             # Technical pipeline documentation  
│   └── index.md
└── public/               # Static assets
```

### Template Variable System
The entire template uses a systematic placeholder approach:
- **Project Identity**: `{{ PROJECT_NAME }}`, `{{ REPO_NAME }}`, `{{ PROJECT_DESCRIPTION }}`
- **Content Variables**: `{{ HERO_TEXT }}`, `{{ PROJECT_TAGLINE }}`, `{{ PROCESSING_DESCRIPTION }}`
- **Technical Variables**: `{{ MIN_NEXTFLOW_VERSION }}`, `{{ DEFAULT_GENOME }}`
- **Parameter Variables**: `{{ PARAM1 }}`, `{{ VALUE1 }}`, etc.

### Content Patterns
1. **Nextflow Focus**: All content assumes Nextflow pipeline architecture
2. **African Genomics**: Optimized for African genomic datasets and populations  
3. **Bioinformatics Workflow**: Standard QC → Processing → Analysis → Reporting pipeline
4. **Container Support**: Docker/Singularity/Conda integration
5. **HPC/Cloud Ready**: Resource allocation and scaling guidance

### Navigation & User Experience
- **Progressive Disclosure**: From getting started to advanced configuration
- **Multi-Modal**: Examples, parameters, technical details, and workflows
- **Search Enabled**: Local search functionality built-in
- **Responsive Design**: Mobile-friendly documentation
- **AfriGen-D Branding**: Consistent organizational theming

## Key Technical Details

### VitePress Features Used
- **Math Support**: Enabled for scientific notation
- **Line Numbers**: Code blocks with line numbering
- **Mermaid Diagrams**: Workflow visualization support
- **Local Search**: Built-in content indexing
- **GitHub Integration**: Direct repository linking

### GitHub Pages Integration
- **Base Path**: Configured for `/{REPO_NAME}/` deployment
- **Automated Deployment**: Triggers on docs changes or workflow updates
- **Build Process**: Node.js 18, npm ci, VitePress build
- **Artifact Handling**: Proper GitHub Pages artifact upload

### AfriGen-D Standards
- **Helpdesk Integration**: Links to helpdesk.afrigen-d.org
- **MIT License**: Standard open source licensing
- **GitHub Organization**: AfriGen-D namespace integration
- **Support Channels**: GitHub Issues + helpdesk routing

## Template Customization Notes

When using this template:
1. **Variable Replacement**: All `{{ VARIABLE }}` placeholders need systematic replacement
2. **Content Adaptation**: Modify workflow descriptions for specific pipeline purposes
3. **Parameter Updates**: Replace generic parameters with actual pipeline parameters
4. **Resource Requirements**: Update memory/CPU requirements based on actual needs
5. **Example Data**: Replace placeholder examples with real sample data paths

## Important File Dependencies

- `.vitepress/config.ts`: Navigation structure mirrors actual content organization
- `index.md`: Hero section drives first-user experience
- `guide/getting-started.md`: Critical onboarding flow
- `api/parameters.md`: Must match actual Nextflow pipeline parameters
- `.github/workflows/deploy.yml`: Assumes docs/ subdirectory (may need adjustment)