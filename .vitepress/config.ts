import { defineConfig } from 'vitepress'

export default defineConfig({
  title: 'CheckRef',
  description: 'A Nextflow pipeline for checking allele switches between target VCF and reference panel legend files',
  base: '/checkref-docs/',
  lang: 'en-US',
  ignoreDeadLinks: true,

  themeConfig: {
    logo: '/logo.png',

    nav: [
      { text: 'Home', link: '/' },
      { text: 'Documentation', link: '/docs/' },
      { text: 'Reference', link: '/reference/' },
      { text: 'Tutorials', link: '/tutorials/' },
      { text: 'Examples', link: '/examples/' }
    ],

    sidebar: {
      '/guide/': [
        {
          text: 'Getting Started',
          collapsed: false,
          items: [
            { text: 'Introduction', link: '/guide/getting-started' },
            { text: 'Installation', link: '/guide/installation' },
            { text: 'Quick Start', link: '/guide/quick-start' },
            { text: 'Configuration', link: '/guide/configuration' }
          ]
        },
        {
          text: 'Advanced Usage',
          collapsed: false,
          items: [
            { text: 'Input Files', link: '/guide/input-files' },
            { text: 'Configuration', link: '/guide/configuration' }
          ]
        }
      ],
      '/reference/': [
        {
          text: 'Command Line Reference',
          collapsed: false,
          items: [
            { text: 'Reference Overview', link: '/reference/' },
            { text: 'Parameters', link: '/reference/parameters' },
            { text: 'Profiles', link: '/reference/profiles' },
            { text: 'Configuration', link: '/reference/configuration' }
          ]
        }
      ],
      '/workflow/': [
        {
          text: 'Workflow Details',
          collapsed: false,
          items: [
            { text: 'Overview', link: '/workflow/' },
            { text: 'Process Flow', link: '/workflow/process-flow' },
            { text: 'Subworkflows', link: '/workflow/subworkflows' },
            { text: 'Resource Usage', link: '/workflow/resources' }
          ]
        }
      ],
      '/tutorials/': [
        {
          text: 'Step-by-Step Tutorials',
          collapsed: false,
          items: [
            { text: 'Tutorial Overview', link: '/tutorials/' },
            { text: 'Quick Start (10 min)', link: '/tutorials/quick-start' },
            { text: 'Multi-File Tutorial', link: '/tutorials/multi-file-tutorial' },
            { text: 'Method Selection', link: '/tutorials/method-selection' }
          ]
        }
      ],
      '/docs/': [
        {
          text: 'Core Concepts',
          collapsed: false,
          items: [
            { text: 'Documentation Overview', link: '/docs/' },
            { text: 'Correction Methods', link: '/docs/correction-methods' },
            { text: 'Single File Analysis', link: '/docs/single-file' },
            { text: 'Understanding Results', link: '/docs/understanding-results' }
          ]
        },
        {
          text: 'Processing Workflows',
          collapsed: false,
          items: [
            { text: 'Multi-File Processing', link: '/docs/multi-file' },
            { text: 'Quality Control', link: '/docs/quality-control' }
          ]
        },
        {
          text: 'Reference Materials',
          collapsed: false,
          items: [
            { text: 'Troubleshooting', link: '/docs/troubleshooting' }
          ]
        }
      ]
    },

    socialLinks: [
      { icon: 'github', link: 'https://github.com/AfriGen-D/checkref' }
    ],

    footer: {
      message: 'Released under the MIT License.',
      copyright: 'Copyright © 2025 AfriGen-D Project'
    },

    search: {
      provider: 'local'
    }
  },

  markdown: {
    lineNumbers: true
  }
})