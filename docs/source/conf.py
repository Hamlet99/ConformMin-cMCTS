# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os, sys

sys.path.insert('../..')

project = 'ConformMin-cMCTS'
copyright = '2024, Hamlet Khachatryan'
author = 'Hamlet Khachatryan'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 'sphinx_mdinclude', 'sphinx_math_dollar', 'sphinx.ext.mathjax']

mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@2/MathJax.js?config=TeX-AMS-MML_HTMLorMML"

mathjax_config = {
    'tex3jax': {
        'inlineMath': [ ["\\(","\\)"] ],
        'displayMath': [["\\[","\\]"] ],
    },
}

MathJax = {
  'tex': {
    'inlineMath': [['$', '$'], ['\\(', '\\)']]
  },
  'svg': {
    'fontCache': 'global'
  }
};

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'classic'
html_static_path = ['_static']
html_theme_options = {
    "rightsidebar": "false",
    "sidebarwidth": "20%",
    "stickysidebar": "true"}
