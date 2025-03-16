# Configuration file for the Sphinx documentation builder

import os
import sys
sys.path.insert(0, os.path.abspath('..'))  # Path to the root of your project

# Project information
project = 'HistoMapTx'
copyright = '2025, Your Name'
author = 'Your Name'

# Add extensions
extensions = [
    'sphinx.ext.autodoc',  # Include documentation from docstrings
    'sphinx.ext.viewcode',  # Add links to source code
    'sphinx.ext.napoleon',  # Support for NumPy and Google style docstrings
    'sphinx_autodoc_typehints',  # Use type hints for docs
    'nbsphinx',  # For including Jupyter notebooks
    'sphinx_gallery.gen_gallery',  # For code examples
]

# Theme
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Include Python files as both source files and executables
nbsphinx_execute = 'always'

# Syntax highlighting
pygments_style = 'sphinx'

# Include private members
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'private-members': False,
}