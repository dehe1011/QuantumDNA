import os
import sys

sys.path.insert(0, os.path.abspath(".."))

import qDNA

project = "QuantumDNA"
copyright = "2024, Dennis Herb"
author = "Dennis Herb"
release = qDNA.__version__


extensions = [
    "sphinx_rtd_theme",
    "numpydoc",
    "sphinxcontrib.bibtex",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    # "nbsphinx",
    # "pandoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
]
autosummary_generate = True

# Indicate BibTex file
bibtex_bibfiles = ["biblio.bib"]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Options for HTML output
html_theme = "sphinx_rtd_theme"
html_static_path = []
