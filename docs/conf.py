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
    "sphinx.ext.autodoc",
    "nbsphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
]
autosummary_generate = True

# Indicate BibTex file
bibtex_bibfiles = ["biblio.bib"]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Options for HTML output
html_theme = "sphinx_rtd_theme"
html_static_path = []

exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
]

# Syntax highlighing for code blocks in the documentation
pygments_style = "sphinx"

# Do not add parentheses to function names in cross-references like :func:
add_function_parentheses = False

# Do not add full module names
add_module_names = False

numpydoc_show_inherited_class_members = False  # Show inherited class members in the documentation
numpydoc_show_property_with_doc = (
    False  # Show properties with documentation in the class documentation
)

autodoc_default_options = {
    "members": True,  # Include all members of the class/module
    "undoc-members": False,  # Include undocumented members
    "inherited-members": False,  # Include inherited members
    "show-inheritance": True,  # Show inheritance in the documentation
}
