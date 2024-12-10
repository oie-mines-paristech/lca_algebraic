# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import os
import sys

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "lca_algebraic"
copyright = "2022, Centre O.I.E - Raphaël Jolivet"
author = "Centre O.I.E - Raphaël Jolivet"

sys.path.append(os.path.join(os.path.dirname(__file__), "_deps"))
sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.autodoc", "sphinx.ext.viewcode", "myst_parser", "sphinxarg.ext", "gitrep2", "nbsphinx"]

myst_enable_extensions = ["attrs_block", "attrs_inline"]

# These folders are copied to the documentation's HTML output
html_static_path = ["_static"]

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
html_css_files = [
    "css/custom.css",
]

html_logo = "_static/img/logo_lca_algebraic.png"
logo_only = True

templates_path = ["_templates"]
exclude_patterns = []
source_suffix = [".rst", ".md"]

gitref_remote_url = "https://github.com/oie-mines-paristech/lca_algebraic.git"
gitref_branch = "main"


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = "pydata_sphinx_theme"
html_theme = "furo"
nbsphinx_execute = "never"
