# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import pygaia

# -- Project information -----------------------------------------------------

project = "PyGaia"
copyright = "2012-2023, Anthony Brown, Leiden University, Gaia Data Processing and Analysis Consortium"
author = "Anthony Brown"

# The full version, including alpha/beta/rc tags
# release = "3.0"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx_automodapi.automodapi",
    "matplotlib.sphinxext.plot_directive",
    "sphinx.ext.napoleon",
]
numpydoc_show_class_members = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"
html_title = "PyGaia"
# Not needed for now
# html_static_path = ["_static"]
html_show_sourcelink = False
html_theme_options = {
    "path_to_docs": "docs",
    "repository_url": "https://github.com/agabrown/PyGaia",
    "repository_branch": "main",
    "launch_buttons": {
        "binderhub_url": "https://mybinder.org",
        "colab_url": "https://colab.research.google.com/",
        "notebook_interface": "jupyterlab",
    },
    "use_edit_page_button": False,
    "use_issues_button": False,
    "use_repository_button": True,
    "use_download_button": True,
}
html_baseurl = "https://pygaia.readthedocs.io/en/latest/"
