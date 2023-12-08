# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
from datetime import datetime

project = 'shxarray'
copyright = str(datetime.now().year)+', Roelof Rietbroek'
author = 'Roelof Rietbroek'

release = 'alpha'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['nbsphinx']
# extensions = [
    # 'sphinxcontrib.apidoc',
    # 'sphinx.ext.viewcode',
    # 'sphinx.ext.githubpages',
    # 'sphinxcontrib.autoprogram',
    # 'nbsphinx','sphinx.ext.napoleon',
# ]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store','**.ipynb_checkpoints']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
