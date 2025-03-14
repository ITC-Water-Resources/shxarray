# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
from datetime import datetime
import os

project = 'shxarray'
copyright = str(datetime.now().year)+', Roelof Rietbroek'
author = 'Roelof Rietbroek'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration


extensions = ['nbsphinx', 'sphinxcontrib.apidoc', 'sphinx.ext.autodoc','sphinx.ext.napoleon']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store','**.ipynb_checkpoints']

apidoc_template_dir='_templates'
#figure out the actual installation directory
import shxarray
apidoc_module_dir=os.path.dirname(shxarray.__file__)

apidoc_output_dir = 'references'
apidoc_separate_modules = True
apidoc_module_first=True
apidoc_toc_file=False

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
#html_theme = 'sphinx_rtd_theme'
html_theme='pydata_sphinx_theme'

html_static_path = ['_static']

napoleon_numpy_docstring = True

nbsphinx_prolog = """
Download this Jupyter notebook from `github <https://github.com/ITC-Water-Resources/shxarray/blob/main/docs/source/notebooks/{{ env.doc2path(env.docname, base=None) }}>`_

----
"""
html_theme_options = {
    "use_edit_page_button": False,
    "icon_links": [
        {
            # Label for this link
            "name": "GitHub",
            # URL where the link will redirect
            "url": "https://github.com/ITC-Water-Resources/shxarray",  # required
            # Icon class (if "type": "fontawesome"), or path to local image (if "type": "local")
            "icon": "fa-brands fa-github",
            # The type of image to be used (see below for details)
            "type": "fontawesome",
        }
   ],
   "logo": {
        "alt_text": "shxarray - Home",
        "image_light": "_static/shxarraylogo_dark.png",
        "image_dark": "_static/shxarraylogo_light.png",
    }
}

html_favicon= '_static/favicon.ico'

# html_context = {
    # # "github_url": "https://github.com", # or your GitHub Enterprise site
    # "github_user": "ITC-Water-Resources",
    # "github_repo": "shxarray",
    # "github_version": "main",
    # "doc_path": "docs/source",
# }


# def linkcode_resolve(domain, info):
    # if domain != 'py':
        # return None
    # if not info['module']:
        # return None
    # filename = info['module'].replace('.', '/')
    # return "https://github.com/ITC-Water-Resources/shxarray/%s.py" % filename

# linkcode_url="https://github.com/ITC-Water-Resources/shxarray"
# linkcode_link_text="view on"
