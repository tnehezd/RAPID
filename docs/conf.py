import os
import sys

# Ensure these paths are correct relative to your conf.py (which is in docs/source/)
# The Doxygen XML output is expected in docs/xml/.
# So, from docs/source/, you go up one level (to docs/) and then into 'xml'.
sys.path.insert(0, os.path.abspath('.'))      # Adds 'docs/source/' to Python path
sys.path.insert(0, os.path.abspath('../xml')) # <--- CORRECTED: This adds 'docs/xml/' to Python path for Breathe

# You don't usually need to add the main 'src' folder to sys.path for C/C++ Doxygen with Breathe,
# as Breathe reads the XML, not the raw C files.
# sys.path.insert(0, os.path.abspath('../../src')) # Remove or comment out if not needed for other Sphinx extensions


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'RAPID code'
copyright = '2025, D. Tarczay-Nehez'
author = 'D. Tarczay-Nehez'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Ensure 'breathe' and 'sphinx_rtd_theme' are the only extensions if you only need Doxygen integration
extensions = [
    'sphinx_rtd_theme',
    'breathe',
    # 'sphinx.ext.autodoc', # Primary for Python codes - keep commented or remove for C projects
]

templates_path = ['_templates']
exclude_patterns = []

# If you chose separate source and build directories during sphinx-quickstart,
# '_static' is usually created inside 'source/'.
# So, html_static_path = ['_static'] is likely correct if the folder is at docs/source/_static.
# If you don't use custom static files, you can keep it commented out.
# For now, let's keep it commented if you don't have custom static assets.
# html_static_path = ['_static']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'


# -- Breathe configuration ---------------------------------------------------

# Tell Breathe where to find the Doxygen XML output
# 'xml' here is the directory name where Doxygen put its XML files.
# This path is relative to the directory *added to sys.path* (in this case, 'docs/xml/').
# So, if sys.path.insert(0, os.path.abspath('../xml')) added 'docs/xml/' to the path,
# then 'xml' here would refer to 'docs/xml/xml', which is wrong.
# The correct way is just to reference the project by its name after adding its path to sys.path.
# HOWEVER, a more direct way for breathe_projects is to give a path relative to 'conf.py'.
breathe_projects = {
    "rapid": "xml" # <--- CORRECTED: This means from 'docs/source/' go up to 'docs/' then into 'xml/'
}
breathe_default_project = "rapid"