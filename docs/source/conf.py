import sys
import os


sys.path.insert(0, os.path.abspath('../..'))

project = 'mhodlr'
copyright = '2024, inEXASCALE'
author = 'Erin Carson, Xinye Chen, and Xiaobo Liu'

# The full version, including alpha/beta/rc tags
release = '0.0.1'


# -- General configuration ---------------------------------------------------

extensions = [
    "sphinxcontrib.matlab",
    "sphinx.ext.autodoc",
    "sphinx_copybutton",
    "myst_parser",
]
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']
source_suffix = '.rst'
exclude_patterns = []
locale_dirs = ['locale/']
gettext_compact = False

exclude_patterns = []

pygments_style = "sphinx"

html_static_path = ['_static']
html_style = 'css/_.css'
html_theme = "sphinx_rtd_theme" # html_theme = 'alabaster'
html_logo = "LOGO2.png"
html_theme_options = {
    'logo_only': False,
    'navigation_depth': 5,
}
