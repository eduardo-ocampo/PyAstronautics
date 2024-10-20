# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'PyAstronautics'
copyright = '2024, Eduardo Ocampo'
author = 'Eduardo Ocampo'
release = '0.0.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx_wagtail_theme",
    "myst_parser"
]

# https://myst-parser.readthedocs.io/en/latest/syntax/optional.html#
myst_enable_extensions = [
    "dollarmath",
    "amsmath",
    "deflist",
    "fieldlist",
    "html_admonition",
    "html_image",
    "colon_fence",
    "smartquotes",
    "replacements",
    # "linkify",
    # "strikethrough",
    "substitution",
    "tasklist",
]
myst_heading_anchors = 4

templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
html_theme = "sphinx_wagtail_theme"
html_static_path = ['_static']

html_css_files = [
    'custom.css',
]
html_theme_options = {
       "project_name": project,
    #    "default_mode": "light"
}
