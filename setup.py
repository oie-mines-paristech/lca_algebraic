import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "lca_algebraic",
    version = read("VERSION").strip(),
    author = "OIE - Mines ParisTech",
    author_email = "raphael.jolivet@mines-paristech.fr",
    description = ("This library provides a layer above brightway2 for defining parametric models and running super fast LCA for monte carlo analysis."),
    license = "BSD",
    keywords = "LCA brightway2 monte-carlo parametric",
    url = "https://github.com/oie-mines-paristech/lca_algebraic/",
    packages=['lca_algebraic'],
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    classifiers=[],
    install_requires=[
        'tabulate',
        'ipywidgets',
        'pandas',
        'seaborn',
        'sympy',
        'matplotlib',
        'brightway2>=2.3',
        'SALib']
)
