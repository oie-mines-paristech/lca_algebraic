import os
import platform
import subprocess
from datetime import datetime

from setuptools import setup

REQUIREMENTS = [
    "tabulate",
    "ipywidgets",
    "pandas",
    "pyarrow",
    "seaborn",
    "sympy",
    "matplotlib",
    "deprecation",
    "brightway2",
    "SALib",
    "pint",
    "typing-extensions",
]

PYPARDISO = "pypardiso"


# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read().strip()


def run(args):
    return subprocess.run(args, stdout=subprocess.PIPE).stdout.decode("utf-8").splitlines()


version = read("VERSION")
name = "lca_algebraic"

if platform.processor() == "x86_64":
    REQUIREMENTS.append(PYPARDISO)

# Try to get branch from git
try:
    branches = run(["git", "branch"])
    curr_branch = next(line for line in branches if "*" in line)
    curr_branch = curr_branch.replace(" ", "").replace("*", "")

    if curr_branch != "main":
        if "/" in curr_branch:
            curr_branch = curr_branch.split("/")[-1]

        name += f"_{curr_branch}"

        # commit = run(["git", "log"])[0].split()[1][0:8]

        start = datetime.strptime("2021-01-01", "%Y-%m-%d")
        now = datetime.now()

        min_diff = int((now - start).total_seconds() // 60)

        version += "." + str(min_diff) + "_dev"

except Exception as e:
    print("Failed to get git branch. Might be in TOX ?.", e)

setup(
    name=name,
    version=version,
    author="OIE - Mines ParisTech",
    author_email="raphael.jolivet@mines-paristech.fr",
    description=(
        "This library provides a layer above brightway2 for defining "
        "parametric models and running super fast LCA for monte carlo analysis."
    ),
    license="BSD",
    keywords="LCA brightway2 monte-carlo parametric",
    url="https://lca-algebraic.readthedocs.io/en/stable/",
    project_urls={
        "Homepage": "https://github.com/oie-mines-paristech/lca_algebraic",
        "Repository": "https://github.com/oie-mines-paristech/lca_algebraic.git",
        "Documentation": "https://lca-algebraic.readthedocs.io/",
        "Changelog": "https://github.com/oie-mines-paristech/lca_algebraic/blob/main/RELEASE_NOTES.md",
    },
    packages=["lca_algebraic"],
    long_description=read("README.md"),
    long_description_content_type="text/markdown",
    classifiers=[],
    install_requires=REQUIREMENTS,
)
