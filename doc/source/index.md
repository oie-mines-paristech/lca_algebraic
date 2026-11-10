# Introduction

![Python 3.10](https://img.shields.io/badge/python-3.10-blue) 
![Python 3.11](https://img.shields.io/badge/python-3.11-blue)
![Python 3.12](https://img.shields.io/badge/python-3.12-blue)
![Python 3.13](https://img.shields.io/badge/python-3.13-blue)
![Python 3.14](https://img.shields.io/badge/python-3.14-blue)

[![Brightway 2.4](https://img.shields.io/badge/brightway-2.4-blue)](https://docs.brightway.dev/en/legacy/index.html)
[![Brightway 2.5](https://img.shields.io/badge/brightway-2.5-blue)](https://docs.brightway.dev/en/latest/)


This library is a layer above [**Brightway**](https://brightway.dev/) designed for the definition of **parametric inventories** 
with fast computation of LCA impacts, suitable for **monte-carlo** / global sensitivity analysis 

It integrates the magic of [Sympy](https://www.sympy.org/en/index.html) in order to write parametric formulas as regular Python expressions.

**lca-algebraic** provides a set of **helper functions** for : 
* **compact** & **human readable** definition of activities :  
    * search background (tech and biosphere) activities 
    * create new foreground activities with parametrized amounts
    * parametrize / update existing background activities (extending the class **Activity**)
* Definition of parameters
* Fast computation of LCAs
* Computation of monte carlo method and global sensitivity analysis (Sobol indices) 
* Support for automatic check of [homogeneity of physical units](https://lca-algebraic.readthedocs.io/en/stable/api/units.html)


# ⚙ Installation

*lca_algebraic* supports both [Brightway 2.4 (legacy)](https://docs.brightway.dev/en/legacy/index.html)
and [Brightway 2.5](https://docs.brightway.dev/en/latest/) through optional pip extras, from a single package.

## 1) Setup separate environment

First create a python environment, with **Python** [>=3.10, <3.15] :

**With Conda (or [mamba](https://mamba.readthedocs.io/en/latest/index.html))**

```bash
conda create -n lca python==3.12
conda activate lca
```

**With virtual env**

```bash
python3.12 -m venv .venv
source .venv/bin/activate
```

## 2) Install lca_algebraic

Install for *Brightway 2.4* :

```bash
pip install lca_algebraic[bw2]
```

Or for *Brightway 2.5* :

```bash
pip install lca_algebraic[bw25]
```

For development (tests, pre-commit, etc.) :

```bash
pip install lca_algebraic[bw2,dev]
# or
pip install lca_algebraic[bw25,dev]
```

## 3) [Optional] Install Jupyter & Activity Browser 

You may also install Jupyter and [Activity Browser](https://github.com/LCA-ActivityBrowser/activity-browser) on the same 
environment.

**Jupyter** :
> pip  install jupyter

**Activity Browser** can only be installed via conda/mamba. Note that it can also be installed on a separate Python env and will 
still be able to access and browse the projects created programmatically with *lca_algebraic* / *Brightway*.  
> conda install activity-browser


>  **NOTE**
> While the inventories created in *lca_algebraic* are stored in the Brightway project, 
> the formulas and parameters are not compatible with **Activity Browser**
> Before computing impacts with vanilla **Brightway2** or **Activity Browser**, 
> you may use the function [freezeParams()](https://lca-algebraic.readthedocs.io/en/stable/api/parameters.html#lca_algebraic.freezeParams) 
> to update the amounts in your database for a given scenario / set of parameter values.     


# Mailing list

Please register to this dedicated mailing list to discuss the evolutions of this library and be informed of future releases :

[lca_algebraic@groupes.mines-paristech.fr](https://groupes.minesparis.psl.eu/wws/subscribe/lca_algebraic)

# Source code

The source code is available of our {gitref}`github </>`

The list of changes for each release is available in the {gitref}`release notes </RELEASE_NOTES.md>`.


# Licence & Copyright

This library has been developed by [OIE - MinesParistech](http://www.oie.mines-paristech.fr), for the project [*INCER-ACV*](https://librairie.ademe.fr/energies-renouvelables-reseaux-et-stockage/4448-incer-acv.html), 
lead by [ADEME](https://www.ademe.fr/). 

It is distributed under the {gitref}`BSD License </LICENSE>`.


# Summary 

```{toctree}
---
maxdepth: 2
---
Introduction <self>
api/index
notebooks/index
```