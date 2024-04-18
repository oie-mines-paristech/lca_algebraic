# Introduction

This library is a layer above [**brightway2**](https://brightway.dev/) designed for the definition of **parametric inventories** 
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

# Installation

We don't provide conda package anymore.

This packages is available via [pip /pypi](https://pypi.org/project/lca-algebraic/)

## 1) Setup separate environement

First create a python environment, with **Python** [>=3.9] :

**With Conda (or [mamba](https://mamba.readthedocs.io/en/latest/index.html))**

```bash
conda env create -n lca python==3.10
conda activate lca
```

**With virtual env**

```bash
python3.10 -m venv .venv
source .venv/bin/activate
```

## 2) Install lca_algebraic

> pip install lca_algebraic


# Licence & Copyright

This library has been developed by [OIE - MinesParistech](http://www.oie.mines-paristech.fr), for the project [*INCER-ACV*](https://librairie.ademe.fr/energies-renouvelables-reseaux-et-stockage/4448-incer-acv.html), 
lead by [ADEME](https://www.ademe.fr/). 

It is distributed under the [BSD License](./LICENSE)

# Documentation

Full documentation and example notebooks are [hosted on **readthedocs**](https://lca-algebraic.readthedocs.io/)