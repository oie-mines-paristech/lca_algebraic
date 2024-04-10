# Introduction

This library is a layer above [**brightway2**](https://brightway.dev/) designed for the definition of **parametric inventories** 
with fast computation of LCA impacts, suitable for **monte-carlo** analyis.

It integrates the magic of [Sympy](https://www.sympy.org/en/index.html) in order 
to write parametric formulas as regular Python expressions.

**lca-algebraic** provides a set of  **helper functions** for : 
* **compact** & **human readable** definition of activites :  
    * search background (tech and biosphere) activities 
    * create new foreground activites with parametrized amounts
    * parametrize / update existing background activities (extending the class **Activity**)
* Definition of parameters
* Fast computation of LCAs
* Computation of monte carlo method and global sensitivity analysis (Sobol indices) 

# Installation

We don't provide conda package anymore.
Please use pip for installation

## 1) Setup separate environement

First create a python environment, with **Python** [>=3.9] :

**With Conda or mamba**

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

## Source code

The source code is available of our {gitref}`github </>`

# Licence & Copyright

This library has been developed by [OIE - MinesParistech](http://www.oie.mines-paristech.fr), for the project *INCER-ACV*, 
lead by [ADEME](https://www.ademe.fr/). 

It is distributed under the {gitref}`BSD License </LICENSE>`.

# Summary 

```{toctree}
---
maxdepth: 2
---
Introduction <self>
api/index
```