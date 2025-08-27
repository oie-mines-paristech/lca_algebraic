

<img src="./doc/source/_static/img/logo_lca_algebraic.png" alt="logo" width="200" style="margin:auto;display:block"/>


![Python 3.10](https://img.shields.io/badge/python-3.10-blue) 
![Python 3.11](https://img.shields.io/badge/python-3.11-blue)
![Python 3.12](https://img.shields.io/badge/python-3.12-blue)

[![Brightway 2.4](https://img.shields.io/badge/brightway-2.4-blue)](https://docs.brightway.dev/en/legacy/index.html)
[![Brightway 2.5](https://img.shields.io/badge/brightway-2.5-blue)](https://docs.brightway.dev/en/latest/)

**lca_algebraic** is a layer above [**Brightway**](https://brightway.dev/) designed for the definition of **parametric inventories** 
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

We support both [Brightway 2.4 (legacy)](https://docs.brightway.dev/en/legacy/index.html)
and [Brightway 2.5](https://docs.brightway.dev/en/latest/) via two separate branches / libraries :

* [lca_algebraic](https://pypi.org/project/lca-algebraic/) (for Brightway 2.4)
* [lca_algebraic_bw25](https://pypi.org/project/lca-algebraic-bw25/) (for  Brightway 2.5)

## 1) Setup separate environment

First create a python environment, with **Python** [>=3.10, <3.13] :

**With Conda (or [mamba](https://mamba.readthedocs.io/en/latest/index.html))**

```bash
conda create -n lca python==3.10
conda activate lca
```

**With virtual env**

```bash
python3.10 -m venv .venv
source .venv/bin/activate
```

## 2) Install lca_algebraic

> pip install lca_algebraic 

Or, for brightway 25 :

> pip install lca_algebraic_bw25

## 3) [Optional] Install Jupyter & Activity Browser 

You may also install Jupyter and [Activity Browser](https://github.com/LCA-ActivityBrowser/activity-browser) on the same 
environment.

**Jupyter** :
> pip  install jupyter

**Activity Browser** can only be installed via conda/mamba. Note that it can also be installed on a separate Python env and will 
still be able to access and browse the projects created programmatically with *lca_algebraic* / *Brightway*.  
> conda install activity-browser

> **NOTE**
> While the inventories created in *lca_algebraic* are stored in the Brightway project, 
> the formulas and parameters are not compatible with **Activity Browser**
> Before computing impacts with vanilla **Brightway2** or **Activity Browser**, 
> you may use the function [freezeParams()](https://lca-algebraic.readthedocs.io/en/stable/api/parameters.html#lca_algebraic.freezeParams) 
> to update the amounts in your database for a given scenario / set of parameter values.     



# 📚 Documentation & resources

Full documentation is [hosted on **readthedocs**](https://lca-algebraic.readthedocs.io/)

We provide some notebooks :
* [Example notebook](./notebooks/example-notebook.ipynb) : Basic functionalities  
* [Handbook](./notebooks/handbook.ipynb) : More examples, also showing the usage of the Brightway functions.
* [Workshop](https://git.sophia.mines-paristech.fr/oie/lca-algebraic-workshop) :
  A "real life" exercise used as a short training on *lca_algebraic*

# 📧 Mailing list

Please register to this dedicated mailing list to discuss the evolutions of this library and be informed of future releases :

[lca_algebraic@groupes.mines-paristech.fr](https://groupes.minesparis.psl.eu/wws/subscribe/lca_algebraic)


# © Licence & Copyright

This library has been developed by [MinesParis - PSL - O.I.E team](https://www.oie.minesparis.psl.eu/), for the project [*INCER-ACV*](https://librairie.ademe.fr/energies-renouvelables-reseaux-et-stockage/4448-incer-acv.html), 
lead by [ADEME](https://www.ademe.fr/). 

It is distributed under the [BSD License](./LICENSE)

# Logo

Please use the following logo to advertise about this librairy.
![](./doc/source/_static/img/logo_lca_algebraic.png)
