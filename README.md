
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


# Usage & documentation 

Please refer to the [online documentation](https://lca_algebraic.readthedocs.org/)

# Licence & Copyright

This library has been developed by [OIE - MinesParistech](http://www.oie.mines-paristech.fr), for the project *INCER-ACV*, 
lead by [ADEME](https://www.ademe.fr/). 

It is distributed under the **[BSD licence](/LICENSE)**.

  
# Principles 

The main idea of this libray is to move from **procedural definition** of models (slow and prone to errors) to a **declarative / purely functionnal** definition of parametric models (models as **pure functions**). 

This enables **fast computation of LCA impacts**. 
We leverage the **power of symbolic calculus** provided by the great libary [SymPy](https://www.sympy.org/en/index.html).

We define our model in a **separate DB**, as a nested combination of : 
* other foreground activities
* background activities :
    * Technical, refering **ecoinvent DB**
    * Biopshere, refering **brightway2** biosphere activities
    
The **amounts** in exchanges are expressed either as **static amounts**, or **symbolic expressions** of pre-defined **parameters**.

Each activity of our **root model** is defined as a **parametrized combination** of the **foreground activities**, which can themselves be expressed by the **background activities**.

When computing LCA for foreground models, the library develops the model as a combination of **only background activities**. It computes **once for all** the impact of **background activities** and compiles a **fast numpy** (vectorial) function for each impact, replacing each background activity by the **static value of the corresponding impact**.

By providing **large vectors** of **parameter values** to those numpy functions, we can compute LCA for **thousands of values** at a time.


# Compatibility with brightway2 

Under the hood, the activities we define with **lca-algebraic** are standard **brightway2** activities. 
The amounts of exchanges are stored as **float values** or **serialized as string** in the property **formula**.

Parameters are also stored in the **brightay2** projets, making it fully compatible with **brightway**.

Thus, a model defined with **lca-algebraic** is stored as a regular **bw2** projet. We can use **bw2** native support for [parametrized dataset](https://2.docs.brightway.dev/intro.html#parameterized-datasets) for computing LCAs, even if much more slower than the method explain here.
