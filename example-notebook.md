<!-- #region -->
# Introduction

This notebook presents **lca-algebraic** a small libray above **brightay2**, designed for the definition of **parametric inventories** with fast computation of LCA impacts, suitable for **monte-carlo** analyis.

**lca-algebraic** provides a set of  **helper functions** for : 
* **compact** & **human readable** definition of activites :  
    * search background (tech and biosphere) activities 
    * create new foreground activites with parametrized amounts
    * parametrize / update existing background activities (extending the class **Activity**)
* Definition of parameters
* Computation of LCAs 
* Computation of statistics (including Sobols indices)
* Compute simplified parametric model by fixing minor input parameters

  
## Principles 

The main idea of this libray is to move from **procedural definition** of models (slow and prone to errors) to a **declarative / purely functionnal** definition of parametric models (models as **pure functions**). 

This enables **fast computation of LCA impacts**, useful for Monte Carlo methods and statistical analysis. 
We leverage the **power of symbolic calculus** provided by the library [SymPy](https://www.sympy.org/en/index.html).

We define our model in a **separate database**, as a nested combination of : 
* other foreground activities
* background activities :
    * Technical, refering **ecoinvent DB**
    * Biosphere, refering **brightway2** biosphere activities
    
The **amounts** in exchanges are expressed either as **static amounts**, or **symbolic expressions** of pre-defined **parameters**.

Each activity of our **root model** is defined as a **parametrized combination** of the **foreground activities**, which can themselves be expressed by the **background activities**.

When computing LCA for foreground models, the library develops the model as a combination of **only background activities**. It computes **once for all** the impact of **all required background activities** and compiles a **fast numpy** (vectorial) function for each impact, replacing each background activity by the **static value of the corresponding impact**.

By providing **large vectors** of **parameter values** to those numpy functions, we can compute LCA for **thousands of values** at a time.

![](https://oie-mines-paristech.github.io/lca_algebraic/doc/lca-algebraic.png)


## Compatiblity with brightway2 

Under the hood, the activities we define with **lca-algebraic** are standard **brightway2** activities. 
The amounts of exchanges are stored as **float values** or **serialized as string** in the property **formula**.

Parameters are also stored in the **brightay2** projets, making it fully compatible with **brightway**.

Thus, a model defined with **lca-algebraic** is stored as a regular **bw2** projet. We can use **bw2** native support for [parametrized dataset](https://2.docs.brightway.dev/intro.html#parameterized-datasets) for computing LCAs, even if much more slower than the method explain here.

## Doc

The followng notebook explores the main functions.
Full documentation of the functions is [available here](https://oie-mines-paristech.github.io/lca_algebraic/doc/)
<!-- #endregion -->

```{python}
# Those two lines are for dev only : they watch imported libraries for changes
# %load_ext autoreload
# %autoreload 2

import pandas as pd
import time
import matplotlib.pyplot as plt
import numpy as np
import brightway2 as bw

# Custom utils defined for inter-acv
from lca_algebraic import *
from lca_algebraic.stats import * 
import lca_algebraic
from sympy import init_printing, simplify

init_printing()
```

# Init brightway2 and databases

```{python}
### Init the brightway2 project :choose any project name
initProject('MyProject')

# Import Ecoinvent DB (if not already done)
# Update the name and path to the location of the ecoinvent database
importDb("ecoinvent 3.4", './ecoinvent 3.4_cutoff_ecoSpold02/datasets')

# We use a separate DB for defining our foreground model / activities
# Choose any name
USER_DB = 'Foreground DB'
```

```{python}
# This is better to cleanup the whole foreground model each time, and redefine it in the notebook
# instead of relying on a state or previous run.
# Any persistent state is prone to errors.
resetDb(USER_DB)

# Parameters are stored at project level : 
# Reset them also
# You may remove this line if you import a project and parameters from an external source (see loadParam(..))
resetParams()

# Overview of the databases
list_databases()
```

<!-- #region -->
# Introduction to Numpy


Numpy is a python libray for symbolic calculus. 

You write Sympy expression as you write **standard python expressions**, using **sympy symbols** in them. 


The result is then a **symbolic expression that can be manipulated**, instead of a **numeric value**.
<!-- #endregion -->

```{python}
from sympy import symbols 

# create sympy symbol
a = symbols("a")

# Expressions are not directly evaluated 
f = a * 2 + 4 
f
```

```{python}
# symbols can be replaced by values afterwards 
f.subs(dict(a=3))
```

In practice, you don't need to care about Sympy. Just remember that : 
* The parameters defined below are **instances of sympy symbols**
* Any **valid python expression** containing a **sympy symbol** will create a **sympy symbolic expression**


# Define input parameters

First, we define the input parameters of the model together with their distribution.

The numeric parameters are **instances of sympy 'Symbol'**. 

Thus, any python arithmetic expression composed of parameters will result in a **symbolic expression** to be used later in the definition of the model, rather than a static numeric result.

```{python}
# Example of 'float' parameters
a = newFloatParam(
    'a', 
    default=0.5, min=0.2, max=2,  distrib=DistributionType.TRIANGLE, # Distribution type, linear by default
    description="hello world",
    label="extended label for a")

b = newFloatParam(
    'b',
    default=0.5, # Fixed if no min /max provided
    description="foo bar")

# Parameters can be attached to a specific DB instead of global project
# We recommand this for team working / working with several DBs
# lca_algebraic supports resolution of parameters with duplicate names attached to different scopes 
b2 = newFloatParam(
    'b',
    default=0.5, # Fixed if no min /max provided
    description="foo bar",
    dbname=USER_DB)

share_recycled_aluminium = newFloatParam(
    'share_recycled_aluminium',  
    default=0.6, min=0, max=1, std=0.2, distrib=DistributionType.NORMAL, # Normal distrib, with std dev
    description="Share of reycled aluminium")

c = newFloatParam(
    'c',  
    default=0.6, std=0.2, distrib=DistributionType.LOGNORMAL)

beta = newFloatParam(
    'beta',  
    default=0.6, std=0.2, a=2, b=5, distrib=DistributionType.BETA)

# You can define boolean parameters, taking only discrete values 0 or 1
bool_param = newBoolParam(
    'bool_param', 
    default=1)

# Example 'enum' parameter, acting like a switch between several possibilities
# Enum parameters are not Symbol themselves
# They are a facility to represent many boolean parameters at once '<paramName>_<enumValue>' 
# and should be used with the 'newSwitchAct' method 
elec_switch_param = newEnumParam(
    'elec_switch_param', 
    values=["us", "eu"], # If provided as list, all possibilities have te same probability
    default="us", 
    description="Switch on electricty mix")

# Another example enum param
techno_param = newEnumParam(
    'techno_param', 
    values={
        "technoA":0.4, 
        "technoB":0.1,
        "technoC":0.5}, # You can provide a statistical weight for each value, by using a dict
    default="technoA", 
    description="Choice of techonoly")
```

```{python}
# List of parameters
list_parameters()
```

## Persistance of parameters

By default, new parameters are kept in memory but also persisted in the project (unless save=False).
You can persist parameters afterwards with `persistParams()`.
You can load also load parameters from an existing database with `loadParams()`.
The persistance of parameters and the distribution is compatible with **Brightway2** and **Activity Browser**  [see documentation of stat_arrays](https://stats-arrays.readthedocs.io/en/latest/)

```{python}
# Load parameters previously  persisted in the dabatase.
loadParams()
```

# Manage several Datasets

lca_algebraic supports several foreground / background datasets. Background datasets are considered static / non parametrized by the library : they use standard LCA method of **Brightway2**. 

Foreground databases are considered parametric and their activities are developped as functions of parameters and background activities.

## Set status of a database

The functions **setForeground(...)** and **setBackground(...)** change the status of a database.

```{python}
setBackground(USER_DB)
list_databases()
```

```{python}
setForeground(USER_DB)
list_databases()
```

## Import / export

`lca_algebraic` extends [BW2Package](https://2.docs.brightway.dev/technical/bw2io.html), adding persistence of parameters.

```{python}
# Save database and parameters as Bzipped JSON
export_db(USER_DB, "db.bw2") 
```

```{python}
# Reimport DB
import_db("db.bw2")
```

## Freeze 

A foreground database can be "frozen" to be used as a background database for a specific scenario : the parametrized amounts in the exhanges are computed for a given configuration of the parameters, and replaced by their value. The formulas are still stored in the database and not lost : the database can still be used as a foreground database until its status is changed with `setBackground(...)`.

This feature is useful for studies requiring several datasets to be used as background by other ones. It also enables to use standard Brightway2 tools, not aware of parametrization. 


```{python}
freezeParams(
    USER_DB, # Name of database to freeze
    
    a=1, b=2) # custom parameter values
```

# Get references to background activities

We provide two functions for easy and fast (indexed) search of activities in reference databases : 
* **findBioAct** : Search activity in **biosphere3** db
* **findTechAct** : Search activity in **ecoinvent** db

Those methods are **faster** and **safer** than using traditionnal "list-comprehension" search : 
They will **fail with an error** if **more than one activity** matches, preventing the model to be based on a random selection of one activity.


```{python}
# Biosphere activities
ground_occupuation = findBioAct('Occupation, industrial area') # Search by name
heat = findBioAct('Heat, waste', categories=['air']) # Add category selector

# Technosphere activities

# You can add an optionnal location selector
alu = findTechAct("aluminium alloy production, AlMg3", loc="RER")
alu_scrap = findTechAct('aluminium scrap, new, Recycled Content cut-off')

# Elec 
eu_elec = findTechAct("market group for electricity, medium voltage", 'ENTSO-E')
us_elec = findTechAct("market group for electricity, medium voltage", 'US')
```

# Define the model

The model is defined as a nested combination of background activities with amounts.

Amounts are defined either as constant float values or algebric formulas implying the parameters defined above.


## Create new activities

```{python}
# Create a new activity
activity1 = newActivity(USER_DB, # We define foreground activities in our own DB
    "first foreground activity", # Name of the activity
    "kg", # Unit
    exchanges= { # We define exhanges as a dictionarry of 'activity : amount'
        ground_occupuation:3 * b, # Amount can be a fixed value 
        heat: b + 0.2  # Amount can be a Sympy expression (any arithmetic expression of Parameters)
    })

# You can create a virtual "switch" activity combining several activities with a switch parameter
elec_mix = newSwitchAct(USER_DB, 
    "elect mix", # Name
    elec_switch_param, # Sith parameter
    { # Dictionnary of enum values / activities
        "us" : us_elec, # By default associated amount is 1
        "eu" : (eu_elec, 0.8)  # You can also provide custom amout or formula with a tuple 
    })
```

## Copy and update existing activity

You can copy and update an existing background activity.

Several new helper methods have been added to the class **Activity** for easy update of exchanges.

```{python}
alu2 = copyActivity(
    USER_DB, # The copy of a background activity is done in our own DB, so that we can safely update it                
    alu, # Initial activity : won't be altered
    "Aluminium 2") # New name

# Update exchanges by their name 
alu2.updateExchanges({
    
    # Update amount : the special symbol *old_amount* references the previous amount of this exchange
    "aluminium, cast alloy": old_amount * (1 - share_recycled_aluminium),
    
    # Update input activity. Note also that you can use '*' wildcard in exchange name
    "electricity*": elec_mix,
    
    # Update both input activity and amount. 
    # Note that you can use '#' for specifying the location of exchange (useful for duplicate exchange names)
    "chromium#GLO" : dict(amount=4.0, input=activity1)
}) 

# Add exchanges 
alu2.addExchanges({alu_scrap :  12})
```

## Final model

Usually, we normalize the final model as the whole LCI divided by the functional value (production of the system)


```{python}
functional_value = a + 5

model = newActivity(USER_DB, "model", "kg", {
    activity1 : b * 5 + a + 1, # Reference the activity we just created
    alu2: 3 * share_recycled_aluminium, 
    alu:0.4 * a})

normalized_model = newActivity(USER_DB, "normalized model", "kg", {
    model : 1 / functional_value})
```

## Or load existing model /activities from database

Alternatively, you may not define the model again, but load it from the USER DB.
We don't recommand it while developping the model as it may lead to unstable code / uncontrolled state of the DB.

```{python}
activity1 = findActivity("first foreground activity", db_name=USER_DB)
model = findActivity("model", db_name=USER_DB)
normalized_model = findActivity("normalized model", db_name=USER_DB)
alu2 = findActivity("Aluminium 2", db_name=USER_DB)
```

## Display activities

**printAct** displays the list of all exchanges of an activity.

Note that symbolic expressions have not been evaluated at this stage

```{python}
# Print_act displays activities as tables
printAct(activity1) 
printAct(model)
printAct(normalized_model)
```

```{python}
# You can also compute amounts by replacing parameters with a float value 
printAct(activity1, b=1.5) 
```

```{python}
# You print several activities at once to compare them
printAct(alu, alu2)
```

# Select the impacts to consider

```{python}
# List of impacts to consider
impacts = [m for m in bw.methods if 'ILCD 1.0.8 2016' in str(m) and 'no LT' in str(m)]
```

```{python}
# You can provide alternate names for some methods
set_custom_impact_labels({
    impacts[0] : 'Resource usage',
    impacts[1]: 'Climate change'})
```

# Compute LCA

We provide two methods  for computing LCA : 
* **multiLCA** : It uses **brightway2** native parametric support. It is **much slower** and kept for **comparing results**.
* **multiLCAAlgebric** : It computes an algebric expression of the model and computes LCA once for all the background activities. Then it express each impact as a function of the parameters. This expression is then compiled into 'numpy' native code, for fast computation on vectors of samples. This version is 1 million time faster.

```{python}
# Uses brightway2 parameters
multiLCA(
    normalized_model, 
    impacts, 
                   
    # Parameters of the model
    a=1, 
    elec_switch_param="us",
    share_recycled_aluminium=0.4)
```

```{python}
# Compute with algebric implementation 
# The values should be the same as above
multiLCAAlgebric(
    normalized_model, # The model 
    impacts, # Impacts
    
    # Parameters of the model
    a=1, 
    elec_switch_param="us",
    share_recycled_aluminium=0.4)
```

```{python}
# You can compute several LCAs at a time and compare them:
multiLCAAlgebric(
    [alu, alu2], # The models
    
    impacts, # Impacts
    
    # Parameters of the model
    share_recycled_aluminium=0.3,
    elec_switch_param="us")
```

## Evaluate the contribution of a subset of activites to global impact

You can extract the contribution of a subset of activites to the global impact of the model

```{python}
# Contribution of impact of aluminium
multiLCAAlgebric(
    normalized_model, # The model 
    impacts, # Impacts
    
    # List of sub activites to consider
    extract_activities = [alu2],
    
    # Parameters of the model
    a=1, 
    elec_switch_param="us",
    share_recycled_aluminium=0.4)
```

## Fast computation of many parameter values

```{python}
# Fast computation for millions of separate samples
multiLCAAlgebric(
    model, # The model 
    impacts, # Impacts
    
    # Parameters of the model
    a=list(range(1, 100000)), # All lists should have the same size
    share_recycled_aluminium=1, # Those parameters are fixed
    elec_switch_param="eu")
```

 # Statistic functions
 
 ## One at a time 
 
 We provide several functions for computing **statistics** for **local variations** of parameters (one at a time).
 
 ### oat_matrix(model, impacts)
 
 Shows a **matrix of impacts x parameters** colored according to the variation of the impact in the bounds of the parameter.

 


```{python}
oat_matrix(model, impacts)
```

### oat_dashboard_matrix

This functions draws a dashboard showing :
* A dropdown list, for choosing a parameter
* Several graphs of evolution of impacts for this parameter
* Full table of data
* A graph of "bars" representing the variation of each impact for this parameter (similar to the information given in oat_matrix) 

```{python}
oat_dashboard_interact(
    model, impacts, 
    
    # Optionnal layout parameters
    figspace=(0.5,0.5),
    figsize=(15, 15),
    sharex=True)
```

<!-- #region -->
## Monte-carlo methods & Sobol indices

Here we leverage fast computation of monte-carlo approches. 

We compute **global sensivity analysis** (GSA).
Not only local ones.


### Sobol Matrix 

Similar to OAT matrix, we compute Sobol indices. they represent the ratio between the variance due to a given parameter and the total variance.

for easier comparison, we translate those relative sobol indices into "deviation / mean" importance :

$$RelativeDeviation = \frac{\sqrt{sobol(param) \times totalVariance(impact))}}{mean(impact)}$$



<!-- #endregion -->

```{python}
# Show sobol indices 
incer_stochastic_matrix(model, impacts)
```

###  Graphs of impacts and their distribution

We provide a dashboard showing **violin graphs** : the exact probabilistic distribution for each impact. Together with medians of the impacts.

```{python}
incer_stochastic_violin(
    model, impacts,
    
    # Optionnal layout parameters
    figspace=(0.5,0.5),
    figsize=(15, 15),
    sharex=True, 
    nb_cols=3)
```

```{python}
##### Alternatively, graphs can be shown horizontally, together with a box of statistical outcomes
distrib(
    model, impacts,
    
    # Optionnal layout parameters
    height=7, width=15,
    nb_cols=2,
    percentiles=[5, 95])
```

### Full dashboard

A dashboard groups all this information in a single interface with tabs.

It also shows total variation of impacts. This last graph could be improved by showing stacked colored bars with the contribution of each parameter to this variation, according to Sobol indices. 

```{python}
incer_stochastic_dashboard(model, impacts)
```

# Producing simplified models 

One of te outcome of the statisticall analysis above would be to identify main input parameters and produce simplidied models, fixing the minor ones.

We provide several functions for doing this.

## Explore initial algrebraic model

```{python}
# First, let's look at the full expression defining our model
expr, _ = actToExpression(normalized_model)
expr
```

## Compute simplified models

We provide some method to automatically select a subset of parameters, based on the **sobol indices**, and then compute simplified models for it.

We also round numerical expression to 3 digits, and we remove terms in sums that are less than 1% of total.

```{python}
simplified = sobol_simplify_model(
    normalized_model, # The model
    impacts, # Impacts to consider
    n=10000, # For large model, you may test other value and ensure ST and sum(S1) are close to 1.0 
    fixed_mode = FixedParamMode.MEDIAN, # We replace minor parameters by median by default,
    min_ratio=0.8, # Min ratio of variability to explain

    num_digits=3)
```

```{python}
# Let's look at the expression for first impact again 
# much simpler ! 
simplified[0].expr
```

## Compare simplified model with full model

Finally, we can compare the distribution of those simplified model against the full model. We provide a function for graphical display of it, and compuation of de R-Square score.


```{python}
compare_simplified(normalized_model, impacts, simplified)
```
