# 1.1 

* Fixed bug in ActivityExtended#getExchanges introduced by merge #14
* Added linear interpolation between activities with function `interpolate_activities`
* Rename multiLCaAlgebraic => compute_impacts
* Added breakdown of impacts by arbitrary attribute with the parameter `axis` of `compute_impacts`
* Added `functional_unit` in compute_impacts : You are not obliged to define a custom activity anymore 
* findActivities() is now case insensitive by default
* Fixed bug #38 : getOutputAmount was wrong with activities having circular input exchanges with themselves
* Fixed bug #37 : Increased number of results in findActivities
* Added disk cache of LCIA results and act Expressions
* Fixed bug when loading params with undefined uncertainty
* Added "formula" in params, to have parameters automatically computed from others

# 1.0.4 

Fix typo in helpers

# 1.0.3

* Custom fixes for lca_publish : 
  * added meta data "interited_from", for Activities copied from other Db
  * added meta data "isProxy", for dummy tech activities used as proxy to bio activities by lca_algebraic
* Fixed [bug #12](https://github.com/oie-mines-paristech/lca_algebraic/issues/12) : 
  functions decorated with @with_db_context now accept named parameters again
* Fixed [bug #13](https://github.com/oie-mines-paristech/lca_algebraic/issues/13) : '<' not supported between instances of 'str' and 'NoneType'
* Added support for dict of models in multLCAAlgebraic

# 1.0.2 

Disabled parallel import of DB : was causing error on windows

# 1.0.1

* Added better formatting of formulas : render enum switch as piecewize functions
* Fixed version of brightway to prevent dependencies issues

# 1.0.0

* Add multi-dataset support, replacing **SET_USER_DB** by **setBackground()**, **setForeground()**
* Added **freezeParams()** to freeze parametric dataset for specific scenario to be used by non parametric tools or other foreground datasets.
* Added dataset parameters : parameters scoped to a specific dataset, with scope / conflict resolution of parameters with same names.
* Added extended import / export of datasets with their parameters : **export_db()**, **import_db()**
* Added non regression tests

# 0.0.15

* Added explicit error when cycle is detected in activities. 
  This usually means a Bg activity has been imported into User DB :
  lca_algebraic does not support cycles / recursion in User Db, 
  since those are developped as litteral formulas. 
* Added support for several background databases. 
  User Db is guessed or can be specified using **SET_USER_DB**()
* Support coproducts in graphtraversal. 
  See [PR #9] by @ntropy-esa
* Add support for loading / saving parameters into Brightway DB, with compatibility with Brightway2 and Activity Browser.
  See functions **loadParams**() and **persistParams()**. See also issue #8.
* Added new function **exploreImpacts()**, like **printAct()** but showing impact of each exchange.  

# 0.0.14

* Internal change for generated webapp 
* Add several display options to graphs()
* Fix overlapping graph bug

# 0.0.13 

* Fixed bug in param **extract_activities** of **multiLCAAlgebric**

# 0.0.12

* Added a step of simplification of models in **sobol_simplify_model** by removing minor terms of sums
* Added support for alternate names for impact methods with the function **set_custom_impact_labels**
* Added parameter **extract_activities** to **multiLCAAlgebric** in order to compute contribution of a subset of activities to the global impact 

# 0.0.11

Added **BETA** distribution type. See [online doc](https://oie-mines-paristech.github.io/lca_algebraic/doc/params.html#lca_algebraic.params.DistributionType)

# 0.0.10

Fixed [issue #4](https://github.com/oie-mines-paristech/lca_algebraic/issues/4) : Added back handling of default params

# 0.0.9 

* Removed cache for **find{Bio/Tech}Act**, in order to lower memory usage.
Instead, we now use core Brightway capabilities for searching activities.
* Added function **sobol_simplify_model** for generating simplified models, 
and **compare_simplified** for comparing their distributions with the initial model.

# 0.0.8

Cleanup dependencies

# 0.0.7

* Added parameters 'var_params' for all stats tools, for providing list of parameters to vary.
  Other parameters will be set at their default value.
* Fixed typo in the function name **incer_stochastic_dahboard** => **incer_stochastic_dashboard**
* Fixed bug of int overload : convert everything to float

# 0.0.6 

* Rename **modelToExpr()** to **simplifiedModel()**
* Improve performance of computation. Added multithread computation for LCA and Sobol.
* Added method **sobol_simplify_model** that computes list of most meaningful parameters, 
based on Sobol analysis, and then generates simplified models based on it.

# 0.0.5

* Add function `modelToExpr(model, impacts)`, returning expression with background activities replaced 
by the value of each impact.


# 0.0.4

* Added button for **exporting** data as csv on any **Dataframe** in stats tools
* **actToExpression** now replaces fixed params by their default values, providing simplified models.
* **newSwitchAct** can take custom amounts (not only '1') by providing tuple (act, amount)

# 0.0.3

First released version