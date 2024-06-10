# Units

*lca_algebraic* provides an optional support of automatic checks and conversion of **physical units**. It can maintain the 
unit of quantities expressed as algebraic expression of parameters, and ensure [dimensional consistency](https://en.wikipedia.org/wiki/Dimensional_analysis) of your 
 inventory.

For this purpose, *lca_algebraic* integrates the great library [Pint](https://pint.readthedocs.io/en/stable/)

## Activation

To activate the support for units, change the global `Settings` at the beginning of your notebook / script :

```python
from lca_algebraic import Settings

Settings.units_enabled = true
```

## Usage

Once units are activated, `newFlotParam(...)` creates a [Pint Quantity](https://pint.readthedocs.io/en/stable/user/defining-quantities.html?highlight=quantity), 
which holds both a **physical unit** (the one defined in the `unit` argument of `newFloatParam`) 
and a **magnitude** (the lca_algebraic parameter itself). 

You also need to specify the `unit` for any static value used in formulas. 
You can specify them using either of the following syntaxes :
* `<value> * <unit>` or 
*  `<value> | <unit>`

Example:

```python
import lca_algebraic as agb
from lca_algebraic import unit_registry as u

agb.Settings.units_enabled = True

# We define a static value, with its unit
ENERGY_PER_MASS = 1500 * u.megajoule / u.kg
# or
ENERGY_PER_MASS = 1500 | u.megajoule / u.kg

# We define a lca_algebraic parameter, with its unit
p1 = agb.newFloatParam("p1", unit="kg", default=10)

# We combine the two in an arithmetic expression
print(p1*ENERGY_PER_MASS)
```

This code should return 
```
1500*p1 megajoule
```

As you see, the physical unit of the whole expression is maintained correctly.

On the other hand, the following expression should fail. You can't add a mass to an energy per mass.
```python
print(p1+ENERGY_PER_MASS)
```
This should fail with the following exception : 
```
pint.errors.DimensionalityError: Cannot convert from 'kilogram' ([mass]) to 'megajoule / kilogram' ([length] ** 2 / [time] ** 2)
```

## Convert between units

You can convert a quantity to a compatible unit (of same physical dimension), using either the method 
`<quantity>.to(<unit>)` or the syntax `<quantity> | <unit>`

Example : 

```python
from lca_algebraic import unit_registry as u

print((2 * u.ton).to(u.kg))
print(2 * u.ton | u.kg)
```

This code should print 
```
2000.0 kilogram
2000.0 kilogram
```

## Unit registry

By default, *lca_algebraic* creates a [Unit Registry](https://pint.readthedocs.io/en/stable/api/base.html#pint.UnitRegistry), 
holding most usual units, plus the ones used in *ecoinvent* database.

Units defined in this registry can be accessed as properties (`u.<unit>`) or dictionnary keys (`u["long_unit"]`)

```python
from lca_algebraic import unit_registry as u

print(u.km)
print(u["square meter"])
print(u.kilometer ** 2 / u.kg)
print(u.kWh) # Beware of upper/lower case
```

### Adding new units

You can add extra units to the registry

**Aliases**

You may define a named alias composed of other units

```python

from lca_algebraic import unit_registry as u
import lca_algebraic as agb

agb.define_alias_unit("gigaton", 10**9 * u.ton)

# Should print '1000000000.0 metric_ton'
print(2 * u.gigaton | u.ton)
``` 

**New unit**

You may also define completely separate unit, as a new physical dimension that can't be converted to anything else.

```python

from lca_algebraic import unit_registry as u
import lca_algebraic as agb

agb.define_separate_unit("person")

# Should print "2 person"
print(2 * u.person)
``` 


## Units in exchanges

When assigning a quantity an exchange, either via `newActivity()` or `updateActivity()`, *lca_algebraic* will check that the 
unit of the quantity is either :

* The unit of the **target activity** of the exchange 
* The unit of the **target activity** of the exchange, divided by the **unit of the output**.

For instance, if you define an activity **act1**, creating **1 kg** of a product, and using some electricity (in 
*kWh*), then the exchange for electricity should be defined either in *kWh* (implicitely "for one 1 kg of product created") or in 
*kWh/kg*.

```python

from lca_algebraic import unit_registry as u
import lca_algebraic as agb

# The unit of electricity is kilowatthour
electricity = agb.findTechAct(...)

act1 = agb.newActivity(
    db_name=USER_DB,
    name="act1",
    unit="kg",
    exchanges={
        electricity: 100 | u.kWh, # Ok
        electricity: 100 | u.kWh / u.kg, # Ok
        electricity: 100 | u.kilometer # Would fail
    })
``` 

## Auto scale

*auto_scale* is a property of the unit registry. It enables the library to silently / automatically transform 
between units of the same dimensions (from *kilometer* to *meter* for instance). 

It is disabled by default, forcing the user to explicitly ask for a conversion of units (using the method `.to(..)` or the 
syntax `| <unit>`). 

Enabling the autoscale :
```python
from lca_algebraic import unit_registry as u 

u.auto_scale = True # False by default

```
We advise to leave it disabled, as it leads to better awarness of physical units to the user and may limit errors even further.

In practice, when disabled (default), you need to explicitly convert between units :

```python

from lca_algebraic import unit_registry as u
import lca_algebraic as agb

# The unit of electricity is kilowatthour
electricity = agb.findTechAct(...)

ENERGY = 105 | u.megajoule

act1 = agb.newActivity(
    db_name=USER_DB,
    name="act1",
    unit="kg",
    exchanges={
        electricity: ENERGY | u.kWh # Megajoule is compatible with kWh. Its conversion should be asked explicitly as shown here.
    })
``` 

## Full API

### check_unit_consistency()

```{eval-rst} 
.. autofunction:: lca_algebraic.check_unit_consistency
```

### define_separate_unit()

```{eval-rst} 
.. autofunction:: lca_algebraic.define_separate_unit
```

### define_alias_unit()

```{eval-rst} 
.. autofunction:: lca_algebraic.define_alias_unit
```



