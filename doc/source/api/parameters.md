# Parameters

Parameters extend *Sympy* **symbols**. They can be used in any python expression and result 
as Sympy Expressions to be used as amounts in **Exchanges**. 

## Create parameters

### Float parameters

A float parameter represents a decimal value. 
The user should provide its *name*, *unit* and *statistical distribution*.

```{eval-rst} 
.. autofunction:: lca_algebraic.newFloatParam
```

#### Distribution types

The following types of float distributions are supported. 

```{eval-rst} 
.. autoclass:: lca_algebraic.DistributionType
   :members:
```

### Boolean parameters

Boolean parameters can take only two values *0* and *1*.

```{eval-rst} 
.. autofunction:: lca_algebraic.newBoolParam
```

### Enum parameters

Enum parameters represent a set of mutually exclusive choices.

```{eval-rst} 
.. autofunction:: lca_algebraic.newEnumParam
```


## Utils

### List of parameters

```{eval-rst} 
.. autofunction:: lca_algebraic.list_parameters
```

```{eval-rst} 
.. autofunction:: lca_algebraic.all_params
```

### Save / load params

```{eval-rst} 
.. autofunction:: lca_algebraic.resetParams
```

```{eval-rst} 
.. autofunction:: lca_algebraic.freezeParams
```

```{eval-rst} 
.. autofunction:: lca_algebraic.loadParams
```

```{eval-rst} 
.. autofunction:: lca_algebraic.persistParams
```

### Misc

```{eval-rst} 
.. autofunction:: lca_algebraic.switchValue
```




