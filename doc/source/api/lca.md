# Impacts

Functions for computing environmental impacts.

## Ventilating impacts by axis

Most impact computations return a single total per model. You often want to know
*where* the impact comes from : which life cycle phase, which sub-module, which
component. This breakdown is called an **axis**, and is enabled by the `axis`
parameter of `compute_impacts` / `compute_impacts_xarray`.

An axis is simply the name of a custom attribute (e.g. `"phase"`, `"subsystem"`)
that you set on your foreground model. You can set it at two levels, which can be
freely combined on the same axis.

### Tagging a whole activity

Tagging an activity attributes **all its outgoing exchanges** to the category. Set
the attribute with `updateMeta()`, or directly when creating / copying the activity :

```python
pv_manufacturing = agb.newActivity(USER_DB, "pv manufacturing", "unit",
                                   {pv_panel: 1},
                                   phase="manufacturing")
# or, on an existing activity :
pv_manufacturing.updateMeta(phase="manufacturing")
```

### Tagging a single exchange

When one activity mixes several categories, tag the individual exchanges instead.
Pass a **dict of attributes** (rather than a bare amount) to `newActivity`,
`addExchanges` or `updateExchanges`. Only that exchange is attributed to the
category ; the other exchanges of the same activity stay unaffected :

```python
model = agb.newActivity(USER_DB, "model", "unit", {
    pv_manufacturing: dict(amount=1, phase="manufacturing"),
    pv_operation:     dict(amount=1, phase="operation"),
    transport:        1,  # not tagged -> counted as "*other*"
})
```

### Reading the results

```python
compute_impacts(model, [climate_change], axis="phase")
# returns one line per phase, plus "*all*" (total) and "*other*" (untagged part)
```

Activities and exchanges left without the attribute are gathered in the `*other*`
line, and `*all*` always holds the total impact.

## Compute impacts

```{eval-rst} 
.. autofunction:: lca_algebraic.compute_impacts
```

## Compute inventory

```{eval-rst} 
.. autofunction:: lca_algebraic.compute_inventory
```













