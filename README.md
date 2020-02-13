# Troubleshoting

* Had to downgrade to conda 3.7.5 with python3.6. see : https://github.com/conda/conda/issues/9367#issuecomment-555156949


# Notes on model 

* Recycling aluminium : why electricity is not reduced ? why updating alu instead of creating linear interpolation between scrap and new alu ?
* Mounting system : scrap steel & alu are equal to alu and steel used ? This is not the case for ground model
* "Market for scrap alu" vs "aluminium scrap, new, Recycled Content cut-off"
* Wierd ground 
* Electrical installation : why by weight ?? Function is not affine ?
* Inverter : 3Kw or 2.5kw ?? Scrap alu : Yet another ? Positive recycling ? Input not changed ?
* Silicon : elec dataset not reset
* Casting dataset : Water = Water, cooling = 5 in Scarlet model ?? 
* Wafer manufacturing : why no elec switch ?
* generate_diamond_powder : Electricy was randomly selected ('in')
* WAfer : 
    * new diamond wiring not used ?
    * TAP water is huge
    * exc['amount'] += : cumulative !!
    * DW = False / True : values stay between computation

* TEG == Silicon ??
* PV cell manufacturing : 
    * 'metallization paste production, back side, copper' does not exist : not created either
    * New copper amount 0.67 * 0.67 ??
    * Silver amount :1e-3 (not others ??)
    * 'aluminium alloy, AlMg3' changed for name of new dataset, but not is selection ?
    
 