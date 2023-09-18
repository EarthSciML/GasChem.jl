

""" 
Arrhenius equation:
``` math
    k = a0 * exp( c0 / T ) * (T/300)**b0
```
"""
function arrhenius(T, a0, b0, c0)
    @constants K_300 = 300 # [unit=u"K"]
    k = a0 * exp(c0 / T) * (K_300 / T)^b0
end