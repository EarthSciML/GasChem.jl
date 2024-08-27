# Chemical Species

Here is a list of the chemical species in the mechanism:

```@example 1
using GasChem, DataFrames, EarthSciMLBase, ModelingToolkit, DynamicQuantities
using ModelingToolkit:t
gc = structural_simplify(GEOSChemGasPhase())
vars = unknowns(gc)
DataFrame(
        :Name => [string(Symbolics.tosymbol(v, escape=false)) for v ∈ vars],
        :Units => [ModelingToolkit.get_unit(v) for v ∈ vars],
        :Description => [ModelingToolkit.getdescription(v) for v ∈ vars],
        :Default => [ModelingToolkit.getdefault(v) for v ∈ vars])
```