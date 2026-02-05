# API Reference

Full API documentation for the Seinfeld & Pandis Chapter 6 tropospheric
chemistry components.

## OH Production (Section 6.1)

```@docs
OHProduction
```

## NOx Photochemistry (Section 6.2)

```@docs
NOxPhotochemistry
PhotostationaryState
```

## CO Oxidation (Section 6.3)

```@docs
COOxidation
OzoneProductionEfficiency
```

## Methane Oxidation (Section 6.4)

```@docs
MethaneOxidation
MethaneOxidationODE
```

## Combined System

```@docs
TroposphericChemistrySystem
get_typical_conditions
get_urban_conditions
get_remote_conditions
GasChem.AtmosphericLifetime
```

```@autodocs
Modules = [GasChem]
```
