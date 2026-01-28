# [NOx Photochemistry](@id nox_photochemistry)

## Overview

The NOx (NO + NO2) photochemical cycle is fundamental to tropospheric ozone
chemistry. In the absence of other species, three reactions cycle rapidly between
NO, NO2, and O3 during the day, establishing a "photostationary state" known
as the Leighton relationship (Eq. 6.6).

This system implements equations 6.5-6.8 from Section 6.2 of Seinfeld & Pandis,
including the ground-state O atom steady state, the photostationary state ozone
concentration, the photostationary state deviation parameter (Phi), and the
net ozone production rate.

Two components are provided:
- `NOxPhotochemistry`: Full diagnostic system with all four equations
- `PhotostationaryState`: Simplified system for analyzing deviations from the Leighton relationship

**Reference**: Seinfeld, J.H. and Pandis, S.N. (2006). *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change*, 2nd Edition. John Wiley & Sons. Section 6.2, pp. 207-212.

```@docs
NOxPhotochemistry
```

```@docs
PhotostationaryState
```

## Implementation

### NOxPhotochemistry: State Variables

```@example nox_phot
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities, GasChem

sys = NOxPhotochemistry()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### NOxPhotochemistry: Parameters

```@example nox_phot
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape=false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### NOxPhotochemistry: Equations

```@example nox_phot
eqs = equations(sys)
```

### PhotostationaryState: State Variables

```@example nox_phot
pss = PhotostationaryState()
vars_pss = unknowns(pss)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars_pss],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars_pss],
    :Description => [ModelingToolkit.getdescription(v) for v in vars_pss]
)
```

### PhotostationaryState: Equations

```@example nox_phot
equations(pss)
```

## Analysis

### Photostationary State O3 vs NO2/NO Ratio

The Leighton relationship (Eq. 6.6) predicts that the photostationary state
ozone concentration is proportional to the NO2/NO ratio:

``[O_3]_{pss} = \frac{j_{NO_2}}{k_{NO+O_3}} \cdot \frac{[NO_2]}{[NO]}``

This plot shows how the predicted O3 varies with the NO2/NO ratio for
typical midday photolysis conditions.

```@example nox_phot
using Plots

# Rate constants
j_NO2 = 8e-3      # NO2 photolysis rate [s^-1], midday value
k_NO_O3 = 1.8e-14 # NO + O3 rate constant at 298 K [cm3/molec/s]

# Vary NO2/NO ratio
ratio_NO2_NO = range(0.1, 10.0, length=200)

# Photostationary state O3 (Eq. 6.6)
# [O3] = j_NO2 * [NO2] / (k_NO_O3 * [NO]) = (j_NO2 / k_NO_O3) * (NO2/NO)
O3_pss = (j_NO2 / k_NO_O3) .* ratio_NO2_NO  # molec/cm3

# Convert to ppb (at STP: 1 ppb = 2.5e10 molec/cm3)
O3_ppb = O3_pss ./ 2.5e10

plot(ratio_NO2_NO, O3_ppb,
    xlabel="[NO₂]/[NO] ratio",
    ylabel="O₃ (ppb)",
    title="Photostationary State O₃ (Eq. 6.6)",
    label="O₃_pss = (j_NO₂/k_NO+O₃) × [NO₂]/[NO]",
    linewidth=2, legend=:topleft, size=(600, 400))
savefig("nox_pss_o3.svg") # hide
```

![Photostationary state O3 vs NO2/NO ratio](nox_pss_o3.svg)

The linear relationship shows that higher NO2/NO ratios lead to higher
photostationary state ozone. For a ratio of 2 (typical of moderately
polluted conditions), the predicted O3 is about 36 ppb. In real urban
environments, Phi > 1 because peroxy radicals (HO2, RO2) provide an
additional pathway for NO-to-NO2 conversion beyond the O3 + NO reaction.

### Photostationary State Parameter (Phi) Interpretation

The deviation of Phi from unity indicates the importance of peroxy radical chemistry:

```@example nox_phot
# Fixed concentrations for illustration
j_NO2_val = 8e-3
k_NO_O3_val = 1.8e-14
NO_val = 2.5e10       # 1 ppb
NO2_val = 5.0e10      # 2 ppb
O3_measured = range(0.5e12, 3e12, length=200)  # 20-120 ppb

# Phi = j_NO2 * [NO2] / (k_NO_O3 * [NO] * [O3])
Phi = j_NO2_val .* NO2_val ./ (k_NO_O3_val .* NO_val .* O3_measured)
O3_ppb_axis = O3_measured ./ 2.5e10

plot(O3_ppb_axis, Phi,
    xlabel="Measured O₃ (ppb)",
    ylabel="Φ",
    title="Photostationary State Parameter vs O₃",
    label="Φ (NO=1 ppb, NO₂=2 ppb)",
    linewidth=2, legend=:topright, size=(600, 400))
hline!([1.0], linestyle=:dash, color=:gray, label="Φ = 1 (exact PSS)")
annotate!([(90, 2.0, text("Φ > 1: net O₃ production\n(peroxy radicals active)", 8, :left)),
           (90, 0.6, text("Φ < 1: net O₃ loss", 8, :left))])
savefig("nox_phi.svg") # hide
```

![Photostationary state parameter Phi](nox_phi.svg)

When measured O3 is lower than the photostationary state value, Phi > 1,
indicating that peroxy radicals are converting NO to NO2 and producing O3.
When measured O3 exceeds the PSS value, Phi < 1. Typical urban daytime
measurements show Phi = 1.5-3, reflecting significant peroxy radical activity.
