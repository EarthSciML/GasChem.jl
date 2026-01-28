# [CO Oxidation](@id co_oxidation)

## Overview

CO oxidation is the simplest hydrocarbon oxidation mechanism and illustrates
the fundamental HOx cycling that drives tropospheric ozone production. The
net reaction is:

    CO + 2O2 + hv -> CO2 + O3

The key feature is that OH is regenerated through the HO2 + NO reaction,
allowing catalytic ozone production. The HOx chain length (number of
OH-HO2 cycles before termination) and the Ozone Production Efficiency
(OPE = moles O3 produced per mole NOx consumed) are key diagnostics.

Two components are provided:
- `COOxidation`: Full CO oxidation diagnostic system (Eqs. 6.9-6.17)
- `OzoneProductionEfficiency`: OPE diagnostic (Eqs. 6.21-6.24)

**Reference**: Seinfeld, J.H. and Pandis, S.N. (2006). *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change*, 2nd Edition. John Wiley & Sons. Section 6.3, pp. 212-219.

```@docs
COOxidation
```

```@docs
OzoneProductionEfficiency
```

## Implementation

### COOxidation: State Variables

```@example co_ox
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities, GasChem

sys = COOxidation()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### COOxidation: Parameters

```@example co_ox
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape=false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### COOxidation: Equations

```@example co_ox
eqs = equations(sys)
```

### OzoneProductionEfficiency: State Variables

```@example co_ox
ope_sys = OzoneProductionEfficiency()
vars_ope = unknowns(ope_sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars_ope],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars_ope],
    :Description => [ModelingToolkit.getdescription(v) for v in vars_ope]
)
```

### OzoneProductionEfficiency: Equations

```@example co_ox
equations(ope_sys)
```

## Analysis

### HOx Chain Length and OPE vs NOx

The HOx chain length and Ozone Production Efficiency (OPE) both depend
strongly on the NOx concentration. At high NOx, the dominant HOx termination
is OH + NO2 -> HNO3, giving short chains and low OPE. At low NOx, the
dominant termination is HO2 + HO2 -> H2O2, giving long chains and high OPE.

```@example co_ox
using Plots

# Rate constants at 298 K
k_CO_OH = 2.4e-13    # CO + OH [cm3/molec/s]
k_HO2_NO = 8.1e-12   # HO2 + NO [cm3/molec/s]
k_HO2_HO2 = 2.9e-12  # HO2 + HO2 [cm3/molec/s]
k_OH_NO2 = 1.0e-11   # OH + NO2 [cm3/molec/s]
k_HO2_O3 = 2.0e-15   # HO2 + O3 [cm3/molec/s]
k_OH_O3 = 7.3e-14    # OH + O3 [cm3/molec/s]

# Fixed conditions
CO = 2.5e12    # 100 ppb
O3 = 1e12      # 40 ppb
OH = 1e6       # typical
P_OH = 1e6     # molec/cm3/s

# Vary NOx (expressed as NO, ppb)
NO_ppb = 10 .^ range(-2, 2, length=200)
NO = NO_ppb .* 2.5e10  # convert to molec/cm3
NO2 = 2 .* NO          # assume NO2/NO = 2

# Estimate HO2 from steady state:
# At high NOx: HO2 ~ k_CO_OH * CO * OH / (k_HO2_NO * NO)
# At low NOx: HO2 ~ sqrt(P_OH / (2 * k_HO2_HO2))
HO2_high = k_CO_OH .* CO .* OH ./ (k_HO2_NO .* NO)
HO2_low = sqrt(P_OH / (2 * k_HO2_HO2))
HO2 = min.(HO2_high, HO2_low)

# HOx loss rate
L_HOx = k_OH_NO2 .* OH .* NO2 .+ 2 .* k_HO2_HO2 .* HO2.^2

# HOx chain length
chain = (k_HO2_NO .* HO2 .* NO) ./ L_HOx

# OPE = P(O3) / L(NOx)
P_O3 = k_HO2_NO .* HO2 .* NO
L_NOx = k_OH_NO2 .* OH .* NO2
OPE = P_O3 ./ L_NOx

p1 = plot(NO_ppb, chain,
    xlabel="NO (ppb)", ylabel="Chain Length",
    title="HOx Chain Length vs NO",
    xscale=:log10, yscale=:log10,
    linewidth=2, label="Chain length",
    legend=:topright, ylims=(0.1, 1000))

p2 = plot(NO_ppb, OPE,
    xlabel="NO (ppb)", ylabel="OPE (mol O₃ / mol NOx)",
    title="Ozone Production Efficiency vs NO",
    xscale=:log10, yscale=:log10,
    linewidth=2, label="OPE",
    legend=:topright, ylims=(0.1, 1000))

plot(p1, p2, layout=(1, 2), size=(800, 350))
savefig("co_chain_ope.svg") # hide
```

![HOx chain length and OPE vs NOx](co_chain_ope.svg)

Both the chain length and OPE decrease as NOx increases, reflecting the
transition from the NOx-limited regime (low NOx, HO2 + HO2 termination,
OPE ~ 10-30) to the VOC-limited regime (high NOx, OH + NO2 termination,
OPE ~ 1-3). This is consistent with Table 6.2 of Seinfeld & Pandis.

### Net O3 Production Rate vs NO

The net O3 production rate initially increases with NO (more HO2 + NO
reactions) but levels off or decreases at very high NO as O3 titration
(NO + O3) and reduced HO2 (due to faster cycling) compete.

```@example co_ox
# Net O3 production = HO2*NO - OH*O3 - HO2*O3 loss terms
P_O3_net = k_HO2_NO .* HO2 .* NO .- k_OH_O3 .* OH .* O3 .- k_HO2_O3 .* HO2 .* O3

plot(NO_ppb, P_O3_net ./ 1e6,
    xlabel="NO (ppb)",
    ylabel="Net P(O₃) (10⁶ molec cm⁻³ s⁻¹)",
    title="Net O₃ Production Rate vs NO",
    xscale=:log10,
    linewidth=2, label="P(O₃)_net",
    legend=:topleft, size=(600, 400))
savefig("co_po3_vs_no.svg") # hide
```

![Net O3 production rate vs NO](co_po3_vs_no.svg)

This figure illustrates the nonlinear dependence of O3 production on NOx,
a fundamental feature of tropospheric photochemistry described in Section 6.3.
