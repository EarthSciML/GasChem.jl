# [OH Production](@id oh_production)

## Overview

The hydroxyl radical (OH) is the primary oxidant in the troposphere, initiating
the oxidation of most trace gases. The dominant source of OH is photolysis of
ozone at wavelengths below 320 nm, producing electronically excited O(1D) atoms
which can react with water vapor to form OH.

This system implements equations 6.1-6.4 from Section 6.1 of Seinfeld & Pandis,
computing the steady-state O(1D) concentration, the OH yield (fraction of O(1D)
that produces OH rather than being quenched), and the OH production rate.

**Reference**: Seinfeld, J.H. and Pandis, S.N. (2006). *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change*, 2nd Edition. John Wiley & Sons. Section 6.1, pp. 204-207.

```@docs
OHProduction
```

## Implementation

The `OHProduction` component is an algebraic system that computes diagnostic
quantities from input concentrations of O3, H2O, and the total air number
density M.

### State Variables

```@example oh_prod
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities, GasChem

sys = OHProduction()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example oh_prod
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape=false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example oh_prod
eqs = equations(sys)
```

## Analysis

### OH Yield vs Water Vapor Concentration

The OH yield (epsilon_OH) represents the fraction of O(1D) atoms that react with H2O
to produce OH rather than being quenched back to O(3P). It increases with humidity
because more water vapor competes with N2/O2 for the O(1D) atoms.

From Equation 6.4:
``\varepsilon_{OH} = \frac{k_4 [H_2O]}{k_3 [M] + k_4 [H_2O]}``

Typical tropospheric values range from about 0.05 (dry, upper troposphere) to
0.15 (humid, lower troposphere).

```@example oh_prod
using Plots

# Rate constants at 298 K
k3_N2 = 2.6e-11  # O(1D) + N2 quenching rate [cm3/molec/s]
k3_O2 = 4.0e-11  # O(1D) + O2 quenching rate [cm3/molec/s]
k4 = 2.2e-10     # O(1D) + H2O rate [cm3/molec/s]
f_N2 = 0.78
f_O2 = 0.21
k3_eff = f_N2 * k3_N2 + f_O2 * k3_O2
M = 2.5e19  # total air [molec/cm3]

# Vary H2O from dry to very humid conditions
H2O = range(1e16, 1.5e18, length=200)  # [molec/cm3]

# Compute OH yield (Eq. 6.4)
eps_OH = [k4 * h / (k3_eff * M + k4 * h) for h in H2O]

# Compute OH production rate (Eq. 6.3) for O3 = 40 ppb
j_O3 = 1e-5  # s^-1
O3 = 1e12    # ~40 ppb [molec/cm3]
P_OH = [2 * j_O3 * O3 * e for e in eps_OH]

p1 = plot(H2O ./ 1e17, eps_OH .* 100,
    xlabel="[H₂O] (10¹⁷ molec cm⁻³)",
    ylabel="OH Yield (%)",
    title="OH Yield vs Water Vapor",
    label="ε_OH (Eq. 6.4)",
    linewidth=2, legend=:bottomright)

p2 = plot(H2O ./ 1e17, P_OH ./ 1e6,
    xlabel="[H₂O] (10¹⁷ molec cm⁻³)",
    ylabel="P(OH) (10⁶ molec cm⁻³ s⁻¹)",
    title="OH Production Rate vs Humidity",
    label="P_OH (Eq. 6.3, O₃=40 ppb)",
    linewidth=2, legend=:bottomright)

plot(p1, p2, layout=(1, 2), size=(800, 350))
savefig("oh_yield_humidity.svg") # hide
```

![OH yield and production rate vs humidity](oh_yield_humidity.svg)

The left panel shows that OH yield increases approximately linearly with water
vapor at typical tropospheric concentrations (where ``k_3[M] \gg k_4[H_2O]``),
ranging from about 3% in dry air to about 12% at high humidity. The right panel
shows the corresponding OH production rate, which is proportional to both the
O3 concentration and the OH yield.
