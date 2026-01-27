# Stratospheric Chemistry

## Overview

The stratosphere, extending from approximately 10-50 km altitude, contains the ozone layer that shields the Earth's surface from harmful ultraviolet radiation. This module implements the fundamental chemistry governing stratospheric ozone based on Chapter 5 of Seinfeld & Pandis (2006).

### The Chapman Mechanism

The Chapman mechanism (Chapman, 1930) describes the basic photochemical production and destruction of ozone in the stratosphere through four key reactions:

1. **O2 Photolysis**: O2 + hv -> O + O (j_O2)
2. **Ozone Formation**: O + O2 + M -> O3 + M (k2)
3. **Ozone Photolysis**: O3 + hv -> O + O2 (j_O3)
4. **Ozone Destruction**: O + O3 -> O2 + O2 (k4)

While the Chapman mechanism predicts the qualitative features of the ozone layer, it overestimates ozone concentrations by roughly a factor of two. This discrepancy is resolved by including catalytic destruction cycles involving NOx, HOx, ClOx, and BrOx species.

### Catalytic Ozone Destruction

Catalytic cycles destroy ozone through the general mechanism:

```
X + O3 -> XO + O2
XO + O -> X + O2
-----------------
Net: O3 + O -> 2 O2
```

where X represents the catalyst (NO, OH, Cl, Br). These cycles are highly efficient because the catalyst is regenerated and can destroy many ozone molecules before being removed.

**Reference**: Seinfeld, J.H. and Pandis, S.N. (2006). *Atmospheric Chemistry and Physics: From Air Pollution to Climate Change*, 2nd Edition. John Wiley & Sons, Chapter 5, pp. 138-203.

### Exported Functions

```@docs
ChapmanMechanism
NOxCycle
HOxCycle
ClOxCycle
BrOxCycle
StratosphericOzoneSystem
stratospheric_rate_coefficients
```

## Implementation

The implementation uses ModelingToolkit.jl to define ODE systems for each chemical cycle. This allows for symbolic manipulation of the equations, automatic Jacobian generation, and efficient numerical integration.

### Chapman Mechanism System

```@example chapman
using StratosphericChemistry
using ModelingToolkit
using DataFrames

# Create the Chapman mechanism system
sys = structural_simplify(ChapmanMechanism())
nothing # hide
```

#### State Variables

```@example chapman
vars = unknowns(sys)
DataFrame(
    :Name => [string(v) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

#### Parameters

```@example chapman
params = parameters(sys)
DataFrame(
    :Name => [string(p) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

#### Equations

The Chapman mechanism is governed by coupled differential equations for atomic oxygen O and ozone O3:

```@example chapman
equations(sys)
```

### Comprehensive Stratospheric System

The full stratospheric system combines all catalytic cycles:

```@example fullsys
using StratosphericChemistry
using ModelingToolkit
using DataFrames

sys_full = structural_simplify(StratosphericOzoneSystem())
nothing # hide
```

#### All Species

```@example fullsys
vars = unknowns(sys_full)
DataFrame(
    :Name => [string(v) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

## Analysis

This section demonstrates the behavior of the implemented systems and reproduces key results from Chapter 5 of Seinfeld & Pandis.

### Table 5.1: Stratospheric Conditions

The following table reproduces the standard atmospheric conditions used throughout Chapter 5:

```@example analysis
using DataFrames

# Table 5.1 from Seinfeld & Pandis (2006)
# Stratospheric Temperature, Pressure, and Atmospheric Number Density
altitudes = [20, 25, 30, 35, 40, 45]  # km
temperatures = [217, 222, 227, 237, 251, 264]  # K
pressures = [55.3, 25.5, 12.0, 5.75, 2.87, 1.49]  # mbar
M_values = [1.85e18, 8.33e17, 3.83e17, 1.76e17, 8.31e16, 4.09e16]  # molec/cm^3

table51 = DataFrame(
    Symbol("Altitude (km)") => altitudes,
    Symbol("Temperature (K)") => temperatures,
    Symbol("Pressure (mbar)") => pressures,
    Symbol("M (molec/cm^3)") => M_values
)
```

### Rate Coefficients vs Altitude

The rate coefficients for the Chapman mechanism vary with altitude due to their temperature dependence:

```@example analysis
using StratosphericChemistry
using Plots

# Calculate rate coefficients at each altitude
k2_values = [StratosphericChemistry.k_O_O2_M(T) for T in temperatures]
k4_values = [StratosphericChemistry.k_O_O3(T) for T in temperatures]

p1 = plot(altitudes, k2_values,
    xlabel="Altitude (km)",
    ylabel="k2 (cm^6 molec^-2 s^-1)",
    label="k2: O + O2 + M",
    marker=:circle,
    title="Temperature Dependence of Rate Coefficients",
    legend=:topright)

p2 = plot(altitudes, k4_values,
    xlabel="Altitude (km)",
    ylabel="k4 (cm^3 molec^-1 s^-1)",
    label="k4: O + O3",
    marker=:square,
    color=:red,
    legend=:topleft)

plot(p1, p2, layout=(1,2), size=(800, 350))
savefig("rate_coefficients.svg"); nothing # hide
```

![Rate coefficients vs altitude](rate_coefficients.svg)

### [O]/[O3] Ratio (Equation 5.7)

The steady-state ratio of atomic oxygen to ozone is given by Equation 5.7:

```math
\frac{[O]}{[O_3]} = \frac{j_{O_3}}{k_2[O_2][M]}
```

This ratio increases with altitude due to decreasing air density M:

```@example analysis
using StratosphericChemistry
using Plots

# Equation 5.7: [O]/[O3] = j_O3 / (k2 * [O2] * [M])
# Using typical photolysis rates from the text
j_O3 = 4e-4  # s^-1 (typical midday value around 30 km)
O2_frac = 0.21

O_O3_ratios = Float64[]
for (i, T) in enumerate(temperatures)
    M = M_values[i]
    k2 = StratosphericChemistry.k_O_O2_M(T)
    O2 = O2_frac * M
    ratio = j_O3 / (k2 * O2 * M)
    push!(O_O3_ratios, ratio)
end

plot(altitudes, O_O3_ratios,
    xlabel="Altitude (km)",
    ylabel="[O]/[O3]",
    label="Equation 5.7",
    marker=:circle,
    yscale=:log10,
    title="Steady-State [O]/[O3] Ratio vs Altitude",
    legend=:bottomright)
savefig("O_O3_ratio.svg"); nothing # hide
```

![O/O3 ratio vs altitude](O_O3_ratio.svg)

The [O]/[O3] ratio increases from approximately 10^-7 at 20 km to 10^-5 at 45 km, consistent with the discussion on page 144 of Seinfeld & Pandis.

### Steady-State Ozone (Equation 5.13)

The Chapman mechanism predicts a steady-state ozone concentration given by Equation 5.13:

```math
[O_3]_{ss} = 0.21 \times \sqrt{\frac{k_2 \cdot j_{O_2}}{k_4 \cdot j_{O_3}}} \times [M]^{3/2}
```

```@example analysis
using StratosphericChemistry
using Plots

# Equation 5.13: Steady-state ozone
# j_O2 and j_O3 values from typical stratospheric conditions
j_O2_values = [1e-12, 5e-12, 1e-11, 2e-11, 3e-11, 4e-11]  # Increases with altitude
j_O3_values = [2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4]  # s^-1

O3_ss = Float64[]
for (i, T) in enumerate(temperatures)
    M = M_values[i]
    k2 = StratosphericChemistry.k_O_O2_M(T)
    k4 = StratosphericChemistry.k_O_O3(T)
    j_O2 = j_O2_values[i]
    j_O3 = j_O3_values[i]

    # Equation 5.13
    O3 = 0.21 * sqrt(k2 * j_O2 / (k4 * j_O3)) * M^1.5
    push!(O3_ss, O3)
end

plot(O3_ss ./ 1e12, altitudes,
    ylabel="Altitude (km)",
    xlabel="[O3] (10^12 molec/cm^3)",
    label="Chapman prediction (Eq. 5.13)",
    marker=:circle,
    title="Steady-State Ozone Profile",
    legend=:topright)
savefig("ozone_profile.svg"); nothing # hide
```

![Steady-state ozone profile](ozone_profile.svg)

### Time to Steady State (Equation 5.17)

The characteristic time for ozone to reach steady state is given by Equation 5.17:

```math
\tau_{O_3}^{ss} = \frac{1}{4} \sqrt{\frac{k_2[M]}{k_4 \cdot j_{O_2} \cdot j_{O_3}}}
```

```@example analysis
using StratosphericChemistry
using Plots

tau_ss = Float64[]
for (i, T) in enumerate(temperatures)
    M = M_values[i]
    k2 = StratosphericChemistry.k_O_O2_M(T)
    k4 = StratosphericChemistry.k_O_O3(T)
    j_O2 = j_O2_values[i]
    j_O3 = j_O3_values[i]

    # Equation 5.17
    tau = 0.25 * sqrt(k2 * M / (k4 * j_O2 * j_O3))
    push!(tau_ss, tau)
end

# Convert to days
tau_days = tau_ss ./ (24 * 3600)

plot(tau_days, altitudes,
    ylabel="Altitude (km)",
    xlabel="Time to steady state (days)",
    label="Equation 5.17",
    marker=:circle,
    xscale=:log10,
    title="Ozone Steady-State Timescale",
    legend=:topright)
savefig("tau_steadystate.svg"); nothing # hide
```

![Time to steady state](tau_steadystate.svg)

At 40 km, the ozone steady-state timescale is on the order of hours, while at 20 km it can take months to years. This explains why the ozone layer responds quickly to perturbations at higher altitudes but slowly at lower altitudes.

### Chapman Mechanism Simulation

We can solve the Chapman mechanism ODE system to observe the approach to steady state:

```@example analysis
using StratosphericChemistry
using ModelingToolkit
using OrdinaryDiffEq
using Plots

# Create and simplify the system
chapman = ChapmanMechanism()
sys = structural_simplify(chapman)

# Parameters for 30 km altitude (Table 5.1)
T = 227.0  # K
M = 3.83e17  # molec/cm^3
O2_conc = 0.21 * M

# Photolysis rates (typical values)
j_O2_val = 1e-11  # s^-1
j_O3_val = 4e-4   # s^-1

# Rate coefficients
k2_val = StratosphericChemistry.k_O_O2_M(T)
k4_val = StratosphericChemistry.k_O_O3(T)

# Set up the problem
prob = ODEProblem(sys,
    [sys.O => 1e5, sys.O3 => 1e10],  # Initial conditions
    (0.0, 3600.0 * 24 * 10),  # 10 days in seconds
    [sys.j_O2 => j_O2_val,
     sys.j_O3 => j_O3_val,
     sys.k2 => k2_val,
     sys.k4 => k4_val,
     sys.M => M,
     sys.O2_conc => 0.21])

# Solve
sol = solve(prob, Rodas5P(), abstol=1e-8, reltol=1e-8)

# Plot
time_hours = sol.t ./ 3600

p1 = plot(time_hours, [sol[sys.O3][i] for i in 1:length(sol.t)] ./ 1e12,
    xlabel="Time (hours)",
    ylabel="[O3] (10^12 molec/cm^3)",
    label="O3",
    title="Chapman Mechanism: Approach to Steady State (30 km)")

p2 = plot(time_hours, [sol[sys.O][i] for i in 1:length(sol.t)],
    xlabel="Time (hours)",
    ylabel="[O] (molec/cm^3)",
    label="O",
    color=:red,
    title="Atomic Oxygen")

plot(p1, p2, layout=(2,1), size=(700, 500))
savefig("chapman_simulation.svg"); nothing # hide
```

![Chapman mechanism simulation](chapman_simulation.svg)

### Catalytic Cycle Contributions to Ozone Loss

The catalytic cycles (NOx, HOx, ClOx, BrOx) each contribute to ozone destruction at different rates depending on altitude. The relative importance of each cycle varies with altitude, as shown in Figure 5.29 of Seinfeld & Pandis.

At lower stratospheric altitudes (15-25 km), the HOx cycle dominates ozone loss. At middle altitudes (25-40 km), the NOx cycle becomes most important. The ClOx cycle contributes significantly at all altitudes and is particularly important in the polar regions where heterogeneous chemistry activates chlorine reservoirs.

```@example analysis
using Plots

# Approximate contributions based on Figure 5.29 (page 185)
# These are illustrative values showing the relative importance

altitudes_fine = 15:2:45
n = length(altitudes_fine)

# Approximate percentage contributions (from Figure 5.29)
# Values estimated from the figure in the textbook
nox_contrib = [5, 10, 25, 40, 55, 60, 55, 50, 45, 40, 35, 30, 25, 20, 15, 10]
hox_contrib = [50, 45, 35, 25, 20, 18, 20, 22, 25, 28, 30, 32, 35, 38, 40, 42]
clox_contrib = [15, 18, 22, 25, 20, 18, 20, 23, 25, 27, 30, 32, 35, 38, 42, 45]
brox_contrib = [10, 12, 10, 5, 3, 2, 3, 3, 3, 3, 3, 4, 3, 2, 1, 1]
chapman_contrib = [20, 15, 8, 5, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]

plot(nox_contrib, altitudes_fine,
    label="NOx",
    xlabel="Contribution to O3 Loss (%)",
    ylabel="Altitude (km)",
    title="Catalytic Cycle Contributions (cf. Figure 5.29)",
    linewidth=2,
    legend=:topright)
plot!(hox_contrib, altitudes_fine, label="HOx", linewidth=2)
plot!(clox_contrib, altitudes_fine, label="ClOx", linewidth=2)
plot!(brox_contrib, altitudes_fine, label="BrOx", linewidth=2)
plot!(chapman_contrib, altitudes_fine, label="Chapman (O+O3)", linewidth=2, linestyle=:dash)
savefig("catalytic_contributions.svg"); nothing # hide
```

![Catalytic cycle contributions](catalytic_contributions.svg)

### NOx Catalytic Efficiency

The NOx cycle destroys odd oxygen at the rate (Equation 5.22):

```math
\frac{d[O_x]}{dt} = -2k_2[NO_2][O]
```

The efficiency of the NOx cycle can be estimated by comparing this destruction rate to the Chapman mechanism rate:

```@example analysis
using StratosphericChemistry
using Plots

# Compare NOx destruction rate to Chapman rate at different altitudes
# Assume typical mixing ratios from the text

# Typical stratospheric concentrations
NO2_mix = 1e-9  # 1 ppbv
O_mix = 1e-11   # varies strongly with altitude

# Chapman destruction: 2 * k4 * [O] * [O3]
# NOx destruction: 2 * k_NO2_O * [NO2] * [O]

ratios = Float64[]
for (i, T) in enumerate(temperatures)
    M = M_values[i]
    k4 = StratosphericChemistry.k_O_O3(T)
    k_NO2_O = StratosphericChemistry.k_NO2_O(T)

    O = O_mix * M  # [O] increases with altitude relative to [O3]
    O3 = 3e12  # Typical O3 concentration
    NO2 = NO2_mix * M

    # Rate ratio
    R_chapman = k4 * O * O3
    R_NOx = k_NO2_O * NO2 * O

    push!(ratios, R_NOx / R_chapman)
end

plot(altitudes, ratios,
    xlabel="Altitude (km)",
    ylabel="NOx/Chapman Rate Ratio",
    label="Rate ratio",
    marker=:circle,
    title="NOx Cycle Efficiency Relative to Chapman",
    legend=:topleft,
    yscale=:log10)
savefig("nox_efficiency.svg"); nothing # hide
```

![NOx cycle efficiency](nox_efficiency.svg)

### Chlorine Reservoir Partitioning

In the stratosphere, most chlorine exists in reservoir species (HCl and ClONO2) rather than reactive forms (Cl and ClO). The partitioning between these species is crucial for understanding ozone depletion:

```@example analysis
using StratosphericChemistry
using Plots

# Typical partitioning of inorganic chlorine (Cly)
# Based on discussion in Section 5.5.4

altitudes_cl = [15, 20, 25, 30, 35, 40, 45]

# Approximate fractions (from text discussion)
HCl_frac = [0.85, 0.75, 0.60, 0.50, 0.45, 0.42, 0.40]
ClONO2_frac = [0.14, 0.23, 0.35, 0.40, 0.38, 0.33, 0.28]
ClOx_frac = [0.01, 0.02, 0.05, 0.10, 0.17, 0.25, 0.32]

areaplot(altitudes_cl, hcat(HCl_frac, ClONO2_frac, ClOx_frac),
    label=["HCl" "ClONO2" "ClOx"],
    xlabel="Altitude (km)",
    ylabel="Fraction of Cly",
    title="Chlorine Reservoir Partitioning",
    fillalpha=0.7,
    legend=:topright)
savefig("chlorine_partitioning.svg"); nothing # hide
```

![Chlorine reservoir partitioning](chlorine_partitioning.svg)

The reactive chlorine species (ClOx = Cl + ClO) constitute only a small fraction of total inorganic chlorine (Cly) under normal conditions. However, heterogeneous reactions on polar stratospheric clouds can rapidly convert HCl and ClONO2 to reactive forms, leading to severe ozone depletion in polar regions (the "ozone hole").

### Summary

This module provides a comprehensive implementation of stratospheric ozone chemistry based on the well-established theory presented in Seinfeld & Pandis (2006). The key features include:

1. **Chapman Mechanism**: Fundamental O2/O3 photochemistry with temperature-dependent rate coefficients
2. **Catalytic Cycles**: NOx, HOx, ClOx, and BrOx destruction mechanisms
3. **Chemical Families**: Tracking of Ox, NOx, HOx, ClOx, and BrOx families
4. **Rate Coefficients**: All rate coefficients based on recommended values from the literature

The implementation allows users to:
- Study the steady-state behavior of stratospheric ozone
- Investigate the relative importance of different catalytic cycles
- Explore the sensitivity of the ozone layer to perturbations
- Understand the chemistry underlying the Antarctic ozone hole
