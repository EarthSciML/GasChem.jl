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

### N2O as the Primary NOx Source

Nitrous oxide (N2O) transported from the troposphere is the primary source of stratospheric NOx through its reaction with excited atomic oxygen O(1D) (Section 5.3.1):

```
N2O + O(1D) -> NO + NO    (k = 6.7 × 10⁻¹¹ cm³/molec/s)
N2O + O(1D) -> N2 + O2    (k = 4.9 × 10⁻¹¹ cm³/molec/s)
```

The NO yield is approximately 58% (6.7/(6.7+4.9)), consistent with Equation 5.18.

**Reference**: Seinfeld, J.H. and Pandis, S.N. (2006). *Atmospheric Chemistry and Physics: From Air Pollution to Climate Change*, 2nd Edition. John Wiley & Sons, Chapter 5, pp. 138-203.

### Exported Functions

```@docs
ChapmanMechanism
```

```@docs
NOxCycle
```

```@docs
HOxCycle
```

```@docs
ClOxCycle
```

```@docs
BrOxCycle
```

```@docs
StratosphericOzoneSystem
```

## Implementation

The implementation uses ModelingToolkit.jl to define ODE systems for each chemical cycle. All rate coefficients are embedded as `@constants` with Arrhenius parameters within each component, ensuring that temperature-dependent rate computation is handled symbolically by ModelingToolkit. This allows for symbolic manipulation of the equations, automatic Jacobian generation, and efficient numerical integration.

### Chapman Mechanism System

```@example chapman
using GasChem
using ModelingToolkit
using ModelingToolkit: t
using DynamicQuantities
using Symbolics
using DataFrames

sys = ChapmanMechanism()
nothing # hide
```

#### State Variables

```@example chapman
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

#### Parameters

```@example chapman
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
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
using GasChem
using ModelingToolkit
using ModelingToolkit: t
using DynamicQuantities
using Symbolics
using DataFrames

sys_full = StratosphericOzoneSystem()
nothing # hide
```

#### All Species

```@example fullsys
vars = unknowns(sys_full)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

#### Parameters

```@example fullsys
params = parameters(sys_full)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

#### Equations

The full stratospheric system consists of 15 ODEs and 5 algebraic family definitions:

```@example fullsys
equations(sys_full)
```

## Analysis

This section demonstrates the behavior of the implemented systems and reproduces key results from Chapter 5 of Seinfeld & Pandis.

### Table 5.1: Stratospheric Conditions

The following table reproduces the standard atmospheric conditions used throughout Chapter 5:

```@example analysis
using DataFrames

# Table 5.1 from Seinfeld & Pandis (2006), page 142
altitudes = [20, 25, 30, 35, 40, 45]  # km
temperatures = [217, 222, 227, 237, 251, 265]  # K
pressures_hPa = [55, 25, 12, 5.6, 2.8, 1.4]  # hPa
p_over_p0 = [0.054, 0.025, 0.012, 0.0055, 0.0028, 0.0014]
M_values = [1.4e18, 6.4e17, 3.1e17, 1.4e17, 7.1e16, 3.6e16]  # molec/cm^3

# SI number densities for use in later analysis sections
M_SI = M_values .* 1e6  # molec/cm^3 → molec/m^3

# CGS to SI conversion factors
CGS_TO_SI_K2 = 1e-12  # cm^6/molec^2/s → m^6/s
CGS_TO_SI_K = 1e-6    # cm^3/molec/s → m^3/s

# Analytical rate coefficient functions from textbook Tables B.1/B.2 (CGS → SI)
k_O_O2_M_si(T) = 6.0e-34 * (T / 300.0)^(-2.4) * CGS_TO_SI_K2
k_O_O3_si(T) = 8.0e-12 * exp(-2060.0 / T) * CGS_TO_SI_K
k_NO2_O_si(T) = 5.6e-12 * exp(180.0 / T) * CGS_TO_SI_K
k_OH_O3_si(T) = 1.7e-12 * exp(-940.0 / T) * CGS_TO_SI_K
k_ClO_O_si(T) = 3.0e-11 * exp(70.0 / T) * CGS_TO_SI_K
k_Br_O3_si(T) = 1.7e-11 * exp(-800.0 / T) * CGS_TO_SI_K

table51 = DataFrame(
    Symbol("z (km)") => altitudes,
    Symbol("T (K)") => temperatures,
    Symbol("p (hPa)") => pressures_hPa,
    Symbol("p/p₀") => p_over_p0,
    Symbol("[M] (molec/cm³)") => M_values
)
```

### Table 5.2: Chemical Families in Stratospheric Chemistry

Table 5.2 from Seinfeld & Pandis (2006) summarizes the chemical families that are important in stratospheric chemistry:

```@example analysis
table52 = DataFrame(
    :Symbol => ["Oₓ", "NOₓ", "NOᵧ", "HOₓ", "Clᵧ", "ClOₓ", "Brᵧ"],
    :Name => ["Odd oxygen", "Nitrogen oxides", "Oxidized nitrogen",
        "Hydrogen radicals", "Inorganic chlorine", "Reactive chlorine",
        "Inorganic bromine"],
    :Components => [
        "O + O3",
        "NO + NO2",
        "NO + NO2 + HNO3 + 2N2O5 + ClONO2 + NO3 + HOONO2 + BrONO2",
        "OH + HO2",
        "Cl + 2Cl2 + ClO + OClO + 2Cl2O2 + HOCl + ClONO2 + HCl + BrCl",
        "Cl + ClO",
        "Br + BrO + HOBr + BrONO2"
    ]
)
```

### Rate Coefficients vs Altitude

The rate coefficients for the Chapman mechanism vary with altitude due to their temperature dependence:

```@example analysis
using Plots

k2_values = [k_O_O2_M_si(T) for T in temperatures]
k4_values = [k_O_O3_si(T) for T in temperatures]

p1 = plot(altitudes, k2_values,
    xlabel = "Altitude (km)",
    ylabel = "k₂ (m⁶ s⁻¹)",
    label = "k2: O + O2 + M",
    marker = :circle,
    title = "Temperature Dependence of Rate Coefficients",
    legend = :topright)

p2 = plot(altitudes, k4_values,
    xlabel = "Altitude (km)",
    ylabel = "k₄ (m³ s⁻¹)",
    label = "k4: O + O3",
    marker = :square,
    color = :red,
    legend = :topleft)

plot(p1, p2, layout = (1, 2), size = (800, 350))
```

### [O]/[O3] Ratio (Equation 5.7)

The steady-state ratio of atomic oxygen to ozone is given by Equation 5.7:

```math
\frac{[O]}{[O_3]} = \frac{j_{O_3}}{k_2[O_2][M]}
```

This ratio increases with altitude due to decreasing air density M:

```@example analysis
using Plots

j_O3 = 4e-4  # s^-1 (typical midday value around 30 km)
O2_frac = 0.21

O_O3_ratios = Float64[]
for (i, T) in enumerate(temperatures)
    M = M_SI[i]
    k2 = k_O_O2_M_si(T)  # m^6/s (SI)
    O2 = O2_frac * M
    ratio = j_O3 / (k2 * O2 * M)
    push!(O_O3_ratios, ratio)
end

plot(altitudes, O_O3_ratios,
    xlabel = "Altitude (km)",
    ylabel = "[O]/[O3]",
    label = "Equation 5.7",
    marker = :circle,
    yscale = :log10,
    title = "Steady-State [O]/[O3] Ratio vs Altitude",
    legend = :bottomright)
```

The [O]/[O3] ratio increases from approximately 10^-7 at 20 km to 10^-5 at 45 km, consistent with the discussion on page 144 of Seinfeld & Pandis.

### Steady-State Ozone (Equation 5.13)

The Chapman mechanism predicts a steady-state ozone concentration given by Equation 5.13:

```math
[O_3]_{ss} = 0.21 \times \sqrt{\frac{k_2 \cdot j_{O_2}}{k_4 \cdot j_{O_3}}} \times [M]^{3/2}
```

```@example analysis
using Plots

j_O2_values = [1e-12, 5e-12, 1e-11, 2e-11, 3e-11, 4e-11]  # Increases with altitude
j_O3_values = [2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4]  # s^-1

O3_ss = Float64[]
for (i, T) in enumerate(temperatures)
    M = M_SI[i]  # SI: m^-3
    k2 = k_O_O2_M_si(T)  # SI: m^6/s
    k4 = k_O_O3_si(T)    # SI: m^3/s
    j_O2 = j_O2_values[i]
    j_O3 = j_O3_values[i]

    O3 = 0.21 * sqrt(k2 * j_O2 / (k4 * j_O3)) * M^1.5  # SI: m^-3
    push!(O3_ss, O3)
end

# Convert SI m^-3 to CGS molec/cm^3 for display (÷ 1e6)
plot(O3_ss ./ 1e6 ./ 1e12, altitudes,
    ylabel = "Altitude (km)",
    xlabel = "[O3] (10¹² molec/cm³)",
    label = "Chapman prediction (Eq. 5.13)",
    marker = :circle,
    title = "Steady-State Ozone Profile",
    legend = :topright)
```

### Time to Steady State (Equation 5.17)

The characteristic time for ozone to reach steady state is given by Equation 5.17:

```math
\tau_{O_3}^{ss} = \frac{1}{4} \sqrt{\frac{k_2[M]}{k_4 \cdot j_{O_2} \cdot j_{O_3}}}
```

```@example analysis
using Plots

tau_ss = Float64[]
for (i, T) in enumerate(temperatures)
    M = M_SI[i]  # SI: m^-3
    k2 = k_O_O2_M_si(T)  # SI: m^6/s
    k4 = k_O_O3_si(T)    # SI: m^3/s
    j_O2 = j_O2_values[i]
    j_O3 = j_O3_values[i]

    tau = 0.25 * sqrt(k2 * M / (k4 * j_O2 * j_O3))
    push!(tau_ss, tau)
end

tau_days = tau_ss ./ (24 * 3600)

plot(tau_days, altitudes,
    ylabel = "Altitude (km)",
    xlabel = "Time to steady state (days)",
    label = "Equation 5.17",
    marker = :circle,
    xscale = :log10,
    title = "Ozone Steady-State Timescale",
    legend = :topright)
```

At 40 km, the ozone steady-state timescale is on the order of hours, while at 20 km it can take months to years. This explains why the ozone layer responds quickly to perturbations at higher altitudes but slowly at lower altitudes.

### Chapman Mechanism Simulation

We can solve the Chapman mechanism ODE system to observe the approach to steady state:

```@example analysis
using GasChem
using ModelingToolkit
using OrdinaryDiffEqDefault
using Plots

chapman = ChapmanMechanism()
sys = mtkcompile(chapman)

# Parameters for 30 km altitude (Table 5.1)
T_val = 227.0  # K
M_30km = 3.1e23  # m^-3 (3.1e17 molec/cm^3 × 1e6)

j_O2_val = 1e-11  # s^-1
j_O3_val = 4e-4   # s^-1

prob = ODEProblem(sys,
    [sys.O => 1e11, sys.O3 => 1e16,  # SI: m^-3
        sys.j_O2 => j_O2_val,
        sys.j_O3 => j_O3_val,
        sys.T => T_val,
        sys.M => M_30km,
        sys.O2_mix => 0.21],
    (0.0, 3600.0 * 24 * 10))  # 10 days

sol = solve(prob, abstol = 1e-8, reltol = 1e-8)

time_hours = sol.t ./ 3600

# Convert SI m^-3 back to CGS molec/cm^3 for display (÷ 1e6)
p1 = plot(time_hours, sol[sys.O3] ./ 1e6 ./ 1e12,
    xlabel = "Time (hours)",
    ylabel = "[O3] (10¹² molec/cm³)",
    label = "O3",
    title = "Chapman Mechanism: Approach to Steady State (30 km)")

p2 = plot(time_hours, sol[sys.O] ./ 1e6,
    xlabel = "Time (hours)",
    ylabel = "[O] (molec/cm³)",
    label = "O",
    color = :red,
    title = "Atomic Oxygen")

plot(p1, p2, layout = (2, 1), size = (700, 500))
```

### Catalytic Cycle Contributions to Ozone Loss

The catalytic cycles (NOx, HOx, ClOx, BrOx) each contribute to ozone destruction at different rates depending on altitude. The relative importance of each cycle varies with altitude, as shown in Figure 5.29 of Seinfeld & Pandis.

At lower stratospheric altitudes (15-25 km), the HOx cycle dominates ozone loss. At middle altitudes (25-40 km), the NOx cycle becomes most important. The ClOx cycle contributes significantly at all altitudes and is particularly important in the polar regions where heterogeneous chemistry activates chlorine reservoirs.

The following figure computes the relative destruction rates at each altitude, assuming typical species mixing ratios from the textbook:

```@example analysis
using Plots

# Compute catalytic destruction rates
# Using typical stratospheric mixing ratios from Chapter 5
# All concentrations in SI (m^-3)

# Typical mixing ratios (from text)
NO2_vmr = 5e-9   # 5 ppbv
OH_vmr = 5e-7    # 0.5 pptv — very low
ClO_vmr = 1e-10  # ~100 pptv
BrO_vmr = 5e-12  # ~5 pptv

# Reference values in SI for scaling (30 km)
M_ref_si = M_SI[3]

nox_rates = Float64[]
hox_rates = Float64[]
clox_rates = Float64[]
brox_rates = Float64[]
chapman_rates = Float64[]

for (i, T) in enumerate(temperatures)
    M = M_SI[i]  # SI: m^-3
    O_conc = 1e13 * (M_ref_si / M)     # Scale O with inverse density (SI: m^-3)
    O3_conc = 3e18 * (M / M_ref_si)    # Scale O3 with density (SI: m^-3)

    NO2 = NO2_vmr * M
    OH = OH_vmr * M
    ClO = ClO_vmr * M
    BrO = BrO_vmr * M

    # Destruction rates (from Seinfeld & Pandis Eqs. 5.22, 5.25, 5.29)
    k_NO2_O_val = k_NO2_O_si(T)   # SI: m^3/s
    k_OH_O3_val = k_OH_O3_si(T)
    k_ClO_O_val = k_ClO_O_si(T)
    k_Br_O3_val = k_Br_O3_si(T)
    k4 = k_O_O3_si(T)

    R_NOx = 2 * k_NO2_O_val * NO2 * O_conc
    R_HOx = 2 * k_OH_O3_val * OH * O3_conc
    R_ClOx = 2 * k_ClO_O_val * ClO * O_conc
    R_BrOx = 2 * k_Br_O3_val * BrO * O3_conc
    R_chapman = 2 * k4 * O_conc * O3_conc

    total = R_NOx + R_HOx + R_ClOx + R_BrOx + R_chapman
    push!(nox_rates, 100 * R_NOx / total)
    push!(hox_rates, 100 * R_HOx / total)
    push!(clox_rates, 100 * R_ClOx / total)
    push!(brox_rates, 100 * R_BrOx / total)
    push!(chapman_rates, 100 * R_chapman / total)
end

plot(nox_rates, altitudes,
    label = "NOx",
    xlabel = "Contribution to O3 Loss (%)",
    ylabel = "Altitude (km)",
    title = "Catalytic Cycle Contributions (cf. Figure 5.29)",
    linewidth = 2,
    legend = :topright)
plot!(hox_rates, altitudes, label = "HOx", linewidth = 2)
plot!(clox_rates, altitudes, label = "ClOx", linewidth = 2)
plot!(brox_rates, altitudes, label = "BrOx", linewidth = 2)
plot!(chapman_rates, altitudes, label = "Chapman (O+O3)", linewidth = 2, linestyle = :dash)
```

### NOx Catalytic Efficiency

The NOx cycle destroys odd oxygen at the rate (Equation 5.22):

```math
\frac{d[O_x]}{dt} = -2k_2[NO_2][O]
```

The efficiency of the NOx cycle can be estimated by comparing this destruction rate to the Chapman mechanism rate:

```@example analysis
using Plots

NO2_mix = 1e-9  # 1 ppbv
O_mix = 1e-11   # varies strongly with altitude

ratios = Float64[]
for (i, T) in enumerate(temperatures)
    M = M_SI[i]  # SI: m^-3
    k4 = k_O_O3_si(T)      # SI: m^3/s
    k_NO2_O_val = k_NO2_O_si(T)  # SI: m^3/s

    O = O_mix * M
    O3 = 3e18  # SI: m^-3 (3e12 molec/cm^3 × 1e6)
    NO2 = NO2_mix * M

    R_chapman = k4 * O * O3
    R_NOx = k_NO2_O_val * NO2 * O

    push!(ratios, R_NOx / R_chapman)
end

plot(altitudes, ratios,
    xlabel = "Altitude (km)",
    ylabel = "NOx/Chapman Rate Ratio",
    label = "Rate ratio",
    marker = :circle,
    title = "NOx Cycle Efficiency Relative to Chapman",
    legend = :topleft,
    yscale = :log10)
```

### Chlorine Reservoir Partitioning

In the stratosphere, most chlorine exists in reservoir species (HCl and ClONO2) rather than reactive forms (Cl and ClO). The partitioning between these species is crucial for understanding ozone depletion.

The following figure uses the ClOx cycle system from the implementation to simulate the approach to steady-state partitioning:

```@example analysis
using GasChem
using ModelingToolkit
using OrdinaryDiffEqDefault
using Plots

# Run ClOx system at different altitudes and examine the partitioning
sys = ClOxCycle()
compiled = mtkcompile(sys)

# Compute partitioning at a representative altitude (30 km)
# CGS→SI: 1 cm^-3 = 1e6 m^-3
Cly_total = 2e15  # Total Cly in SI m^-3 (2e9 molec/cm^3 × 1e6)

prob = ODEProblem(compiled,
    [compiled.Cl => 1e10, compiled.ClO => 1e13,
        compiled.HCl => 1e15, compiled.ClONO2 => Cly_total - 1e10 - 1e13 - 1e15],
    (0.0, 3600.0 * 24)  # 1 day
)

sol = solve(prob, abstol = 1e-8, reltol = 1e-8)

time_hours = sol.t ./ 3600

Cly = sol[compiled.Cl] .+ sol[compiled.ClO] .+ sol[compiled.HCl] .+ sol[compiled.ClONO2]
HCl_frac = sol[compiled.HCl] ./ Cly
ClONO2_frac = sol[compiled.ClONO2] ./ Cly
ClOx_frac = (sol[compiled.Cl] .+ sol[compiled.ClO]) ./ Cly

plot(time_hours, HCl_frac, label = "HCl/Cly", linewidth = 2,
    xlabel = "Time (hours)", ylabel = "Fraction of Cly",
    title = "Chlorine Reservoir Partitioning (30 km)",
    legend = :right)
plot!(time_hours, ClONO2_frac, label = "ClONO2/Cly", linewidth = 2)
plot!(time_hours, ClOx_frac, label = "ClOx/Cly", linewidth = 2)
```

The reactive chlorine species (ClOx = Cl + ClO) constitute only a small fraction of total inorganic chlorine (Cly) under normal conditions. However, heterogeneous reactions on polar stratospheric clouds can rapidly convert HCl and ClONO2 to reactive forms, leading to severe ozone depletion in polar regions (the "ozone hole").

### Summary

This module provides a comprehensive implementation of stratospheric ozone chemistry based on the well-established theory presented in Seinfeld & Pandis (2006). The key features include:

 1. **Chapman Mechanism**: Fundamental O2/O3 photochemistry with temperature-dependent rate coefficients
 2. **Catalytic Cycles**: NOx, HOx, ClOx, and BrOx destruction mechanisms
 3. **Chemical Families**: Tracking of Ox, NOx, HOx, ClOx, and BrOx families
 4. **N2O Source Chemistry**: N2O + O(1D) as the primary source of stratospheric NOx (Section 5.3.1)
 5. **O(1D) Photochemistry**: Including reactions with H2O, CH4, and N2O
 6. **Rate Coefficients**: All Arrhenius parameters embedded as `@constants` within ModelingToolkit components
