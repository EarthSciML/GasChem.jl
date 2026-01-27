# Atmospheric Lifetime

## Overview

The AtmosphericLifetime module provides a ModelingToolkit.jl implementation of atmospheric constituent lifetime equations from Chapter 2 of Seinfeld and Pandis (2006). These equations form the foundation for understanding the fate of trace species in the atmosphere.

Atmospheric lifetime is a fundamental concept in atmospheric chemistry that characterizes how long a species remains in a given atmospheric reservoir before being removed by chemical reactions, deposition, or transport. Understanding lifetimes is essential for:

- Predicting the accumulation of pollutants and greenhouse gases
- Assessing the response time of the atmosphere to emission changes
- Understanding the spatial scale of pollution (local vs. global)
- Evaluating the effectiveness of emission control strategies

The module implements several key concepts:

| Component | Equations | Description |
|:----------|:----------|:------------|
| [`AtmosphericBudget`](@ref) | Eq. 2.1-2.2 | Mass conservation and steady-state conditions |
| [`SpeciesLifetime`](@ref) | Eq. 2.3-2.6 | Lifetime calculations for various scenarios |
| [`MultipleRemovalLifetime`](@ref) | Eq. 2.7-2.9 | Combined lifetime from parallel removal pathways |
| [`OHReactionLifetime`](@ref) | Eq. 2.12 | Lifetime due to OH radical reaction |
| [`TroposphericBudget`](@ref) | Eq. 2.13-2.17 | Complete tropospheric budget with multiple sources and sinks |

### Reference

> Seinfeld, J. H. and Pandis, S. N.: *Atmospheric Chemistry and Physics: From Air Pollution to Climate Change*, 2nd Edition, Chapter 2: Atmospheric Trace Constituents, John Wiley & Sons, Inc., Hoboken, New Jersey, 2006, ISBN: 978-0-471-72017-1.

---

## Implementation

This section provides detailed documentation of each equation system implemented in the AtmosphericLifetime module.

### AtmosphericBudget (Eq. 2.1-2.2)

The fundamental mass conservation equation for an atmospheric species states that the rate of change of mass in a reservoir equals the sum of transport and chemistry contributions.

**Equation 2.1** (Mass Conservation):
```math
\frac{dQ}{dt} = (F_{in} - F_{out}) + (P - R)
```

where:
- ``Q``: Total mass (or moles) of the species in the reservoir
- ``F_{in}``: Rate of mass inflow from outside the reservoir
- ``F_{out}``: Rate of mass outflow from the reservoir
- ``P``: Rate of chemical production within the reservoir
- ``R``: Rate of chemical removal within the reservoir

**Equation 2.2** (Steady State):

At steady state, ``dQ/dt = 0``, implying:
```math
F_{in} + P = F_{out} + R
```

#### State Variables

```@example budget
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using GasChem.AtmosphericLifetime: AtmosphericBudget

@named sys = AtmosphericBudget()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

#### Equations

```@example budget
equations(sys)
```

---

### SpeciesLifetime (Eq. 2.3-2.6)

The SpeciesLifetime component calculates atmospheric lifetime using various formulations depending on the scenario.

**Equation 2.3** (General Lifetime):
```math
\tau = \frac{Q}{R + F_{out}}
```

This general definition applies to any reservoir with both chemical removal and outflow.

**Equation 2.4** (Global Atmosphere):

For the global atmosphere as a whole, there is no outflow to other reservoirs (``F_{out} = 0``):
```math
\tau = \frac{Q}{R}
```

**Equation 2.5** (Steady State):

At steady state with no external transport (``F_{in} = F_{out} = 0``), production equals removal:
```math
\tau = \frac{Q}{R} = \frac{Q}{P}
```

**Equation 2.6** (First-Order Removal):

When removal follows first-order kinetics (``R = \lambda Q``):
```math
\tau = \frac{1}{\lambda}
```

#### State Variables

```@example lifetime
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using GasChem.AtmosphericLifetime: SpeciesLifetime

@named sys = SpeciesLifetime()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

#### Parameters

```@example lifetime
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [string(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

#### Equations

```@example lifetime
equations(sys)
```

---

### MultipleRemovalLifetime (Eq. 2.7-2.9)

When multiple independent removal pathways operate in parallel, the combined lifetime is calculated from the individual pathway lifetimes.

**Equation 2.7** (Two First-Order Processes):
```math
\tau = \frac{1}{k_1 + k_2}
```

**Equation 2.8** (Inverse Lifetime Sum):
```math
\frac{1}{\tau} = \frac{1}{\tau_1} + \frac{1}{\tau_2}
```

This shows that inverse lifetimes (removal rates) are additive for parallel processes.

**Equation 2.9** (Combined Lifetime Formula):
```math
\tau = \frac{\tau_1 \cdot \tau_2}{\tau_1 + \tau_2}
```

This can be extended to N pathways:
```math
\frac{1}{\tau} = \sum_{i=1}^{N} \frac{1}{\tau_i}
```

#### State Variables

```@example multiple
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using GasChem.AtmosphericLifetime: MultipleRemovalLifetime

@named sys = MultipleRemovalLifetime()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

#### Parameters

```@example multiple
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [string(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

#### Equations

```@example multiple
equations(sys)
```

---

### OHReactionLifetime (Eq. 2.12)

The hydroxyl radical (OH) is the primary oxidant in the troposphere for many trace species. The lifetime due to OH reaction is calculated from the second-order rate constant and the ambient OH concentration.

**Equation 2.12**:
```math
\tau = \frac{1}{k_{OH} \cdot [OH]}
```

where:
- ``k_{OH}``: Second-order rate constant for reaction with OH (cm^3/(molecule s))
- ``[OH]``: Concentration of OH radicals (molecule/cm^3)

The product ``k_{OH} \cdot [OH]`` gives an effective first-order rate constant consistent with Eq. 2.6.

**Typical Values**:
- ``k_{OH}``: varies widely, e.g., ``10^{-14}`` to ``10^{-10}`` cm^3/(molecule s)
- ``[OH]``: ~``10^6`` molecule/cm^3 (global tropospheric average)
- ``\tau``: seconds to years depending on species

#### State Variables

```@example oh
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using GasChem.AtmosphericLifetime: OHReactionLifetime

@named sys = OHReactionLifetime()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

#### Parameters

```@example oh
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [string(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

#### Equations

```@example oh
equations(sys)
```

---

### TroposphericBudget (Eq. 2.13-2.17)

The TroposphericBudget component implements the complete mass balance for a trace species in the troposphere, including multiple source and sink terms.

**Equations 2.13-2.14** (Tropospheric Budget):
```math
\frac{dQ_i}{dt} = P_i^n + P_i^a + P_i^c - (k_i^d + k_i^w + k_i^c + k_i^t) Q_i
```

**Source terms** (production rates):
- ``P_n``: Natural emissions (biogenic, volcanic, oceanic, etc.)
- ``P_a``: Anthropogenic emissions (combustion, industrial, agricultural, etc.)
- ``P_c``: Chemical production within the troposphere

**Removal terms** (first-order rate constants):
- ``k_d``: Dry deposition (uptake by surfaces)
- ``k_w``: Wet deposition (scavenging by precipitation)
- ``k_c``: Chemical loss (reaction with OH, O3, NO3, etc.)
- ``k_t``: Transport to stratosphere (cross-tropopause flux)

**Equation 2.15** (Lifetime from Removal):
```math
\tau_i = \frac{1}{k_d + k_w + k_c + k_t}
```

**Equations 2.16-2.17** (Lifetime from Production at Steady State):
```math
\tau_i = \frac{Q_i}{P_n + P_a + P_c}
```

#### State Variables

```@example trop
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using GasChem.AtmosphericLifetime: TroposphericBudget

@named sys = TroposphericBudget()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

#### Parameters

```@example trop
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [string(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

#### Equations

```@example trop
equations(sys)
```

---

## Analysis

This section demonstrates the physical behavior of the atmospheric lifetime equations through numerical examples and visualizations.

### Time Evolution of Species Mass

The following example shows how a species mass evolves over time in the troposphere, starting from zero and approaching steady state. This demonstrates the fundamental behavior of Eq. 2.13.

```@example time_evolution
using Plots
using OrdinaryDiffEqDefault
using ModelingToolkit
using ModelingToolkit: t
using GasChem.AtmosphericLifetime: TroposphericBudget

# Create the tropospheric budget system
@named trop = TroposphericBudget()

# Prepare for compilation
trop_nns = toggle_namespacing(trop, false)
inputs = [trop_nns.P_n, trop_nns.P_a, trop_nns.P_c]
compiled = mtkcompile(trop; inputs)

# Define parameters for a hypothetical reactive species
# Rate constants (s^-1)
k_d_val = 1e-6    # Dry deposition (~11.6 day lifetime)
k_w_val = 5e-7    # Wet deposition (~23.1 day lifetime)
k_c_val = 2e-6    # Chemical loss (~5.8 day lifetime)
k_t_val = 1e-7    # Stratospheric transport (~115.7 day lifetime)

# Source rates (mol/s)
P_n_val = 1e5     # Natural emissions
P_a_val = 5e4     # Anthropogenic emissions
P_c_val = 2e4     # Chemical production

# Calculate expected steady state and lifetime
P_total_val = P_n_val + P_a_val + P_c_val
k_total_val = k_d_val + k_w_val + k_c_val + k_t_val
Q_ss = P_total_val / k_total_val
tau = 1 / k_total_val
tau_days = tau / 86400

# Set up and solve ODE
Q_0 = 0.0  # Start with zero mass
t_end = 10 * tau  # Simulate for 10 lifetimes

prob = ODEProblem(
    compiled,
    Dict(
        compiled.Q => Q_0,
        compiled.k_d => k_d_val,
        compiled.k_w => k_w_val,
        compiled.k_c => k_c_val,
        compiled.k_t => k_t_val,
        compiled.P_n => P_n_val,
        compiled.P_a => P_a_val,
        compiled.P_c => P_c_val,
    ),
    (0.0, t_end)
)

sol = solve(prob, saveat=tau/10)

# Convert time to days for plotting
t_days = sol.t ./ 86400
Q_values = sol[compiled.Q]

# Create the plot
p1 = plot(
    t_days, Q_values ./ 1e10,
    xlabel = "Time (days)",
    ylabel = "Species Mass (10^10 mol)",
    title = "Approach to Steady State in the Troposphere",
    label = "Q(t)",
    linewidth = 2,
    color = :blue,
    legend = :bottomright,
    size = (700, 450),
    grid = true
)

# Add steady state line
hline!(p1, [Q_ss / 1e10],
    label = "Steady State (Q_ss)",
    linestyle = :dash,
    color = :red,
    linewidth = 2)

# Add e-folding time markers
vline!(p1, [tau_days],
    label = "tau = $(round(tau_days, digits=1)) days",
    linestyle = :dot,
    color = :green,
    linewidth = 2)

# Add annotation
annotate!(p1, tau_days * 5, Q_ss * 0.5 / 1e10,
    text("63% of Q_ss reached\nat t = tau", 10, :left))

savefig(p1, "time_evolution.png")
p1
```

**Physical Interpretation**: The species mass grows exponentially at first, then levels off as removal processes balance the sources. The characteristic timescale for this approach to equilibrium is the atmospheric lifetime ``\tau``. After one lifetime, the mass reaches 63% (``1 - 1/e``) of its steady-state value.

---

### Exponential Decay from Initial Burden

When sources are removed, a species decays exponentially with the characteristic lifetime. This demonstrates Eq. 2.6.

```@example decay
using Plots
using OrdinaryDiffEqDefault
using ModelingToolkit
using ModelingToolkit: t
using GasChem.AtmosphericLifetime: TroposphericBudget

@named trop = TroposphericBudget()

trop_nns = toggle_namespacing(trop, false)
inputs = [trop_nns.P_n, trop_nns.P_a, trop_nns.P_c]
compiled = mtkcompile(trop; inputs)

# Pure decay scenario (no sources)
lambda = 1e-5  # s^-1, corresponds to tau = 10^5 s ~ 1.16 days
Q_0 = 1e12     # Initial burden (mol)

prob = ODEProblem(
    compiled,
    Dict(
        compiled.Q => Q_0,
        compiled.k_d => lambda,
        compiled.k_w => 0.0,
        compiled.k_c => 0.0,
        compiled.k_t => 0.0,
        compiled.P_n => 0.0,
        compiled.P_a => 0.0,
        compiled.P_c => 0.0,
    ),
    (0.0, 5 / lambda)
)

sol = solve(prob, saveat=(1/lambda)/10)

# Calculate analytical solution for comparison
tau = 1 / lambda
t_values = sol.t
Q_analytical = Q_0 .* exp.(-t_values ./ tau)

# Convert to normalized units
t_normalized = t_values ./ tau
Q_normalized = sol[compiled.Q] ./ Q_0

p2 = plot(
    t_normalized, Q_normalized,
    xlabel = "Time (t/tau)",
    ylabel = "Normalized Mass (Q/Q_0)",
    title = "Exponential Decay of Atmospheric Species",
    label = "Numerical Solution",
    linewidth = 3,
    color = :blue,
    legend = :topright,
    size = (700, 450),
    grid = true
)

# Add analytical solution
plot!(p2, t_normalized, Q_analytical ./ Q_0,
    label = "Analytical: exp(-t/tau)",
    linestyle = :dash,
    linewidth = 2,
    color = :red)

# Mark e-folding times
for n in 1:4
    vline!(p2, [n], label = n == 1 ? "e-folding times" : "",
        linestyle = :dot, color = :gray, linewidth = 1)
    Q_at_n = exp(-n)
    scatter!(p2, [n], [Q_at_n],
        label = "", markersize = 6, color = :green)
    annotate!(p2, n + 0.1, Q_at_n + 0.05,
        text("$(round(Q_at_n*100, digits=1))%", 8, :left))
end

savefig(p2, "exponential_decay.png")
p2
```

**Key Results**:
- After 1 lifetime: 36.8% remains (``e^{-1}``)
- After 2 lifetimes: 13.5% remains (``e^{-2}``)
- After 3 lifetimes: 5.0% remains (``e^{-3}``)
- After 5 lifetimes: 0.7% remains (``e^{-5}``)

---

### Effect of Multiple Removal Pathways

This analysis demonstrates how parallel removal pathways combine to reduce the total atmospheric lifetime (Eq. 2.7-2.9).

```@example multiple_pathways
using Plots

# Compare combined lifetime to individual pathway lifetimes
tau_1_range = 1:1:30  # days
tau_2_fixed = 10.0     # days

# Calculate combined lifetimes using Eq. 2.9
tau_combined = [(t1 * tau_2_fixed) / (t1 + tau_2_fixed) for t1 in tau_1_range]

p3 = plot(
    tau_1_range, tau_combined,
    xlabel = "Lifetime from Pathway 1, tau_1 (days)",
    ylabel = "Combined Lifetime, tau (days)",
    title = "Combined Lifetime from Two Removal Pathways (tau_2 = 10 days)",
    label = "tau_combined (Eq. 2.9)",
    linewidth = 3,
    color = :blue,
    legend = :bottomright,
    size = (700, 450),
    grid = true,
    ylims = (0, 12)
)

# Add reference lines
hline!(p3, [tau_2_fixed], label = "tau_2 = 10 days",
    linestyle = :dash, color = :red, linewidth = 2)
plot!(p3, tau_1_range, collect(tau_1_range), label = "tau_1",
    linestyle = :dot, color = :green, linewidth = 2)

# Mark special cases
scatter!(p3, [tau_2_fixed], [(tau_2_fixed * tau_2_fixed) / (2 * tau_2_fixed)],
    label = "Equal lifetimes: tau = tau_1/2",
    markersize = 8, color = :orange)

# Add annotation
annotate!(p3, 20, 3,
    text("Combined lifetime is always\nless than either individual lifetime", 10, :left))

savefig(p3, "multiple_pathways.png")
p3
```

**Physical Interpretation**:
- When ``\tau_1 = \tau_2``, the combined lifetime is exactly half of either individual lifetime
- The combined lifetime is always less than the smaller of the two individual lifetimes
- A fast removal pathway dominates the total lifetime even if a slow pathway exists

---

### Contribution of Removal Pathways to Total Lifetime

This figure shows how each removal pathway contributes to the total inverse lifetime for a typical tropospheric species.

```@example pathway_contributions
using Plots

# Typical rate constants for a reactive species (s^-1)
# Convert to day^-1 for visualization
k_d = 1e-6 * 86400   # Dry deposition
k_w = 5e-7 * 86400   # Wet deposition
k_c = 2e-6 * 86400   # Chemical loss
k_t = 1e-7 * 86400   # Stratospheric transport

k_total = k_d + k_w + k_c + k_t

# Calculate individual lifetimes (days)
tau_dry = 1 / k_d
tau_wet = 1 / k_w
tau_chem = 1 / k_c
tau_transport = 1 / k_t
tau_total = 1 / k_total

# Calculate fractional contributions to 1/tau
frac_dry = k_d / k_total * 100
frac_wet = k_w / k_total * 100
frac_chem = k_c / k_total * 100
frac_transport = k_t / k_total * 100

# Create bar chart of contributions
labels = ["Dry Dep.", "Wet Dep.", "Chemical", "Transport"]
fractions = [frac_dry, frac_wet, frac_chem, frac_transport]
lifetimes = [tau_dry, tau_wet, tau_chem, tau_transport]
colors = [:brown, :blue, :orange, :purple]

p4 = bar(
    labels, fractions,
    xlabel = "",
    ylabel = "Contribution to 1/tau (%)",
    title = "Removal Pathway Contributions\n(Total lifetime = $(round(tau_total, digits=1)) days)",
    label = "",
    color = colors,
    size = (700, 450),
    ylims = (0, 60),
    bar_width = 0.6
)

# Add lifetime labels on bars
for (i, (frac, tau)) in enumerate(zip(fractions, lifetimes))
    annotate!(p4, i, frac + 2,
        text("tau = $(round(tau, digits=1)) d", 9, :center))
end

savefig(p4, "pathway_contributions.png")
p4
```

**Key Observations**:
- Chemical loss typically dominates for reactive species
- Dry deposition is significant for species that deposit readily to surfaces
- Wet deposition is important for soluble species
- Stratospheric transport is usually a minor sink for most tropospheric species

---

### OH Reaction Lifetimes for Common Atmospheric Species

This analysis compares the OH-reaction lifetimes for several important atmospheric trace gases, demonstrating Eq. 2.12.

```@example oh_lifetimes
using Plots, DataFrames

# OH rate constants at 298 K (cm^3 molecule^-1 s^-1) from JPL evaluation
# and typical global mean OH concentration
OH_conc = 1e6  # molecule/cm^3

species_data = [
    ("Methane (CH4)", 6.3e-15, "Years"),
    ("Ethane (C2H6)", 2.4e-13, "Months"),
    ("Propane (C3H8)", 1.1e-12, "Days"),
    ("Carbon Monoxide (CO)", 2.4e-13, "Months"),
    ("Formaldehyde (HCHO)", 8.5e-12, "Hours"),
    ("Isoprene (C5H8)", 1.0e-10, "Hours"),
    ("Nitrogen Dioxide (NO2)", 1.2e-11, "Hours"),
]

# Calculate lifetimes
species_names = String[]
k_OH_values = Float64[]
tau_values = Float64[]
tau_units = String[]

for (name, k_OH, unit) in species_data
    tau = 1 / (k_OH * OH_conc)  # seconds
    push!(species_names, name)
    push!(k_OH_values, k_OH)
    push!(tau_units, unit)

    # Convert to appropriate units
    if unit == "Years"
        push!(tau_values, tau / (365.25 * 86400))
    elseif unit == "Months"
        push!(tau_values, tau / (30 * 86400))
    elseif unit == "Days"
        push!(tau_values, tau / 86400)
    else  # Hours
        push!(tau_values, tau / 3600)
    end
end

# Create summary table
df = DataFrame(
    Species = species_names,
    k_OH_cm3_molec_s = k_OH_values,
    Lifetime = tau_values,
    Units = tau_units
)
df
```

```@example oh_lifetimes
# Create bar chart of lifetimes (all converted to days for comparison)
tau_days = Float64[]
for (name, k_OH, _) in species_data
    tau = 1 / (k_OH * OH_conc)  # seconds
    push!(tau_days, tau / 86400)
end

p5 = bar(
    1:length(species_names),
    log10.(tau_days),
    xlabel = "",
    ylabel = "log10(Lifetime in days)",
    title = "OH Reaction Lifetimes at [OH] = 10^6 molecule/cm^3",
    label = "",
    color = :steelblue,
    size = (800, 500),
    xticks = (1:length(species_names),
        ["CH4", "C2H6", "C3H8", "CO", "HCHO", "Isoprene", "NO2"]),
    bar_width = 0.6
)

# Add horizontal lines for reference timescales
hline!(p5, [log10(365)], label = "1 year", linestyle = :dash, color = :red, linewidth = 2)
hline!(p5, [log10(30)], label = "1 month", linestyle = :dash, color = :orange, linewidth = 2)
hline!(p5, [log10(1)], label = "1 day", linestyle = :dash, color = :green, linewidth = 2)
hline!(p5, [log10(1/24)], label = "1 hour", linestyle = :dash, color = :purple, linewidth = 2)

savefig(p5, "oh_lifetimes.png")
p5
```

**Physical Interpretation**:
- **Long-lived species** (CH4, CO): Global distribution, well-mixed in the troposphere
- **Intermediate species** (C2H6, C3H8): Regional to hemispheric influence
- **Short-lived species** (HCHO, isoprene, NO2): Local influence, highly variable concentrations

---

### Sensitivity to OH Concentration

The effective lifetime depends strongly on the ambient OH concentration, which varies with location, season, and time of day.

```@example oh_sensitivity
using Plots

# OH concentration range (molecule/cm^3)
OH_range = 10 .^ range(5, 7, length=50)

# k_OH for methane and CO
k_OH_CH4 = 6.3e-15  # cm^3/(molecule s)
k_OH_CO = 2.4e-13   # cm^3/(molecule s)

# Calculate lifetimes
tau_CH4_years = @. 1 / (k_OH_CH4 * OH_range) / (365.25 * 86400)
tau_CO_days = @. 1 / (k_OH_CO * OH_range) / 86400

p6 = plot(
    OH_range, tau_CH4_years,
    xlabel = "[OH] (molecule/cm^3)",
    ylabel = "Methane Lifetime (years)",
    title = "Dependence of Lifetime on OH Concentration",
    label = "CH4",
    linewidth = 2,
    color = :blue,
    xscale = :log10,
    yscale = :log10,
    legend = :topright,
    size = (700, 500),
    grid = true
)

# Create second y-axis plot
p6b = twinx(p6)
plot!(p6b, OH_range, tau_CO_days,
    ylabel = "CO Lifetime (days)",
    label = "CO",
    linewidth = 2,
    color = :red,
    xscale = :log10,
    yscale = :log10,
    legend = :bottomleft)

# Mark typical global mean OH
vline!(p6, [1e6], label = "Global mean [OH]",
    linestyle = :dash, color = :gray, linewidth = 2)

savefig(p6, "oh_sensitivity.png")
p6
```

**Physical Interpretation**:
- Lifetime scales inversely with OH concentration (``\tau \propto 1/[OH]``)
- A factor of 10 increase in OH leads to a factor of 10 decrease in lifetime
- OH concentrations vary from ~10^5 molecule/cm^3 (remote marine) to ~10^7 molecule/cm^3 (polluted urban)

---

### Methane Lifetime Budget

Methane has multiple removal pathways, making it a good example of applying Eq. 2.7-2.9. The total atmospheric lifetime combines:
- Tropospheric OH oxidation (~10 years)
- Soil uptake (~150 years)
- Stratospheric loss (~120 years)

```@example ch4_budget
using Plots

# Methane removal pathways (lifetimes in years)
tau_OH = 10.0       # Tropospheric OH reaction
tau_soil = 150.0    # Soil uptake
tau_strat = 120.0   # Stratospheric loss

# Calculate combined lifetime using Eq. 2.8
inverse_tau_total = 1/tau_OH + 1/tau_soil + 1/tau_strat
tau_total = 1 / inverse_tau_total

# Calculate fractional contributions
frac_OH = (1/tau_OH) / inverse_tau_total * 100
frac_soil = (1/tau_soil) / inverse_tau_total * 100
frac_strat = (1/tau_strat) / inverse_tau_total * 100

println("Methane Atmospheric Lifetime Analysis")
println("=" ^ 40)
println("Individual pathway lifetimes:")
println("  OH reaction:     $(tau_OH) years")
println("  Soil uptake:     $(tau_soil) years")
println("  Stratospheric:   $(tau_strat) years")
println()
println("Combined total lifetime: $(round(tau_total, digits=1)) years")
println()
println("Removal pathway contributions:")
println("  OH reaction:     $(round(frac_OH, digits=1))%")
println("  Soil uptake:     $(round(frac_soil, digits=1))%")
println("  Stratospheric:   $(round(frac_strat, digits=1))%")

# Create pie chart of removal contributions
labels_pie = ["OH Reaction\n($(round(frac_OH, digits=1))%)",
          "Soil Uptake\n($(round(frac_soil, digits=1))%)",
          "Stratospheric\n($(round(frac_strat, digits=1))%)"]
fracs = [frac_OH, frac_soil, frac_strat]

p7 = pie(
    fracs,
    label = labels_pie,
    title = "Methane Removal Budget\n(Total lifetime = $(round(tau_total, digits=1)) years)",
    color = [:orange, :brown, :purple],
    size = (600, 500)
)

savefig(p7, "ch4_budget.png")
p7
```

**Key Result**: Despite individual pathway lifetimes ranging from 10 to 150 years, the combined atmospheric lifetime of methane is approximately 8-9 years. This is because the fastest removal process (OH oxidation) dominates the total removal rate.

---

## Summary

The AtmosphericLifetime module provides a comprehensive implementation of the fundamental atmospheric lifetime equations from Seinfeld and Pandis Chapter 2. Key takeaways:

1. **Mass Conservation** (Eq. 2.1): The rate of change of atmospheric burden equals sources minus sinks

2. **Lifetime Definition** (Eq. 2.3-2.6): Lifetime characterizes the persistence of a species in the atmosphere

3. **Multiple Pathways** (Eq. 2.7-2.9): Parallel removal processes combine via their inverse lifetimes; the fastest process dominates

4. **OH Oxidation** (Eq. 2.12): The OH radical is the primary oxidant determining the lifetime of most reduced species

5. **Tropospheric Budget** (Eq. 2.13-2.17): A complete budget includes natural and anthropogenic emissions, chemical production, and multiple removal pathways

These concepts are essential for understanding air quality, climate forcing by greenhouse gases, and the design of emission control strategies.

---

## API Reference

```@docs
AtmosphericBudget
SpeciesLifetime
MultipleRemovalLifetime
OHReactionLifetime
TroposphericBudget
```
