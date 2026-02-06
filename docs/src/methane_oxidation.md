# [Methane Oxidation](@id methane_oxidation)

## Overview

Methane (CH4) is the most abundant hydrocarbon in the atmosphere (~1.8 ppm).
Its oxidation is the archetypal VOC oxidation mechanism, proceeding through
several intermediates:

    CH4 -> CH3O2 -> CH3O -> HCHO -> HCO -> CO -> CO2

Each step can produce ozone when NOx is present through peroxy + NO reactions.
At high NOx, complete CH4 oxidation can produce 3-5 O3 molecules.

Two components are provided:

  - `MethaneOxidation`: Algebraic system computing individual reaction rates and diagnostics from Table 6.1
  - `MethaneOxidationODE`: Full ODE system for time evolution of all species

**Reference**: Seinfeld, J.H. and Pandis, S.N. (2006). *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change*, 2nd Edition. John Wiley & Sons. Section 6.4, Table 6.1, pp. 219-227.

```@docs
MethaneOxidation
```

```@docs
MethaneOxidationODE
```

## Implementation

### MethaneOxidation: State Variables

```@example ch4_ox
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities, GasChem

sys = MethaneOxidation()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### MethaneOxidation: Parameters

```@example ch4_ox
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### MethaneOxidation: Equations

```@example ch4_ox
eqs = equations(sys)
```

### MethaneOxidationODE: State Variables

```@example ch4_ox
ode_sys = MethaneOxidationODE()
vars_ode = unknowns(ode_sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars_ode],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars_ode],
    :Description => [ModelingToolkit.getdescription(v) for v in vars_ode]
)
```

### Table 6.1: Rate Constants

Table 6.1 of Seinfeld & Pandis lists the 17 reactions of the methane oxidation
mechanism with their rate constants at 298 K. The table below reproduces
these values from the implementation parameters (converted from SI back to
cm³ molecule⁻¹ s⁻¹ for comparison with the textbook):

```@example ch4_ox
using DataFrames

# Reproduce Table 6.1 from the model parameters
reactions = [
    "1. CH₄ + OH → CH₃ + H₂O",
    "2. CH₃ + O₂ + M → CH₃O₂ + M",
    "3. CH₃O₂ + NO → CH₃O + NO₂",
    "4. CH₃O₂ + HO₂ → CH₃OOH + O₂",
    "5. CH₃O₂ + CH₃O₂ → products",
    "6. CH₃O + O₂ → HCHO + HO₂",
    "7. CH₃OOH + OH → CH₃O₂ + H₂O",
    "8. CH₃OOH + OH → HCHO + OH + H₂O",
    "9. CH₃OOH + hν → CH₃O + OH",
    "10. HCHO + OH → HCO + H₂O",
    "11. HCHO + hν → HCO + H",
    "12. HCHO + hν → H₂ + CO",
    "13. HCO + O₂ → CO + HO₂",
    "14. H + O₂ + M → HO₂ + M",
    "15. HO₂ + NO → OH + NO₂",
    "16. NO₂ + hν → NO + O",
    "17. O + O₂ + M → O₃ + M"
]

# Rate constants from Table 6.1 at 298 K
k_values = [
    "6.3 × 10⁻¹⁵",
    "1.0 × 10⁻³⁰ [M]",
    "7.7 × 10⁻¹²",
    "5.2 × 10⁻¹²",
    "3.5 × 10⁻¹³",
    "1.9 × 10⁻¹⁵",
    "3.8 × 10⁻¹²",
    "1.9 × 10⁻¹²",
    "j ≈ 5 × 10⁻⁶ s⁻¹",
    "9.0 × 10⁻¹²",
    "j ≈ 3 × 10⁻⁵ s⁻¹",
    "j ≈ 4 × 10⁻⁵ s⁻¹",
    "5.2 × 10⁻¹²",
    "5.7 × 10⁻³² [M]",
    "8.1 × 10⁻¹²",
    "j ≈ 8 × 10⁻³ s⁻¹",
    "6.0 × 10⁻³⁴ [M]"
]

units = [
    "cm³ molec⁻¹ s⁻¹",
    "cm⁶ molec⁻² s⁻¹",
    "cm³ molec⁻¹ s⁻¹",
    "cm³ molec⁻¹ s⁻¹",
    "cm³ molec⁻¹ s⁻¹",
    "cm³ molec⁻¹ s⁻¹",
    "cm³ molec⁻¹ s⁻¹",
    "cm³ molec⁻¹ s⁻¹",
    "s⁻¹",
    "cm³ molec⁻¹ s⁻¹",
    "s⁻¹",
    "s⁻¹",
    "cm³ molec⁻¹ s⁻¹",
    "cm⁶ molec⁻² s⁻¹",
    "cm³ molec⁻¹ s⁻¹",
    "s⁻¹",
    "cm⁶ molec⁻² s⁻¹"
]

DataFrame(:Reaction => reactions, :k_298K => k_values, :Units => units)
```

## Analysis

### Methane Oxidation Chain: O3 Yield at High vs Low NOx

The fate of the methylperoxy radical (CH3O2) depends on the NOx level.
At high NOx, CH3O2 reacts with NO to produce NO2 (and subsequently O3),
while at low NOx, CH3O2 reacts with HO2 to form CH3OOH, which terminates
the radical chain. This determines the overall O3 yield per CH4 molecule
oxidized.

This analysis uses the `MethaneOxidation` system to compute the competing
reaction rates R3 (CH3O2 + NO) and R4 (CH3O2 + HO2) across a range of
NO concentrations.

```@example ch4_ox
using Plots, NonlinearSolve

sys_nns = ModelingToolkit.toggle_namespacing(sys, false)
input_vars = [sys_nns.CH4, sys_nns.CH3, sys_nns.CH3O2, sys_nns.CH3O, sys_nns.CH3OOH,
    sys_nns.HCHO, sys_nns.HCO, sys_nns.OH, sys_nns.HO2, sys_nns.H,
    sys_nns.NO, sys_nns.NO2, sys_nns.O, sys_nns.O2, sys_nns.M]
compiled = mtkcompile(sys; inputs = input_vars)

# Fixed conditions (SI: m⁻³)
M_val = 2.5e25
O2_val = 5.25e24
HO2_val = 1e14     # m⁻³ (= 1e8 cm⁻³)
CH3O2_val = 1e14   # m⁻³
CH4_val = 4.5e19   # ~1800 ppb
OH_val = 1e12      # typical daytime

# Vary NO from 10 ppt to 100 ppb (m⁻³)
NO_ppb = 10 .^ range(-2, 2, length = 200)
NO_vals = NO_ppb .* 2.5e16  # m⁻³

# Compute R3 (CH3O2+NO) and R4 (CH3O2+HO2) using the system
R3_vals = Float64[]
R4_vals = Float64[]

# Set all species to typical values
base_dict = Dict(
    compiled.CH4 => CH4_val, compiled.CH3 => 1e10, compiled.CH3O2 => CH3O2_val,
    compiled.CH3O => 1e10, compiled.CH3OOH => 1e14, compiled.HCHO => 1e15,
    compiled.HCO => 1e10, compiled.OH => OH_val, compiled.HO2 => HO2_val,
    compiled.H => 1e8, compiled.NO => NO_vals[1], compiled.NO2 => 1e16,
    compiled.O => 1e8, compiled.O2 => O2_val, compiled.M => M_val)

prob = NonlinearProblem(compiled, base_dict; build_initializeprob = false)

for no in NO_vals
    newprob = remake(prob, p = [compiled.NO => no])
    sol = solve(newprob)
    push!(R3_vals, sol[compiled.R3])
    push!(R4_vals, sol[compiled.R4])
end

# Fraction through NO pathway
f_NO = R3_vals ./ (R3_vals .+ R4_vals)

# Approximate O3 yield (up to ~4 O3 per CH4 at full NO pathway)
O3_yield = 4.0 .* f_NO

p1 = plot(NO_ppb, f_NO .* 100,
    xlabel = "NO (ppb)", ylabel = "Fraction via NO pathway (%)",
    title = "CH₃O₂ Fate: NO vs HO₂",
    xscale = :log10, linewidth = 2,
    label = "% through CH₃O₂ + NO",
    legend = :bottomright)

p2 = plot(NO_ppb, O3_yield,
    xlabel = "NO (ppb)", ylabel = "O₃ molecules per CH₄",
    title = "O₃ Yield from CH₄ Oxidation",
    xscale = :log10, linewidth = 2,
    label = "Approximate O₃ yield",
    legend = :bottomright, ylims = (0, 5))

plot(p1, p2, layout = (1, 2), size = (800, 350))
savefig("ch4_o3_yield.svg") # hide
```

![CH3O2 fate and O3 yield vs NOx](ch4_o3_yield.svg)

The left panel shows the fraction of CH3O2 that reacts with NO (the O3-producing
pathway) vs HO2 (the chain-terminating pathway). At 1 ppb NO, nearly all
CH3O2 reacts with NO. The right panel shows that the approximate O3 yield per
CH4 molecule increases from near zero at very low NOx (where peroxide
formation dominates) to about 4 at high NOx, consistent with the discussion
in Section 6.4 of Seinfeld & Pandis.

### Formaldehyde Branching Pathways

Formaldehyde (HCHO) is a key intermediate with three loss pathways:
reaction with OH, radical photolysis (HCO + H), and molecular photolysis
(H2 + CO). The branching ratio affects the HOx budget.

This analysis uses the `MethaneOxidation` system to compute the three HCHO
loss rates (R10, R11, R12) at typical daytime conditions.

```@example ch4_ox
# Compute HCHO loss rates using the system
HCHO_val = 1e15   # m⁻³
OH_val_hcho = 1e12  # m⁻³ (= 1e6 cm⁻³)

hcho_prob = remake(prob, p = [
    compiled.HCHO => HCHO_val, compiled.OH => OH_val_hcho])
sol = solve(hcho_prob)

R10_val = sol[compiled.R10]  # HCHO + OH
R11_val = sol[compiled.R11]  # HCHO → HCO + H
R12_val = sol[compiled.R12]  # HCHO → H₂ + CO

L_total = R10_val + R11_val + R12_val
fractions = [R10_val / L_total * 100, R11_val / L_total * 100, R12_val / L_total * 100]
labels_bar = [
    "HCHO + OH\n(radical)", "HCHO + hν → HCO + H\n(radical)", "HCHO + hν → H₂ + CO\n(molecular)"]

bar(labels_bar, fractions,
    ylabel = "Fraction of HCHO Loss (%)",
    title = "Formaldehyde Loss Pathways",
    label = false, color = [:steelblue :orange :green],
    size = (600, 400), ylims = (0, 100))
savefig("hcho_branching.svg") # hide
```

![HCHO loss pathway branching](hcho_branching.svg)

At typical daytime conditions (OH = 10^6 molec/cm3), the molecular photolysis
channel (producing H2 + CO) is the largest single loss pathway. Both the OH
reaction and radical photolysis channels produce HOx radicals that contribute
to further O3 production, while the molecular channel does not produce radicals.
