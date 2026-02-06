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
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### COOxidation: Parameters

```@example co_ox
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
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
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars_ope],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars_ope],
    :Description => [ModelingToolkit.getdescription(v) for v in vars_ope]
)
```

### OzoneProductionEfficiency: Equations

```@example co_ox
equations(ope_sys)
```

## Analysis

### Figure 6.3: OPE vs NOx for CO Oxidation

Figure 6.3 of Seinfeld & Pandis shows the ozone production efficiency (OPE)
for atmospheric CO oxidation at 298 K at the Earth's surface as a function
of the NOx (NO + NO2) level. The conditions are: ``P_{HO_x} = 1`` ppt s⁻¹,
``[NO]/[NO_2] = 0.1``, and CO mixing ratio of 200 ppb.

OPE is computed from Eqs. 6.10, 6.23, and 6.22 by solving the quadratic
equation for [HO2] (Eq. 6.23 combined with 6.10) and then computing
OPE = P(O3) / L(NOx). The `COOxidation` and `OzoneProductionEfficiency`
systems are used to compute the diagnostics.

```@example co_ox
using Plots, NonlinearSolve

# Use the OPE system for computing OPE from given concentrations
ope_nns = ModelingToolkit.toggle_namespacing(ope_sys, false)
ope_inputs = [ope_nns.OH, ope_nns.HO2, ope_nns.RO2, ope_nns.NO, ope_nns.NO2]
ope_compiled = mtkcompile(ope_sys; inputs = ope_inputs)

# Rate constants at 298 K (cm³ molecule⁻¹ s⁻¹) for the quadratic solve
# (These are the same values as in the COOxidation system parameters)
k_CO_OH = 2.4e-13    # CO + OH
k_HO2_NO = 8.1e-12   # HO₂ + NO
k_HO2_HO2 = 2.9e-12  # HO₂ + HO₂
k_OH_NO2 = 1.0e-11   # OH + NO₂

# Conditions from Figure 6.3 caption
M = 2.5e19             # total air at surface [molec/cm³]
CO_val = 200e-9 * M    # 200 ppb CO [molec/cm³]
NO_NO2_ratio = 0.1     # [NO]/[NO₂] = 0.1
P_HOx = 1e-12 * M      # 1 ppt/s [molec/cm³/s]

# Vary NOx from 1 ppt to 10⁶ ppt (= 1 ppm)
NOx_ppt = 10 .^ range(0, 6, length = 300)
NOx = NOx_ppt .* 1e-12 .* M

# Partition NOx
NO2_vals = NOx ./ (1 + NO_NO2_ratio)
NO_vals = NO_NO2_ratio .* NO2_vals

# Solve quadratic for [HO2] from Eqs. 6.10 and 6.23
a_vals = 2 .* k_HO2_HO2 .* (1 .+ k_OH_NO2 .* NO2_vals ./ (k_CO_OH .* CO_val))
b_vals = k_HO2_NO .* k_OH_NO2 .* NO2_vals .* NO_vals ./ (k_CO_OH .* CO_val)
c_val = -P_HOx

HO2_vals = (-b_vals .+ sqrt.(b_vals .^ 2 .- 4 .* a_vals .* c_val)) ./ (2 .* a_vals)
OH_vals = (P_HOx .- 2 .* k_HO2_HO2 .* HO2_vals .^ 2) ./ (k_OH_NO2 .* NO2_vals)

# Compute OPE using the OPE system (converting cgs to SI: × 1e6)
OPE_vals = Float64[]
prob = NonlinearProblem(ope_compiled,
    Dict(ope_compiled.OH => OH_vals[1] * 1e6, ope_compiled.HO2 => HO2_vals[1] * 1e6,
         ope_compiled.RO2 => 0.0, ope_compiled.NO => NO_vals[1] * 1e6,
         ope_compiled.NO2 => NO2_vals[1] * 1e6);
    build_initializeprob = false)

for i in eachindex(NOx_ppt)
    newprob = remake(prob, p = [
        ope_compiled.OH => OH_vals[i] * 1e6, ope_compiled.HO2 => HO2_vals[i] * 1e6,
        ope_compiled.NO => NO_vals[i] * 1e6, ope_compiled.NO2 => NO2_vals[i] * 1e6])
    sol = solve(newprob)
    push!(OPE_vals, sol[ope_compiled.OPE])
end

plot(NOx_ppt, OPE_vals,
    xlabel = "NOₓ (ppt)",
    ylabel = "Ozone Production Efficiency",
    title = "Figure 6.3: OPE for CO Oxidation",
    xscale = :log10,
    linewidth = 2, label = "OPE (Eqs. 6.10, 6.22, 6.23)",
    legend = :topright,
    ylims = (0, 9),
    size = (600, 400))
annotate!([(
    1e4, 7, text("P_HOx = 1 ppt/s\n[NO]/[NO₂] = 0.1\nCO = 200 ppb\nT = 298 K", 8, :left))])
savefig("co_fig6_3.svg") # hide
```

![Figure 6.3: OPE vs NOx](co_fig6_3.svg)

The OPE is largest at the lowest NOx concentrations; at these low levels,
NOx termination by OH + NO₂ is suppressed and each NOx molecule participates
in more O₃ production cycles. At 100 ppb NOx, OPE approaches zero as the
OH + NO₂ reaction occurs preferentially relative to propagation of the cycle.
This reproduces Figure 6.3 of Seinfeld & Pandis.

### Figure 6.4: CO Oxidation Characteristics vs NO

Figure 6.4 shows [HO₂], P(O₃), and the HOx loss terms HHL and NHL as functions
of [NO] for three values of ``P_{HO_x}`` (0.1, 0.6, and 1.2 ppt s⁻¹).
Conditions: 298 K, ``[NO_2]/[NO] = 7``.

The `COOxidation` system is used to compute P(O₃) and the HOx loss diagnostics
from the computed OH and HO₂ concentrations.

  - HHL = ``2 k_{HO_2+HO_2} [HO_2]^2`` (HOx loss via HO₂ self-reaction, Eq. 6.11)
  - NHL = ``k_{OH+NO_2} [OH] [NO_2]`` (HOx loss via OH+NO₂, Eq. 6.12)

```@example co_ox
sys_nns = ModelingToolkit.toggle_namespacing(sys, false)
co_inputs = [sys_nns.CO, sys_nns.OH, sys_nns.HO2, sys_nns.NO, sys_nns.NO2, sys_nns.O3]
co_compiled = mtkcompile(sys; inputs = co_inputs)

# Conditions from Figure 6.4 caption
M = 2.5e19
NO2_NO_ratio = 7.0
P_HOx_ppt = [0.1, 0.6, 1.2]
P_HOx_mks = P_HOx_ppt .* 1e-12 .* M  # molec/cm³/s

# Vary NO from 0 to 6e10 molec/cm³
NO_range = range(1e8, 6e10, length = 500)
NO2_range = NO2_NO_ratio .* NO_range

# CO = 200 ppb
CO_cgs = 200e-9 * M

p_ho2 = plot(title = "(a) [HO₂]", xlabel = "NO (molec cm⁻³)", ylabel = "HO₂ (molec cm⁻³)")
p_po3 = plot(title = "(b) P(O₃)", xlabel = "NO (molec cm⁻³)", ylabel = "P_O₃ (molec cm⁻³ s⁻¹)")
p_hhl = plot(title = "(c) HHL and NHL", xlabel = "NO (molec cm⁻³)",
    ylabel = "HHL, NHL (molec cm⁻³ s⁻¹)")

for (i, P_HOx) in enumerate(P_HOx_mks)
    lbl = "P_HOx = $(P_HOx_ppt[i]) ppt/s"

    # Solve quadratic for HO2 (cgs units)
    a_v = 2 .* k_HO2_HO2 .* (1 .+ k_OH_NO2 .* NO2_range ./ (k_CO_OH .* CO_cgs))
    b_v = k_HO2_NO .* k_OH_NO2 .* NO2_range .* NO_range ./ (k_CO_OH .* CO_cgs)
    c_v = -P_HOx

    HO2_v = (-b_v .+ sqrt.(b_v .^ 2 .- 4 .* a_v .* c_v)) ./ (2 .* a_v)
    OH_v = (P_HOx .- 2 .* k_HO2_HO2 .* HO2_v .^ 2) ./ (k_OH_NO2 .* NO2_range)

    # Use COOxidation system to compute diagnostics (SI units)
    PO3_v = Float64[]
    HHL_v = Float64[]
    NHL_v = Float64[]

    co_prob = NonlinearProblem(co_compiled,
        Dict(co_compiled.CO => CO_cgs * 1e6, co_compiled.OH => OH_v[1] * 1e6,
             co_compiled.HO2 => HO2_v[1] * 1e6, co_compiled.NO => NO_range[1] * 1e6,
             co_compiled.NO2 => NO2_range[1] * 1e6, co_compiled.O3 => 1e18);
        build_initializeprob = false)

    for j in eachindex(NO_range)
        newprob = remake(co_prob, p = [
            co_compiled.OH => OH_v[j] * 1e6, co_compiled.HO2 => HO2_v[j] * 1e6,
            co_compiled.NO => NO_range[j] * 1e6, co_compiled.NO2 => NO2_range[j] * 1e6])
        sol = solve(newprob)
        # Convert back to cgs for plotting (m⁻³ s⁻¹ → cm⁻³ s⁻¹ = ×1e-6)
        push!(PO3_v, sol[co_compiled.P_O3] * 1e-6)
        push!(HHL_v, (sol[co_compiled.L_HOx] - k_OH_NO2 * 1e-6 * OH_v[j] * 1e6 * NO2_range[j] * 1e6) * 1e-6)
        push!(NHL_v, k_OH_NO2 * OH_v[j] * NO2_range[j])
    end

    plot!(p_ho2, NO_range, HO2_v, label = lbl, linewidth = 2)
    plot!(p_po3, NO_range, PO3_v, label = lbl, linewidth = 2)
    plot!(p_hhl, NO_range, HHL_v, label = "HHL " * lbl, linewidth = 2, linestyle = :solid)
    plot!(p_hhl, NO_range, NHL_v, label = "NHL " * lbl, linewidth = 2, linestyle = :dash)
end

plot(p_ho2, p_po3, p_hhl, layout = (3, 1), size = (600, 900), left_margin = 5 * Plots.mm)
savefig("co_fig6_4.svg") # hide
```

![Figure 6.4: CO oxidation characteristics vs NO](co_fig6_4.svg)

Panel (a) shows that [HO₂] decreases with [NO] because the HO₂ + NO reaction
consumes HO₂ more efficiently at higher NO. Panel (b) shows that P(O₃) achieves
a maximum at an intermediate [NO]; at low [NO], HO₂ is abundant but there is
insufficient NO for the O₃-producing HO₂ + NO reaction, while at high [NO],
HO₂ is depleted. Panel (c) shows the crossover from HHL-dominated (low NOx,
HO₂ self-reaction) to NHL-dominated (high NOx, OH + NO₂) HOx termination.
The maximum in P(O₃) occurs at a larger [NO] than the HHL/NHL crossover.
This reproduces Figure 6.4 of Seinfeld & Pandis.
