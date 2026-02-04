"""
    AtmosphericLifetime

ModelingToolkit.jl implementation of atmospheric constituent lifetime equations from:
Seinfeld, J.H. and Pandis, S.N. (2006). Atmospheric Chemistry and Physics,
2nd Edition, Chapter 2: Atmospheric Trace Constituents.

This module implements:
- Mass conservation equation (Eq. 2.1)
- Steady-state condition (Eq. 2.2)
- Atmospheric lifetime calculations (Eqs. 2.3-2.6)
- Multiple removal pathway lifetimes (Eqs. 2.7-2.11)
- OH reaction lifetime (Eq. 2.12)
- Tropospheric budget equations (Eqs. 2.13-2.17)

Reference:
Seinfeld, J. H. and Pandis, S. N.: Atmospheric Chemistry and Physics: From Air
Pollution to Climate Change, 2nd Edition, John Wiley & Sons, Inc., Hoboken,
New Jersey, 2006, ISBN: 978-0-471-72017-1.
"""
module AtmosphericLifetime

using ModelingToolkit
using ModelingToolkit: t, D
using DynamicQuantities
using DynamicQuantities: @register_unit

# Register the "molecule" unit as a dimensionless quantity for use in
# number density calculations (molec/cm^3) and rate constants (cm^3/molec/s)
@register_unit molec 1

export AtmosphericBudget, SpeciesLifetime, MultipleRemovalLifetime
export OHReactionLifetime, TroposphericBudget

#=============================================================================
# AtmosphericBudget: Mass Conservation (Eq. 2.1)
=============================================================================#

"""
    AtmosphericBudget(; name=:AtmosphericBudget)

Implements the fundamental mass conservation equation for an atmospheric species.

Equation 2.1:
    dQ/dt = (F_in - F_out) + (P - R)

where:
- Q: Total mass (or moles) of the species in the reservoir [mol]
- F_in: Rate of mass inflow from outside the reservoir [mol/s]
- F_out: Rate of mass outflow from the reservoir [mol/s]
- P: Rate of chemical production within the reservoir [mol/s]
- R: Rate of chemical removal within the reservoir [mol/s]

The net rate of change equals the sum of transport (inflow minus outflow) and
chemistry (production minus removal) contributions.

Equation 2.2 (Steady State):
At steady state, dQ/dt = 0, implying:
    F_in + P = F_out + R

This component also computes the net transport and net chemistry terms for
diagnostic purposes.
"""
@component function AtmosphericBudget(; name=:AtmosphericBudget)
    @variables begin
        Q(t), [description = "Total mass of species in reservoir", unit = u"mol"]
        F_in(t), [description = "Rate of mass inflow from outside reservoir", unit = u"mol/s"]
        F_out(t), [description = "Rate of mass outflow from reservoir", unit = u"mol/s"]
        P(t), [description = "Rate of chemical production in reservoir", unit = u"mol/s"]
        R(t), [description = "Rate of chemical removal in reservoir", unit = u"mol/s"]
        net_transport(t), [description = "Net transport contribution (F_in - F_out)", unit = u"mol/s"]
        net_chemistry(t), [description = "Net chemistry contribution (P - R)", unit = u"mol/s"]
    end

    eqs = [
        # Eq. 2.1 - Mass conservation equation
        D(Q) ~ (F_in - F_out) + (P - R),
        # Diagnostic: net transport
        net_transport ~ F_in - F_out,
        # Diagnostic: net chemistry
        net_chemistry ~ P - R,
    ]

    return System(eqs, t; name)
end

#=============================================================================
# SpeciesLifetime: Lifetime Calculations (Eqs. 2.3-2.6)
=============================================================================#

"""
    SpeciesLifetime(; name=:SpeciesLifetime)

Calculates atmospheric lifetime for a species using various formulations.

Implements:

Equation 2.3 (General form):
    tau = Q / (R + F_out)

This is the general definition of atmospheric lifetime, applicable to any
reservoir with both chemical removal R and outflow F_out.

Equation 2.4 (Global atmosphere, F_out = 0):
    tau = Q / R

For the global atmosphere as a whole, there is no outflow to other
reservoirs, so the lifetime is determined solely by chemical removal.

Equation 2.5 (Steady state):
    tau = Q / R = Q / P

At steady state, R = P (when F_in = F_out = 0 for the global case), so
the lifetime can be computed from either the removal or production rate.

Equation 2.6 (First-order removal):
    tau = 1 / lambda

When removal follows first-order kinetics (R = lambda * Q), the lifetime
is simply the inverse of the first-order rate constant.

Notes:
- tau_general: General lifetime accounting for both removal and outflow
- tau_removal: Lifetime based on removal only (global atmosphere)
- tau_production: Lifetime based on production (steady state)
- tau_first_order: Lifetime from first-order rate constant
"""
@component function SpeciesLifetime(; name=:SpeciesLifetime)
    @parameters begin
        lambda, [description = "First-order removal rate constant", unit = u"s^-1"]
    end

    @variables begin
        Q(t), [description = "Total mass of species in reservoir", unit = u"mol"]
        R(t), [description = "Rate of chemical removal", unit = u"mol/s"]
        P(t), [description = "Rate of chemical production", unit = u"mol/s"]
        F_out(t), [description = "Rate of mass outflow", unit = u"mol/s"]
        tau_general(t), [description = "General atmospheric lifetime (Eq. 2.3)", unit = u"s"]
        tau_removal(t), [description = "Lifetime from removal rate (Eq. 2.4)", unit = u"s"]
        tau_production(t), [description = "Lifetime from production rate (Eq. 2.5)", unit = u"s"]
        tau_first_order(t), [description = "Lifetime from first-order rate constant (Eq. 2.6)", unit = u"s"]
    end

    eqs = [
        # Eq. 2.3 - General atmospheric lifetime
        tau_general ~ Q / (R + F_out),
        # Eq. 2.4 - Lifetime from removal rate only (global atmosphere)
        tau_removal ~ Q / R,
        # Eq. 2.5 - Lifetime from production rate (steady state)
        tau_production ~ Q / P,
        # Eq. 2.6 - Lifetime from first-order rate constant
        tau_first_order ~ 1 / lambda,
    ]

    return System(eqs, t; name)
end

#=============================================================================
# MultipleRemovalLifetime: Parallel Removal Pathways (Eqs. 2.7-2.11)
=============================================================================#

"""
    MultipleRemovalLifetime(; name=:MultipleRemovalLifetime)

Calculates atmospheric lifetime when multiple removal pathways operate in parallel.

Implements:

Equation 2.7 (Two first-order processes):
    tau = 1 / (k_1 + k_2)

When two independent first-order removal processes operate simultaneously
with rate constants k_1 and k_2, the overall lifetime is determined by
their sum.

Equation 2.8 (Inverse lifetime sum):
    1/tau = 1/tau_1 + 1/tau_2

This is the general form showing that inverse lifetimes (removal rates)
are additive for parallel processes.

Equation 2.9 (Combined lifetime formula):
    tau = (tau_1 * tau_2) / (tau_1 + tau_2)

This follows directly from Eq. 2.8 and gives the combined lifetime in
terms of individual pathway lifetimes.

The component takes two removal rate constants as inputs and computes:
- tau_1, tau_2: Individual pathway lifetimes
- tau_combined: Combined lifetime from all pathways
- inverse_tau: Sum of inverse lifetimes (total removal rate)

This can be extended to N pathways where:
    1/tau = sum(1/tau_i) for i = 1 to N
"""
@component function MultipleRemovalLifetime(; name=:MultipleRemovalLifetime)
    @parameters begin
        k_1, [description = "First-order rate constant for removal pathway 1", unit = u"s^-1"]
        k_2, [description = "First-order rate constant for removal pathway 2", unit = u"s^-1"]
    end

    @variables begin
        tau_1(t), [description = "Lifetime for removal pathway 1", unit = u"s"]
        tau_2(t), [description = "Lifetime for removal pathway 2", unit = u"s"]
        tau_combined(t), [description = "Combined lifetime from all pathways (Eq. 2.9)", unit = u"s"]
        inverse_tau(t), [description = "Sum of inverse lifetimes (1/tau_1 + 1/tau_2)", unit = u"s^-1"]
    end

    eqs = [
        # Individual pathway lifetimes (Eq. 2.6 applied to each)
        tau_1 ~ 1 / k_1,
        tau_2 ~ 1 / k_2,
        # Eq. 2.7 & 2.8 - Inverse lifetime sum
        inverse_tau ~ k_1 + k_2,
        # Eq. 2.9 - Combined lifetime formula
        tau_combined ~ (tau_1 * tau_2) / (tau_1 + tau_2),
    ]

    return System(eqs, t; name)
end

#=============================================================================
# OHReactionLifetime: OH Reaction Lifetime (Eq. 2.12)
=============================================================================#

"""
    OHReactionLifetime(; name=:OHReactionLifetime)

Calculates atmospheric lifetime due to reaction with OH radicals.

Implements Equation 2.12:
    tau = 1 / (k_OH * [OH])

where:
- k_OH: Second-order rate constant for reaction with OH [cm^3/(molec*s)]
- [OH]: Concentration of OH radicals [molec/cm^3]

The OH radical is the primary oxidant in the troposphere for many trace
species. The lifetime due to OH reaction depends on both the intrinsic
reactivity of the species (k_OH) and the ambient OH concentration.

Typical values:
- k_OH: varies widely, e.g., 1e-14 to 1e-10 cm^3/(molec*s)
- [OH]: ~1e6 molec/cm^3 (global tropospheric average)
- tau: seconds to years depending on species

Note: The product k_OH * [OH] gives an effective first-order rate constant
with units of s^-1, consistent with Eq. 2.6.
"""
@component function OHReactionLifetime(; name=:OHReactionLifetime)
    @parameters begin
        k_OH, [description = "Second-order rate constant for OH reaction", unit = u"cm^3/(molec*s)"]
    end

    @variables begin
        OH_conc(t), [description = "OH radical concentration", unit = u"molec/cm^3"]
        k_eff(t), [description = "Effective first-order rate constant (k_OH * [OH])", unit = u"s^-1"]
        tau_OH(t), [description = "Lifetime due to OH reaction (Eq. 2.12)", unit = u"s"]
    end

    eqs = [
        # Effective first-order rate constant
        k_eff ~ k_OH * OH_conc,
        # Eq. 2.12 - Lifetime due to OH reaction
        tau_OH ~ 1 / k_eff,
    ]

    return System(eqs, t; name)
end

#=============================================================================
# TroposphericBudget: Complete Tropospheric Budget (Eqs. 2.13-2.17)
=============================================================================#

"""
    TroposphericBudget(; name=:TroposphericBudget)

Implements the complete tropospheric budget equation for a trace species.

Equations 2.13-2.14:
    dQ_i/dt = P_i^n + P_i^a + P_i^c - (k_i^d + k_i^w + k_i^c + k_i^t) * Q_i

This is the full mass balance for species i in the troposphere, including:

Source terms (production rates in mol/s):
- P_n: Natural emissions (biogenic, volcanic, oceanic, etc.)
- P_a: Anthropogenic emissions (combustion, industrial, agricultural, etc.)
- P_c: Chemical production within the troposphere

Removal terms (first-order rate constants in s^-1):
- k_d: Dry deposition (uptake by surfaces)
- k_w: Wet deposition (scavenging by precipitation)
- k_c: Chemical loss (reaction with OH, O3, NO3, etc.)
- k_t: Transport to stratosphere (cross-tropopause flux)

The total removal rate constant is:
    k_total = k_d + k_w + k_c + k_t

Equation 2.15 (Lifetime from removal):
    tau_i = 1 / (k_d + k_w + k_c + k_t)

Equation 2.16-2.17 (Lifetime from production at steady state):
    tau_i = Q_i / (P_n + P_a + P_c)

At steady state, the removal-based and production-based lifetimes are equal.

The component computes individual pathway lifetimes for diagnostic purposes:
- tau_dry: Lifetime against dry deposition only
- tau_wet: Lifetime against wet deposition only
- tau_chem: Lifetime against chemical loss only
- tau_transport: Lifetime against stratospheric transport only
- tau_total: Total lifetime (Eq. 2.15)
- tau_production: Lifetime from production (Eq. 2.17)
"""
@component function TroposphericBudget(; name=:TroposphericBudget)
    @parameters begin
        k_d, [description = "Dry deposition rate constant", unit = u"s^-1"]
        k_w, [description = "Wet deposition rate constant", unit = u"s^-1"]
        k_c, [description = "Chemical loss rate constant", unit = u"s^-1"]
        k_t, [description = "Transport to stratosphere rate constant", unit = u"s^-1"]
    end

    @variables begin
        Q(t), [description = "Total mass of species in troposphere", unit = u"mol"]
        P_n(t), [description = "Natural emission rate", unit = u"mol/s"]
        P_a(t), [description = "Anthropogenic emission rate", unit = u"mol/s"]
        P_c(t), [description = "Chemical production rate", unit = u"mol/s"]
        P_total(t), [description = "Total production rate", unit = u"mol/s"]
        R_total(t), [description = "Total removal rate", unit = u"mol/s"]
        k_total(t), [description = "Total first-order removal rate constant", unit = u"s^-1"]
        tau_dry(t), [description = "Lifetime against dry deposition", unit = u"s"]
        tau_wet(t), [description = "Lifetime against wet deposition", unit = u"s"]
        tau_chem(t), [description = "Lifetime against chemical loss", unit = u"s"]
        tau_transport(t), [description = "Lifetime against transport to stratosphere", unit = u"s"]
        tau_total(t), [description = "Total atmospheric lifetime (Eq. 2.15)", unit = u"s"]
        tau_production(t), [description = "Lifetime from production (Eq. 2.17)", unit = u"s"]
    end

    eqs = [
        # Total production rate
        P_total ~ P_n + P_a + P_c,
        # Total removal rate constant (Eq. 2.14)
        k_total ~ k_d + k_w + k_c + k_t,
        # Total removal rate
        R_total ~ k_total * Q,
        # Eq. 2.13 - Tropospheric budget equation
        D(Q) ~ P_total - R_total,
        # Individual pathway lifetimes
        tau_dry ~ 1 / k_d,
        tau_wet ~ 1 / k_w,
        tau_chem ~ 1 / k_c,
        tau_transport ~ 1 / k_t,
        # Eq. 2.15 - Total lifetime from removal
        tau_total ~ 1 / k_total,
        # Eq. 2.17 - Lifetime from production (valid at steady state)
        tau_production ~ Q / P_total,
    ]

    return System(eqs, t; name)
end

end # module AtmosphericLifetime
