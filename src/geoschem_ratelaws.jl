# GEOS-Chem rate laws
# Adapted from https://github.com/geoschem/geos-chem/tree/4722f288e90291ba904222f4bbe4fc216d17c34a/KPP/fullchem
# The GEOS-Chem license applies: https://github.com/geoschem/geos-chem/blob/main/LICENSE.txt

const N_A = 6.02214076e23 # Avogadro's number
const cm3_m3 = 1e6

"""
Function to create a constant rate coefficient for second-order reactions
"""
function constant_k(t, T, num_density, c; unit = u"ppb^-1*s^-1", name = :constant_k)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        c = c * N_A / cm3_m3, [unit = u"mol^-1*m^3*s^-1"]
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
    end
    @variables(k(t), [unit = unit],)

    C = num_density * ppb_unit
    ODESystem([k ~ c * C], t, [k], consts; name = name)
end

"""
Function to create a constant rate coefficient for first-order reactions
"""
function constant_k_1(t, c; unit = u"s^-1", name = :constant_k_1)
    consts = @constants begin
        c = c, [unit = u"s^-1"]
    end
    @variables(k(t), [unit = unit],)
    ODESystem([k ~ c], t, [k], consts; name = name)
end

function regress_T(t, T, num_density, a_0, b_0, T_0; name = :acet_oh)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        T_0 = T_0, [unit = u"K"]
        a_0 = a_0
        b_0 = b_0
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3*mol^-1*s^-1"]
    end
    @variables k(t) [unit = u"ppb^-1*s^-1"]

    C = num_density * ppb_unit * unit_conv
    ODESystem([k ~ (a_0 + b_0 * exp(T_0 / T)) * C], t, [k], consts; name = name)
end

"""
Arrhenius equation:

```math
    k = a0 * exp( c0 / T ) * (T/300)^b0
```
"""
function arrhenius_ppb(t, T, num_density, a0, b0, c0; unit = u"ppb^-1*s^-1", name = :arrhenius_ppb)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    t = ParentScope(t)
    consts = @constants begin
        K_300 = 300, [unit = u"K"]
        a0 = a0
        b0 = b0
        c0 = c0, [unit = u"K"]
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3*mol^-1*s^-1"]
    end
    @variables k(t) [unit = unit]
    C = num_density * ppb_unit * unit_conv
    ODESystem([k ~ a0 * exp(c0 / T) * (K_300 / T)^b0 * C], t, [k], consts; name = name)
end

"""
Compute a third-order reaction rate constant using a modified Arrhenius equation,
with output units in `ppb⁻²·s⁻¹`.

This function is appropriate when reactant concentrations are expressed in ppb (mol/mol_air × 10⁹).

The formula used is:

```math
k(t) = a₀ * exp(c₀ / T) * (300 / T)^b₀
```
"""
function arrhenius_ppb_3(
        t, T, num_density, a0, b0, c0; unit = u"ppb^-2*s^-1", name = :arrhenius_ppb_3)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    t = ParentScope(t)
    consts = @constants begin
        K_300 = 300, [unit = u"K"]
        a0 = a0
        b0 = b0
        c0 = c0, [unit = u"K"]
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        unit_conv = (N_A / cm3_m3)^2, [unit = u"m^6*mol^-2*s^-1"]
    end
    @variables k(t) [unit = unit]
    C = num_density^2 * ppb_unit^2 * unit_conv
    ODESystem([k ~ a0 * exp(c0 / T) * (K_300 / T)^b0 * C], t, [k], consts; name = name)
end

"""
Compute a second-order reaction rate constant using a modified Arrhenius equation,
with output in SI units `m³·molec⁻¹·s⁻¹`.

The formula used is:

```math
k(t) = a₀ * exp(c₀ / T) * (300 / T)^b₀
```
"""
function arrhenius_mlc_SI(
        t, T, a0, b0, c0; unit = u"m^3*molec^-1*s^-1", name = :arrhenius_mlc_SI)
    T = ParentScope(T)
    t = ParentScope(t)
    consts = @constants begin
        K_300 = 300, [unit = u"K"]
        a0 = a0
        b0 = b0
        c0 = c0, [unit = u"K"]
        unit_conv = 1e-6, [unit = u"m^3*molec^-1*s^-1"]
    end
    @variables k(t) [unit = unit]
    ODESystem(
        [k ~ a0 * (exp(c0 / T) * (K_300 / T)^b0) * unit_conv], t, [k], consts; name = name)
end

"""
Compute a first-order reaction rate constant using a modified Arrhenius equation,
with output units in `s⁻¹`.

The formula used is:

```math
k(t) = a₀ * exp(c₀ / T) * (300 / T)^b₀
```
"""
function arrhenius_mlc_1(t, T, a0, b0, c0; unit = u"s^-1", name = :arrhenius_mlc_1)
    T = ParentScope(T)
    t = ParentScope(t)
    consts = @constants begin
        K_300 = 300, [unit = u"K"]
        a0 = a0, [unit = unit]
        b0 = b0
        c0 = c0, [unit = u"K"]
    end
    @variables k(t) [unit = unit]
    ODESystem([k ~ a0 * exp(c0 / T) * (K_300 / T)^b0], t, [k], consts; name = name)
end

"""
Third body effect for pressure dependence of rate coefficients.
a1, b1, c1 are the Arrhenius parameters for the lower-limit rate.
a2, b2, c2 are the Arrhenius parameters for the upper-limit rate.
fv         is the falloff curve paramter, (see ATKINSON ET. AL (1992)
J. Phys. Chem. Ref. Data 21, P. 1145). Usually fv = 0.6.
"""
function arr_3rdbody(t, T, num_density, a1, b1, c1, a2, b2, c2, fv; unit = u"ppb^-1*s^-1",
        name = :arr_3rdbody)
    t = ParentScope(t)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        num_density_unit_inv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol", description = "multiply by num_density to obtain SI-based unitless value of num_density"]
        k_unit_inv = 1e6, [unit = u"m^-3*molec*s", description = "multiply by k to obtain cm^3*molec^-1*s^-1-based unitless value of k"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    @named alow = arrhenius_mlc_SI(t, T, a1, b1, c1)
    @named ahigh = arrhenius_mlc_SI(t, T, a2, b2, c2)
    rlow = alow.k * k_unit_inv * num_density * num_density_unit_inv
    rhigh = ahigh.k * k_unit_inv
    xyrat = rlow / rhigh
    blog = log10(xyrat)
    fexp = 1.0 / (1.0 + (blog * blog))

    @variables k(t) [unit = unit]
    C = num_density * ppb_unit * unit_conv
    ODESystem(
        [k ~ rlow * (fv^fexp) / (1.0 + xyrat) * C],
        t, [k],
        consts;
        systems = [alow, ahigh],
        name = name
    )
end

"""
Third body effect for pressure dependence of first-order rate coefficients.
a1, b1, c1 are the Arrhenius parameters for the lower-limit rate.
a2, b2, c2 are the Arrhenius parameters for the upper-limit rate.
fv         is the falloff curve paramter, (see ATKINSON ET. AL (1992)
J. Phys. Chem. Ref. Data 21, P. 1145). Usually fv = 0.6.
"""
function arr_3rdbody_1(t, T, num_density, a1, b1, c1, a2, b2, c2, fv;
        unit = u"s^-1", name = :arr_3rdbody_1)
    t = ParentScope(t)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        k_unit_inv = 1e6, [unit = u"m^-3*molec*s", description = "multiply by k to obtain cm^3*molec^-1*s^-1-based unitless value of k"]
        num_density_unit_inv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol", description = "multiply by num_density to obtain SI-based unitless value of num_density"]
        unit_conv = 1.0, [unit = u"s^-1"]
    end
    @named alow = arrhenius_mlc_SI(t, T, a1, b1, c1)
    @named ahigh = arrhenius_mlc_SI(t, T, a2, b2, c2)
    rlow = alow.k * k_unit_inv * num_density * num_density_unit_inv
    rhigh = ahigh.k * k_unit_inv
    xyrat = rlow / rhigh
    blog = log10(xyrat)
    fexp = 1.0 / (1.0 + (blog * blog))
    @variables k(t) [unit = unit]
    ODESystem(
        [k ~ rlow * (fv^fexp) / (1.0 + xyrat) * unit_conv],
        t,
        [k],
        consts;
        systems = [alow, ahigh],
        name = name
    )
end

"""
Used to compute the rate for this reactions:
HO2 + HO2 = H2O2 + O2
"""
function rate_HO2HO2(t, T, num_density, H2O, a0, c0, a1, c1; name = :rate_HO2HO2)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    H2O = ParentScope(H2O)
    @named k0 = arrhenius_mlc_SI(t, T, a0, 0.0, c0)
    @named k1 = arrhenius_mlc_SI(t, T, a1, 0.0, c1)
    consts = @constants begin
        T_0 = 2200.0, [unit = u"K"]
        k_unit_inv = 1e6, [unit = u"m^-3*molec*s", description = "multiply by k to obtain cm^3*molec^-1*s^-1-based unitless value of k"]
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        num_density_unit_inv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol", description = "multiply by num_density to obtain SI-based unitless value of num_density"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    @variables k(t) [unit = u"ppb^-1*s^-1"]
    C_H2O = num_density * num_density_unit_inv * ppb_unit
    C = num_density * ppb_unit * unit_conv
    ODESystem(
        [
            k ~ (k0.k * k_unit_inv + k1.k * k_unit_inv * num_density * num_density_unit_inv) *
                ((1.0 + 1.4E-21 * (H2O * C_H2O) * exp(T_0 / T))) * C
        ],
        t,
        [k],
        consts;
        systems = [k0, k1],
        name = name
    )
end

"""
Reaction rate for:
OH + CO = HO2 + CO2 (cf. JPL 15-10)
"""
function rate_OHCO(t, T, num_density; name = :rate_OHCO)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        k_unit_inv = 1e6, [unit = u"m^-3*molec*s", description = "multiply by k to obtain cm^3*molec^-1*s^-1-based unitless value of k"]
        num_density_unit_inv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol", description = "multiply by num_density to obtain SI-based unitless value of num_density"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    @named klo1 = arrhenius_mlc_SI(t, T, 5.9E-33, 1, 0)
    @named khi1 = arrhenius_mlc_SI(t, T, 1.1E-12, -1.3, 0)
    xyrat1 = klo1.k * num_density * num_density_unit_inv / khi1.k
    blog1 = log10(xyrat1)
    fexp1 = 1.0 / (1.0 + blog1 * blog1)
    kco1 = klo1.k * k_unit_inv * num_density * num_density_unit_inv * 0.6^fexp1 /
           (1.0 + xyrat1) 
    @named klo2 = arrhenius_mlc_SI(t, T, 1.5E-13, 0, 0)
    @named khi2 = arrhenius_mlc_SI(t, T, 2.1E+09, -6.1, 0)
    xyrat2 = klo2.k * num_density * num_density_unit_inv / khi2.k
    blog2 = log10(xyrat2)
    fexp2 = 1.0 / (1.0 + blog2 * blog2)
    kco2 = klo2.k * k_unit_inv * 0.6^fexp2 / (1.0 + xyrat2) 

    @variables k(t) [unit = u"ppb^-1*s^-1"]
    C = num_density * ppb_unit * unit_conv
    ODESystem(
        [k ~ kco1 * C + kco2 * C],
        t,
        [k],
        consts;
        systems = [klo1, khi1, klo2, khi2],
        name = name
    )
end

"""
Reaction rate for the "A" branch of these RO2 + NO reactions:
MO2  + NO = MENO3
in which the "a1" parameter equals exactly 1.

For these reactions, these Arrhenius law terms evaluate to 1:
(300/T)^b0
(300/T)^b1 * exp(c1/T)
because b0 = b1 = c1 = 0.

Special treatment for methyl nitrate based on observations
as Carter and Atkinson formulation does not apply to C1.
Value based on upper limit of Flocke et al. 1998 as applied
in Fisher et al. 2018
"""
function rate_RO2NO_a1(t, T, num_density, a0, c0; name = :rate_RO2NO_a1)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        a0 = a0
        c0 = c0, [unit = u"K"]
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    C = num_density * ppb_unit * unit_conv
    @variables k(t) [unit = u"ppb^-1*s^-1"]
    ODESystem([k ~ a0 * exp(c0 / T) * 3.0e-4 * C], t, [k], consts; name = name)
end

"""
Reaction rate for the "B" branch of these RO2 + NO reactions:
MO2 + NO = CH2O + NO2 + HO2
in which the "a1" parameter equals exactly 1.

For these reactions, these Arrhenius law terms evaluate to 1:
(300/T)^b0
(300/T)^b1 * exp(c1/T)
because b0 = c0 = c1 = 0.
"""
function rate_RO2NO_b1(t, T, num_density, a0, c0; name = :rate_RO2NO_b1)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        a0 = a0
        c0 = c0, [unit = u"K"]
        fyrno3 = 3.0e-4
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    @variables k(t) [unit = u"ppb^-1*s^-1"]
    C = num_density * ppb_unit * unit_conv
    ODESystem([k ~ a0 * exp(c0 / T) * (1 - fyrno3) * C], t, [k], consts; name = name)
end

"""
Reaction rate for the "A" branch of these RO2 + NO reactions,
ETO2 + NO = ETNO3
A3O2 + NO = NPRNO3
R4O2 + NO = R4N2
B3O2 + NO = IPRNO3
in which the "a1" parameter is greater than 1.0.
"""
function rate_RO2NO_a2(t, T, num_density, a0, c0, a1; name = :rate_RO2NO_a2)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        num_density_unit_inv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol", description = "multiply by num_density to obtain SI-based unitless value of num_density"]
        k_unit_inv = 1e6, [unit = u"m^-3*molec*s", description = "multiply by k to obtain cm^3*molec^-1*s^-1-based unitless value of k"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end

    @named k0 = arrhenius_mlc_SI(t, T, a0, 0, c0)
    @named yyyn = arrhenius_mlc_SI(t, T, 0.826, 8.1, 0.0)
    @variables(xxyn(t), 
        aaa(t),
        zzyn(t),
        rarb(t),
        fyrno3(t),
        k(t), [unit = u"ppb^-1*s^-1"])
    C = num_density * ppb_unit * unit_conv
    eqs = [xxyn ~ 1.94e-22 * exp(0.97 * a1) * num_density * num_density_unit_inv
           aaa ~ log10(xxyn / (yyyn.k * k_unit_inv))
           zzyn ~ (1.0 / (1.0 + (aaa * aaa)))
           rarb ~ (xxyn / (1.0 + (xxyn / (yyyn.k * k_unit_inv)))) *
                  (0.411^zzyn)
           fyrno3 ~ (rarb / (1.0 + rarb))
           k ~ k0.k * k_unit_inv * fyrno3 * C]
    ODESystem(
        eqs,
        t,
        [xxyn, aaa, zzyn, rarb, fyrno3, k],
        consts;
        systems = [k0, yyyn],
        name = name
    )
end

"""
Reaction rate for the "B" branch of these RO2 + NO reactions:
ETO2 + NO = NO2 +     HO2 + ...
A3O2 + NO = NO2 +     HO2 + ...
R4O2 + NO = NO2 + 0.27HO2 + ...
B3O2 + NO = NO2 +     HO2 + ...
in which the "a1" parameter is greater than 1.0.

Use this function when a1 input argument is greater than 1.0.
"""
function rate_RO2NO_b2(t, T, num_density, a0, c0, a1; name = :rate_RO2NO_b2)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        num_density_unit_inv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol", description = "multiply by num_density to obtain SI-based unitless value of num_density"]
        k_unit_inv = 1e6, [unit = u"m^-3*molec*s", description = "multiply by k to obtain cm^3*molec^-1*s^-1-based unitless value of k"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    @named k0 = arrhenius_mlc_SI(t, T, a0, 0.0, c0)
    @named yyyn = arrhenius_mlc_SI(t, T, 0.826, 8.1, 0.0)
    vars = @variables begin
        xxyn(t)
        aaa(t)
        zzyn(t)
        rarb(t)
        fyrno3(t)
        k(t), [unit = u"ppb^-1*s^-1"]
    end
    C = num_density * ppb_unit * unit_conv
    eqs = [xxyn ~ 1.94e-22 * exp(0.97 * a1) * num_density * num_density_unit_inv
           aaa ~ log10(xxyn / (yyyn.k * k_unit_inv))
           zzyn ~ (1.0 / (1.0 + (aaa * aaa)))
           rarb ~ (xxyn / (1.0 + (xxyn / (yyyn.k * k_unit_inv )))) * (0.411^zzyn)
           fyrno3 ~ (rarb / (1.0 + rarb))
           k ~ k0.k * k_unit_inv * (1.0 - fyrno3) * C]
    ODESystem(
        eqs,
        t,
        vars,
        consts;
        systems = [k0, yyyn],
        name = name
    )
end

"""
Temperature Dependent Branching Ratio
"""
function tbranch(t, T, num_density, a0, b0, c0, a1, b1, c1; name = :tbranch)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
        k_unit_inv = 1e6, [unit = u"m^-3*molec*s", description = "multiply by k to obtain cm^3*molec^-1*s^-1-based unitless value of k"]
    end
    @named k0 = arrhenius_mlc_SI(t, T, a0, b0, c0)
    @named k1 = arrhenius_mlc_SI(t, T, a1, b1, c1)
    @variables k(t) [unit = u"ppb^-1*s^-1"]
    C = num_density * ppb_unit * unit_conv
    ODESystem(
        [k ~ (k0.k * k_unit_inv / (1.0 + k1.k * k_unit_inv)) * C],
        t,
        [k],
        consts;
        systems = [k0, k1],
        name = name
    )
end

"""
Used to compute the rate for these reactions:
HNO3  + OH = H2O + NO3
HONIT + OH = NO3 + HAC
"""
function rate_OHHNO3(t, T, num_density, a0, c0, a1, c1, a2, c2; name = :rate_OHHNO3)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    # ---  OH + HNO3:   K = K0 + K3[M] / (1 + K3[M]/K2)  ------
    @named k0 = arrhenius_mlc_SI(t, T, a0, 0, c0)
    @named k1 = arrhenius_mlc_SI(t, T, a1, 0, c1)
    @named k1_5 = arrhenius_mlc_SI(t, T, a2, 0, c2)
    consts = @constants begin
        num_density_unit_inv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol", description = "multiply by num_density to obtain SI-based unitless value of num_density"]
        k_unit_inv = 1e6, [unit = u"m^-3*molec*s", description = "multiply by k to obtain cm^3*molec^-1*s^-1-based unitless value of k"]
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    k2 = (num_density * num_density_unit_inv) * (k1_5.k * k_unit_inv)
    @variables k(t) [unit = u"ppb^-1*s^-1"]
    C = num_density * ppb_unit * unit_conv
    ODESystem(
        [k ~ k0.k * k_unit_inv * C + k2 * C / (1.0 + k2 / (k1.k * k_unit_inv))],
        t,
        [k],
        consts;
        systems = [k0, k1, k1_5],
        name = name
    )
end

"""
Calculates the equilibrium constant  for second-order reactions
Find the backwards reaction by K=kforward/kbackwards
Calculates the rate constant of the forward reaction

Used to compute the rate for these reactions:
PPN        = RCO3 + NO2
PAN        = MCO3 + NO2
"""
function eq_const(t, T, num_density, a0, c0, a1, b1, a2, b2, fv; unit = u"ppb^-1*s^-1", name = :eq_const)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    @named k0 = arrhenius_ppb(t, T, num_density, a0, 0, c0)  # backwards rxn rate
    @named k1 = arr_3rdbody(t, T, num_density, a1, b1, 0, a2, b2, 0, fv)  # forwards rxn rate
    consts = @constants unit_conv = 1.0, [unit = unit]
    @variables k(t) [unit = unit]
    ODESystem(
        [k ~ k1.k / k0.k * unit_conv], t, [k], consts; systems = [k0, k1], name = name)
end
"""
Calculates the equilibrium constant for first-order reactions
Find the backwards reaction by K=kforward/kbackwards
Calculates the rate constant of the forward reaction

Used to compute the rate for these reactions:
ClOO  {+M} = Cl   + O2 {+M}
Cl2O2 {+M} = 2ClO      {+M}
"""
function eq_const_1(t, T, num_density, a0, c0, a1, b1, a2, b2, fv;
        unit = u"s^-1", name = :eq_const)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    @named k0 = arrhenius_ppb(t, T, num_density, a0, 0, c0)  # backwards rxn rate
    @named k1 = arr_3rdbody(t, T, num_density, a1, b1, 0, a2, b2, 0, fv)  # forwards rxn rate
    consts = @constants unit_conv = 1.0, [unit = unit]
    @variables k(t) [unit = unit]
    ODESystem(
        [k ~ k1.k / k0.k * unit_conv], t, [k], consts; systems = [k0, k1], name = name)
end

"""
Carbon Dependence of RO2+HO2, used in these reactions:
A3O2 + HO2 = RA3P
PO2  + HO2 = PP
KO2  + HO2 = 0.150OH + 0.150ALD2 + 0.150MCO3 + 0.850ATOOH
B3O2 + HO2 = RB3P
PRN1 + HO2 = PRPN
"""
function rate_RO2HO2(t, T, num_density, a0, c0, a1; name = :rate_RO2HO2)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    @named k0 = arrhenius_ppb(t, T, num_density, a0, 0, c0)
    @variables k(t) [unit = u"ppb^-1*s^-1"]
    ODESystem(
        [k ~ k0.k * (1.0 - exp(-0.245 * a1))],
        t,
        [k],
        [];
        systems = [k0],
        name = name
    )
end

"""
Used to compute the rate for this reaction:
GLYC + OH = 0.732CH2O + 0.361CO2  + 0.505CO    + 0.227OH

  - 0.773HO2  + 0.134GLYX + 0.134HCOOH
    which is the "A" branch of GLYC + OH.

For this reaction, these Arrhenius law terms evaluate to 1:
(300/T)^b0 * exp(c0/T)
Because b0 = c0 = 0.
"""
function rate_GLYCOH_a(t, T, num_density, a0; name = :rate_GLYCOH_a)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        a0 = a0
        exp_arg = -1.0 / 73.0, [unit = u"K^-1"]
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    glyc_frac = 1.0 - 11.0729 * exp(exp_arg * T)
    glyc_frac = max(glyc_frac, 0.0)
    @variables k(t) [unit = u"ppb^-1*s^-1"]
    C = num_density * ppb_unit * unit_conv
    ODESystem([k ~ a0 * glyc_frac * C], t, [k], consts; name = name)
end

"""
Used to compute the rate for this reaction:
GLYC + OH = HCOOH + OH + CO
which is the "B" branch of GLYC + OH.

For this reaction, these Arrhenius law terms evaluate to 1:
(300/T)^b0 * exp(c0/T)
Because b0 = c0 = 0.
"""
function rate_GLYCOH_b(t, T, num_density, a0; name = :rate_GLYCOH_b)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        a0 = a0
        exp_arg = -1.0 / 73.0, [unit = u"K^-1"]
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    glyc_frac = 1.0 - 11.0729 * exp(exp_arg * T)
    glyc_frac = max(glyc_frac, 0.0)
    @variables k(t) [unit = u"ppb^-1*s^-1"]
    C = num_density * ppb_unit * unit_conv
    ODESystem([k ~ a0 * (1.0 - glyc_frac) * C], t, [k], consts; name = name)
end

"""
Used to compute the rate for this reaction:
HAC + OH = MGLY + HO2
which is the "A" branch of HAC + OH.
"""
function rate_HACOH_a(t, T, num_density, a0, c0; name = :rate_HACOH_a)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        exp_arg = -1.0 / 60.0, [unit = u"K^-1"]
        k_unit_inv = 1e6, [unit = u"m^-3*molec*s", description = "multiply by k to obtain cm^3*molec^-1*s^-1-based unitless value of k"]
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    k0 = arrhenius_mlc_SI(t, T, a0, 0, c0)
    hac_frac = 1.0 - 23.7 * exp(exp_arg * T)
    hac_frac = max(hac_frac, 0.0)
    @variables k(t) [unit = u"ppb^-1*s^-1"]
    C = num_density * ppb_unit * unit_conv
    ODESystem(
        [k ~ k0.k * k_unit_inv * hac_frac * C], t, [k], consts; systems = [k0], name = name)
end

"""
Used to compute the rate for this reaction:
HAC + OH = 0.5HCOOH + OH + 0.5ACTA + 0.5CO2 + 0.5CO + 0.5MO2
which is the "B" branch of HAC + OH.
"""
function rate_HACOH_b(t, T, num_density, a0, c0; name = :rate_HACOH_b)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        exp_arg = -1.0 / 60.0, [unit = u"K^-1"]
        k_unit_inv = 1e6, [unit = u"m^-3*molec*s", description = "multiply by k to obtain cm^3*molec^-1*s^-1-based unitless value of k"]
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    k0 = arrhenius_mlc_SI(t, T, a0, 0, c0)
    hac_frac = 1.0 - 23.7 * exp(exp_arg * T)
    hac_frac = max(hac_frac, 0.0)
    @variables k(t) [unit = u"ppb^-1*s^-1"]
    C = num_density * ppb_unit * unit_conv
    ODESystem([k ~ k0.k * k_unit_inv * (1.0 - hac_frac) * C],
        t, [k], consts; systems = [k0], name = name)
end

"""
Reaction rate for:
DMS + OH = 0.750SO2 + 0.250MSA + MO2
"""
function rate_DMSOH(t, T, num_density, a0, c0, a1, c1; name = :rate_DMSOH)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        c2 = 0.2095e0
        num_density_unit_inv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol", description = "multiply by num_density to obtain SI-based unitless value of num_density"]
        k_unit_inv = 1e6, [unit = u"m^-3*molec*s", description = "multiply by k to obtain cm^3*molec^-1*s^-1-based unitless value of k"]
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    @named k0 = arrhenius_mlc_SI(t, T, a0, 0, c0)
    @named k1 = arrhenius_mlc_SI(t, T, a1, 0, c1)
    @variables k(t) [unit = u"ppb^-1*s^-1"]
    C = num_density * ppb_unit * unit_conv
    ODESystem(
        #[k ~ unit_conv * (k0.k * num_density * c2) / (1 / one_s + k1.k * c2)],
        [k ~ (k0.k * k_unit_inv * num_density * num_density_unit_inv * c2) /
             (1.0 + k1.k * k_unit_inv * c2) * C],
        t,
        [k],
        consts;
        systems = [k0, k1],
        name = name
    )
end

"""
Reaction rate for:
GLYX + NO3 = HNO3 + HO2 + 2CO
i.e. the HO2 + 2*CO branch
"""
function rate_GLYXNO3(t, T, num_density, a0, c0; name = :rate_GLYXNO3)
    # ---  K = K1*([O2]+3.5D18)/(2*[O2]+3.5D18)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    @named k0 = arrhenius_mlc_SI(t, T, a0, 0, c0)
    consts = @constants begin
        num_density_unit_inv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol", description = "multiply by num_density to obtain SI-based unitless value of num_density"]
        k_unit_inv = 1e6, [unit = u"m^-3*molec*s", description = "multiply by k to obtain cm^3*molec^-1*s^-1-based unitless value of k"]
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    @variables k(t) [unit = u"ppb^-1*s^-1"]
    O2 = num_density * num_density_unit_inv * 0.2095
    C = num_density * ppb_unit * unit_conv
    ODESystem(
        [k ~ k0.k * k_unit_inv * (O2 + 3.5E+18) / (2.0 * O2 + 3.5E+18) * C],
        t,
        [k],
        consts;
        systems = [k0],
        name = name
    )
end

"""
Modified Arrhenius law with output in units `s⁻¹`.
"""
function arrplus_mlc_1(t, T, a0, b0, c0, d0, e0; unit = u"s^-1", name = :arrplus_mlc_1)
    T = ParentScope(T)
    consts = @constants begin
        K_300 = 300, [unit = u"K"]
        a0 = a0, [unit = unit]
        b0 = b0, [unit = u"K"]
        e0 = e0, [unit = u"K^-1"]
        zero = 0.0, [unit = unit]
    end
    vars = @variables begin
        k(t), [unit = unit]
        kx(t), [unit = unit]
    end
    ODESystem(
        [kx ~ a0 * (d0 + (T * e0)) * exp(-b0 / T) * (T / K_300)^c0
         k ~ max(kx, zero)],
        t,
        vars,
        consts;
        name = name
    )
end

"""
Modified Arrhenius law with output in units `ppb·s⁻¹`.
"""
function arrplus_ppb(
        t, T, num_density, a0, b0, c0, d0, e0; unit = u"ppb^-1*s^-1", name = :arrplus_ppb)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        K_300 = 300, [unit = u"K"]
        a0 = a0
        b0 = b0, [unit = u"K"]
        e0 = e0, [unit = u"K^-1"]
        zero = 0.0, [unit = unit]
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    @variables(k(t), [unit = unit], kx(t), [unit = unit],)
    C = num_density * ppb_unit * unit_conv
    ODESystem(
        [kx ~ a0 * (d0 + (T * e0)) * exp(-b0 / T) * (T / K_300)^c0 * C
         k ~ max(kx, zero)],
        t,
        [k, kx],
        consts;
        name = name
    )
end


"""
Computes a temperature-dependent reaction rate constant using a modified
Arrhenius expression with additional tunneling effect terms.
Used to compute the rate for these reactions with output in units `s⁻¹`:
IHOO1 = 1.5OH + ...
IHOO4 = 1.5OH + ...
"""
function tunplus_mlc_1(t, T, a0, b0, c0, d0, e0; unit = u"s^-1", name = :tunplus_mlc_1)
    T = ParentScope(T)
    consts = @constants begin
        a0 = a0, [unit = unit]
        b0 = b0, [unit = u"K"]
        c0 = c0, [unit = u"K^3"]
        e0 = e0, [unit = u"K^-1"]
        zero = 0.0, [unit = unit]
    end
    @variables(k0(t), [unit = unit],
        k(t), [unit = unit],)
    ODESystem(
        [k0 ~ a0 * (d0 + (T * e0)) * exp(b0 / T) * exp(c0 / T^3)
         k ~ max(k0, zero)],
        t,
        [k0, k],
        consts;
        name = name
    )
end

"""
Computes a temperature-dependent reaction rate constant using a modified
Arrhenius expression with additional tunneling effect terms.
Used to compute the rate for these reactions with output in units `ppb·s⁻¹`:
IHOO1 = 1.5OH + ...
IHOO4 = 1.5OH + ...
"""
function tunplus_ppb(
        t, T, num_density, a0, b0, c0, d0, e0; unit = u"ppb^-1*s^-1", name = :tunplus_ppb)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        a0 = a0
        b0 = b0, [unit = u"K"]
        c0 = c0, [unit = u"K^3"]
        e0 = e0, [unit = u"K^-1"]
        zero = 0.0, [unit = unit]
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    @variables(k0(t), [unit = unit], k(t), [unit = unit],)
    C = num_density * ppb_unit * unit_conv
    ODESystem(
        [k0 ~ a0 * (d0 + (T * e0)) * exp(b0 / T) * exp(c0 / T^3) * C
         k ~ max(k0, zero)],
        t,
        [k0, k],
        consts;
        name = name
    )
end

"""
Used to compute the rate for these reactions:
ISOP + OH = LISOPOH + IHOO1
ISOP + OH = LISOPOH + IHOO4
"""
function rate_ISO1(t, T, num_density, a0, b0, c0, d0, e0, f0, g0; name = :rate_ISO1)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        ct = 1.0E8, [unit = u"K^3"]
        a0 = a0
        b0 = b0, [unit = u"K"]
        e0 = e0, [unit = u"K"]
        g0 = g0, [unit = u"K"]
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    k0 = d0 * exp(e0 / T) * exp(ct / T^3)
    k1 = f0 * exp(g0 / T)
    k2 = c0 * k0 / (k0 + k1)
    @variables k(t) [unit = u"ppb^-1*s^-1"]
    C = num_density * ppb_unit * unit_conv
    ODESystem([k ~ a0 * exp(b0 / T) * (1.0 - k2) * C], t, [k], consts; name = name)
end

"""
Used to compute the rate for these reactions:
ISOP + OH = 0.3MCO3 + 0.3MGLY + 0.3CH2O

  - 0.15HPALD3 + 0.25HPALD1 + 0.4HO2
  - 0.6CO + 1.5OH + 0.3HPETHNL + LISOPOH
    ISOP + OH = 0.3CH2O + 0.15HPALD4 + 0.25HPALD2
  - 1.5OH + 0.9CO + 0.7HO2 + 0.3MGLY
  - 0.3ATOOH + LISOPOH
"""
function rate_ISO2(t, T, num_density, a0, b0, c0, d0, e0, f0, g0; name = :rate_ISO2)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        ct = 1.0E8, [unit = u"K^3"]
        a0 = a0
        b0 = b0, [unit = u"K"]
        e0 = e0, [unit = u"K"]
        g0 = g0, [unit = u"K"]
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    k0 = d0 * exp(e0 / T) * exp(ct / T^3)
    k1 = f0 * exp(g0 / T)
    k2 = c0 * k0 / (k0 + k1)
    @variables k(t) [unit = u"ppb^-1*s^-1"]
    C = num_density * ppb_unit * unit_conv
    ODESystem([k ~ a0 * exp(b0 / T) * k2 * C], t, [k], consts; name = name)
end

"""
Used to compute the rate for these reactions:
RIPA   + OH = 0.67IEPOXA   + 0.33IEPOXB   + OH + 0.005LVOC
RIPB   + OH = 0.68IEPOXA   + 0.321IEPOB   + OH + 0.005LVOC
IEPOXA + OH = 0.67IEPOXA00 + 0.33IEPOXB00
IEPOXB + OH = 0.81IEPOXA00 + 0.19IEPOXB00
IHN2   + OH = 0.67IEPOXA   + 0.33IEPOXB   + NO2
IHN3   + OH = 0.67IEPOXA   + 0.33IEPOXB   + NO2
IHN1   + OH = IEPOXD       + NO2
IHN4   + OH = IEPOXD       + NO2
INPB   + OH = OH           + ITHN
INPD   + OH = OH           + ITHN
INPD   + OH = NO2          + ICHE
ICN    + OH = NO2          + ICHE
"""
function rate_EPO(t, T, num_density, a1, e1, m1; name = :rate_EPO)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        e1 = e1, [unit = u"K"]
        a1 = a1
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        num_density_unit_inv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol", description = "multiply by num_density to obtain SI-based unitless value of num_density"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    k1 = 1.0 / (m1 * num_density * num_density_unit_inv + 1.0)
    @variables k(t) [unit = u"ppb^-1*s^-1"]
    C = num_density * ppb_unit * unit_conv
    ODESystem([k ~ a1 * exp(e1 / T) * k1 * C], t, [k], consts; name = name)
end

"""
Used to compute the rate for reaction:
BZPAN --> BZCO3 + NO2
"""
function rate_PAN_abab(t, T, num_density, a0, b0, a1, b1, cf; unit = u"s^-1", name = :rate_PAN_abab)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        a0 = a0
        a1 = a1
        b0 = b0, [unit = u"K"]
        b1 = b1, [unit = u"K"]
        num_density_unit_inv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol", description = "multiply by num_density to obtain SI-based unitless value of num_density"]
        unit_conv = 1.0 , [unit = u"s^-1"]
    end
    k0 = a0 * exp(b0 / T)
    k1 = a1 * exp(b1 / T)
    k0 = k0 * num_density * num_density_unit_inv
    kr = k0 / k1
    nc = 0.75 - 1.27 * (log10(cf))
    f = 10.0^(log10(cf) / (1.0 + (log10(kr) / nc)^2))
    @variables k(t) [unit = unit]
    ODESystem([k ~ k0 * k1 * f / (k0 + k1) * unit_conv], t, [k], consts; name = name)
end

"""
Used to compute the rate for reaction:
MACR1OO + NO2 --> MPAN
MACRNO2 + NO2 --> MPAN + NO2
BZCO3 + NO2 --> BZPAN
"""
function rate_PAN_acac(t, T, num_density, a0, c0, a1, c1, cf; name = :rate_PAN_acac)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        K_300 = 300, [unit = u"K"]
        a0 = a0
        a1 = a1
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        num_density_unit_inv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol", description = "multiply by num_density to obtain SI-based unitless value of num_density"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    k0 = a0 * (T / K_300)^c0
    k1 = a1 * (T / K_300)^c1
    k0 = k0 * num_density * num_density_unit_inv
    kr = k0 / k1 # no unit
    nc = 0.75 - 1.27 * (log10(cf))# no unit
    f = 10.0^(log10(cf) / (1.0 + (log10(kr) / nc)^2)) # no unit
    @variables k(t) [unit = u"ppb^-1*s^-1"]
    C = num_density * ppb_unit * unit_conv
    ODESystem([k ~ k0 * k1 * f / (k0 + k1) * C], t, [k], consts; name = name)
end

"""
Used to compute the rate for these reactions:
IHOO1    + NO = IHN2
IHOO4    + NO = IHN4
IHPOO1   + NO = IHTN
IHPOO2   + NO = IHTN
IHPOO2   + NO = IHTN
IEPOXAOO + NO = IHTN
IEPOXBOO + NO = IHTN
IHCOO    + NO = IHTN
ISOPNOO1 + NO = IDN
ISOPNOO2 + NO = IDN
IDHNDOO1 + NO = IDN
IDHNDOO2 + NO = IDN
INO2B    + NO = IDN
INO2D    + NO = IDN
IHPNBOO  + NO = IDN
IHPNDOO  + NO = IDN
MVK0HOO  + NO = 0.438MVKN
MCROHOO  + NO = MCRHN
"""
function rate_NIT(t, T, num_density, a0, b0, c0, n, x0, y0; name = :rate_NIT)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        T_298 = 298.0, [unit = u"K"]
        a0 = a0 #[unit=u"ppb^-1*s^-1"]
        b0 = b0, [unit = u"K"]
        y0 = y0, [unit = u"K^-1"]
        zero = 0.0, [unit = u"ppb^-1*s^-1"]
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        num_density_unit_inv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol", description = "multiply by num_density to obtain SI-based unitless value of num_density"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    vars = @variables begin
        k0(t)
        k1(t)
        k2_(t)
        k2(t)
        k3(t)
        k4(t)
        kx(t), [unit = u"ppb^-1*s^-1"]
        k(t), [unit = u"ppb^-1*s^-1"]
    end
    C = num_density * ppb_unit * unit_conv
    eqs = [k0 ~ 2.0E-22 * exp(n) * (num_density * num_density_unit_inv)
           k1 ~ k0 / (4.3E-1 * (T / T_298)^(-8))
           k2_ ~ (k0 / (1.0 + k1)) * 4.1E-1^(1.0 / (1.0 + (log10(k1))^2))
           k3 ~ k2_ / (k2_ + c0)
           k4 ~ a0 * (x0 - T * y0)
           kx ~ k4 * exp(b0 / T) * k3 * C
           k ~ max(kx, zero)]
    ODESystem(eqs, t, vars, consts; name = name)
end

"""
Used to compute the rate for these reactions:
IHOO1    + NO =      NO2 + ...
IHOO4    + NO =      NO2 + ...
IHP001   + NO =      NO2 + ...
IHP002   + NO =      NO2 + ...
IHP003   + NO =      NO2 + ...
IEPOXAOO + NO =      NO2 + ...
IEPOXBOO + NO =      NO2 + ...
ICHOO    + NO =      NO2 + ...
ISOPNOO1 + NO = 1.728NO2 + ...
ISOPNOO2 + NO =      NO2 + ...
IDHNDOO1 + NO =      NO2 + ...
IDHNDOO2 + NO =      NO2 + ...
IDHNBOO  + NO =      NO2 + ...
IDHNDOO  + NO =      NO2 + ...
INO2B    + NO = 2.000NO2 + ...
INO2D    + NO =      NO2 + ...
IHPNBOO  + NO = 1.065NO2 + ...
IHPNDOO  + NO =      NO2 + ...
MVKOHOO  + NO =      NO2 + ...
MCROHOO  + NO =      NO2 + ...
"""
function rate_ALK(t, T, num_density, a0, b0, c0, n, x0, y0; name = :rate_ALK)
    T = ParentScope(T)
    num_density = ParentScope(num_density)
    consts = @constants begin
        T_298 = 298.0, [unit = u"K"]
        a0 = a0 #[unit=u"ppb^-1*s^-1"]
        b0 = b0, [unit = u"K"]
        y0 = y0, [unit = u"K^-1"]
        zero = 0, [unit = u"ppb^-1*s^-1"]
        ppb_unit = 1e-9, [unit = u"ppb^-1", description = "Convert from mol/mol_air to ppb"]
        num_density_unit_inv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol", description = "multiply by num_density to obtain SI-based unitless value of num_density"]
        unit_conv = 1.0 / cm3_m3 * N_A, [unit = u"m^3/mol/s"]
    end
    vars = @variables begin
        k0(t)
        k1(t)
        k2(t)
        k3(t)
        k4(t)
        kx(t), [unit = u"ppb^-1*s^-1"]
        k(t), [unit = u"ppb^-1*s^-1"]
    end
    C = num_density * ppb_unit * unit_conv
    eqs = [k0 ~ 2.0E-22 * exp(n) * num_density * num_density_unit_inv
           k1 ~ k0 / 4.3E-1 * (T / T_298)^(-8)
           k2 ~ (k0 / (1.0 + k1)) * 4.1E-1^(1.0 / (1.0 + (log10(k1))^2))
           k3 ~ c0 / (k2 + c0)
           k4 ~ a0 * (x0 - T * y0)
           kx ~ k4 * exp(b0 / T) * k3 * C
           k ~ max(kx, zero)]
    ODESystem(eqs, t, vars, consts; name = name)
end
