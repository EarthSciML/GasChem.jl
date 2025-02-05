export SuperFast

""" 
Convert reaction rate value in unit of cm^3/molec/s to (s*ppb)^-1
"""
function rate_toppb(t, T, P, a0; name=:rateconvert)
    T = ParentScope(T)
    P = ParentScope(P)
    @constants(
        A = 6.02e23, [unit = u"molec/mol", description = "Avogadro's number"],
        R = 8.314e6, [unit = u"(Pa*cm^3)/(K*mol)", description = "universal gas constant"],
        ppb_unit = 1e-9, [unit = u"ppb", description = "Convert from mol/mol_air to ppb"],
        num_density_unit_inv = 1, [unit = u"cm^3/molec", description = "multiply by num_density to obtain the unitless value of num_density"],
        a0 = a0, [unit = u"cm^3/molec/s"],
    )
    air_volume = R*T/P #unit in cm^3/mol
    c = A / air_volume * ppb_unit #Convert the second reaction rate value, which corresponds to species with units of molec/cm³, to ppb.
    @variables k(t) [unit = u"(s*ppb)^-1"]
    ODESystem([k ~ a0 * c], t, [k], []; name=name)
end

""" 
Arrhenius equation:
``` math
    k = a0 * exp( c0 / T ) * (T/300)^b0
```
"""
function arrh(t, T, P, a0, b0, c0; name=:arrhenius)
    T = ParentScope(T)
    t = ParentScope(t)
    P = ParentScope(P)
    @constants(
        A = 6.02e23, [unit = u"molec/mol", description = "Avogadro's number"],
        R = 8.314e6, [unit = u"(Pa*cm^3)/(K*mol)", description = "universal gas constant"],
        ppb_unit = 1e-9, [unit = u"ppb", description = "Convert from mol/mol_air to ppb"],
        num_density_unit_inv = 1, [unit = u"cm^3/molec", description = "multiply by num_density to obtain the unitless value of num_density"],
        K_300 = 300, [unit = u"K"],
        a0 = a0, [unit = u"cm^3/molec/s"],
        b0 = b0,
        c0 = c0, [unit = u"K"],
    )
    @variables k(t) [unit = u"(s*ppb)^-1"]
    air_volume = R*T/P #unit in cm^3/mol
    c = A / air_volume * ppb_unit #Convert the second reaction rate value, which corresponds to species with units of molec/cm³, to ppb.
    ODESystem([k ~ a0 * exp(c0 / T) * (K_300 / T)^b0 * c], t, [k], []; name=name)
end

"""
Third body effect for pressure dependence of rate coefficients.
a1, b1, c1 are the Arrhenius parameters for the lower-limit rate.
a2, b2, c2 are the Arrhenius parameters for the upper-limit rate.
fv         is the falloff curve paramter, (see ATKINSON ET. AL (1992)
           J. Phys. Chem. Ref. Data 21, P. 1145). Usually fv = 0.6.
"""
function arr_3rd(t, T, P, a1, b1, c1, a2, b2, c2, fv; name=:arr_3rdbody)
    t = ParentScope(t)
    T = ParentScope(T)
    P = ParentScope(P)
    @constants(
        A = 6.02e23, [unit = u"molec/mol", description = "Avogadro's number"],
        R = 8.314e6, [unit = u"(Pa*cm^3)/(K*mol)", description = "universal gas constant"],
        ppb_unit = 1e-9, [unit = u"ppb", description = "Convert from mol/mol_air to ppb"],
        num_density_unit_inv = 1, [unit = u"cm^3/molec", description = "multiply by num_density to obtain the unitless value of num_density"],
    )
    num_density_unitless = A*P/(R*T)*num_density_unit_inv

    @named alow = arrh(t, T, P, a1, b1, c1)
    @named ahigh = arrh(t, T, P, a2, b2, c2)
    rlow = alow.k * num_density_unitless
    rhigh = ahigh.k
    xyrat = rlow / rhigh
    blog = log10(xyrat)
    fexp = 1.0 / (1.0 + (blog * blog))
    @variables k(t) [unit = u"(s*ppb)^-1"]
    ODESystem([k ~ rlow * (fv^fexp) / (1.0 + xyrat)], t, [k], []; systems=[alow, ahigh], name=name)
end

"""
 Used to compute the rate for this reactions:
    HO2 + HO2 = H2O2 + O2
"""
function rate_2HO2(t, T, P, H2O, a0, c0, a1, c1; name=:rate_HO2HO2)
    T = ParentScope(T)
    P = ParentScope(P)
    H2O = ParentScope(H2O)
    @constants(
        T_0 = 2200.0, [unit = u"K"],
        A = 6.02e23, [unit = u"molec/mol", description = "Avogadro's number"],
        R = 8.314e6, [unit = u"(Pa*cm^3)/(K*mol)", description = "universal gas constant"],
        ppb_unit = 1e-9, [unit = u"ppb", description = "Convert from mol/mol_air to ppb"],
        num_density_unit_inv = 1, [unit = u"cm^3/molec", description = "multiply by num_density to obtain the unitless value of num_density"],
        ppb_inv = 1, [unit = u"ppb^-1"],
    )
    num_density_unitless = A*P/(R*T)*num_density_unit_inv
    H2O_ppb_molec_cm3 = H2O*ppb_inv*1e-9*num_density_unitless # convert value of H2O concentration in unit of ppb to unit of molec/cm3, but here is unitless
    @named k0 = arrh(t, T, P, a0, 0.0, c0)
    @named k1 = arrh(t, T, P, a1, 0.0, c1)

    @variables k(t) [unit = u"(s*ppb)^-1"]
    ODESystem([k ~ (k0.k + k1.k * num_density_unitless) * (1.0 + 1.4e-21 * H2O * H2O_ppb_molec_cm3 * exp(T_0 / T))], t, [k], [];
        systems=[k0, k1], name=name)
end

"""
Reaction rate for:
   OH + CO = HO2 + CO2 (cf. JPL 15-10)
"""
function rate_OH_CO(t, T, P; name=:rate_OHCO)
    T = ParentScope(T)
    P = ParentScope(P)
    @constants(
        A = 6.02e23, [unit = u"molec/mol", description = "Avogadro's number"],
        R = 8.314e6, [unit = u"(Pa*cm^3)/(K*mol)", description = "universal gas constant"],
        ppb_unit = 1e-9, [unit = u"ppb", description = "Convert from mol/mol_air to ppb"],
        num_density_unit_inv = 1, [unit = u"cm^3/molec", description = "multiply by num_density to obtain the unitless value of num_density"],
        ppb_inv = 1, [unit = u"ppb^-1"],
    )
    num_density_unitless = A*P/(R*T)*num_density_unit_inv
    @named klo1 = arrh(t, T, P, 5.9e-33, 1, 0)
    @named khi1 = arrh(t, T, P, 1.1e-12, -1.3, 0)
    xyrat1 = klo1.k * num_density_unitless / khi1.k
    blog1 = log10(xyrat1)
    fexp1 = 1.0 / (1.0 + blog1 * blog1)
    kco1 = klo1.k * num_density_unitless * 0.6^fexp1 / (1.0 + xyrat1)
    @named klo2 = arrh(t, T, P, 1.5e-13, 0, 0)
    @named khi2 = arrh(t, T, P, 2.1e9, -6.1, 0)
    xyrat2 = klo2.k * num_density_unitless / khi2.k
    blog2 = log10(xyrat2)
    fexp2 = 1.0 / (1.0 + blog2 * blog2)
    kco2 = klo2.k * 0.6^fexp2 / (1.0 + xyrat2)
    @variables k(t) [unit = u"(s*ppb)^-1"]
    ODESystem([k ~ kco1 + kco2], t, [k], []; systems=[klo1, khi1, klo2, khi2], name=name)
end

struct SuperFastCoupler
    sys
end

"""
    SuperFast()

This atmospheric chemical system model is built based on the Super Fast Chemical Mechanism, which is one of the simplest representations of atmospheric chemistry. It can efficiently simulate background tropheric ozone chemistry and perform well for those species included in the mechanism. The chemical equations we used is included in the supporting table S2 of the paper:

"Evaluating simplified chemical mechanisms within present-day simulations of the Community Earth System Model version 1.2 with CAM4 (CESM1.2 CAM-chem):
MOZART-4 vs. Reduced Hydrocarbon vs. Super-Fast chemistry" (2018), Benjamin Brown-Steiner, Noelle E. Selin, Ronald G. Prinn, Simone Tilmes, Louisa Emmons, Jean-François Lamarque, and Philip Cameron-Smith.

The input of the function is Temperature, concentrations of all chemicals, and reaction rates of photolysis reactions 

If the keyword argument `rxn_sys` is set to `true`, the function will return a reaction system instead of an ODE system.

# Example
```
using GasChem, EarthSciMLBase, DifferentialEquations, Plots
rs = SuperFast()
sol = solve(ODEProblem(structural_simplify(rs), [], (0,360), [], combinatoric_ratelaws=false), AutoTsit5(Rosenbrock23()), saveat=10.0)
plot(sol)
```
"""
function SuperFast(;name=:SuperFast, rxn_sys=false)
    # Create reaction rate constant system constructors
    rate_systems = []
    i = 1
    function cons(c, T, P)
        sys = rate_toppb(t, T, P, c, name=Symbol(:rateconvert, i))
        i+=1
        push!(rate_systems, sys)
        return sys.k
    end
    function arr(T, P, a0, b0, c0)
        sys = arrh(t, T, P, a0, b0, c0, name=Symbol(:arrhenius, i))
        i+=1
        push!(rate_systems, sys)
        return sys.k
    end
    function arr3(T, P, a1, b1, c1, a2, b2, c2, fv)
        sys = arr_3rd(t, T, P, a1, b1, c1, a2, b2, c2, fv; name=Symbol(:arr_3rdbody, i))
        i+=1
        push!(rate_systems, sys)
        return sys.k
    end
    function rHO2HO2(T, P, H2O, a0, c0, a1, c1)
        sys = rate_2HO2(t, T, P, H2O, a0, c0, a1, c1, name=Symbol(:rate_HO2HO2, i))
        i+=1
        push!(rate_systems, sys)
        return sys.k
    end
    function rOHCO(T, P)
        sys = rate_OH_CO(t, T, P; name=Symbol(:rate_OHCO, i))
        i+=1
        push!(rate_systems, sys)
        return sys.k
    end

    rx_sys = @reaction_network SuperFast begin
        @ivs t [unit = u"s"]
        @parameters(
            jO32OH = 2.27e-4, [unit = u"s^-1"],
            jH2O2 = 1.0097e-5, [unit = u"s^-1"],
            jNO2 = 0.0149, [unit = u"s^-1"],
            jCH2Oa = 0.00014, [unit = u"s^-1"],
            jCH2Ob = 0.00014, [unit = u"s^-1"],
            jCH3OOH = 8.9573e-6, [unit = u"s^-1"],
            T = 280.0, [unit = u"K", description = "Temperature"],
            P = 101325, [unit = u"Pa", description = "Pressure"],
            O2 = 2.1e8, [isconstantspecies=true,unit = u"ppb"],
            CH4 = 1700.0, [isconstantspecies=true, unit = u"ppb"],
            H2O = 450.0, [isconstantspecies=true, unit = u"ppb"],
        )

        @species(
            O3(t) = 20.0, [unit = u"ppb"],
            OH(t) = 0.01, [unit = u"ppb"],
            HO2(t) = 0.01, [unit = u"ppb"],
            NO(t) = 10.0, [unit = u"ppb"],
            NO2(t) = 10.0, [unit = u"ppb"],
            CH3O2(t) = 0.01, [unit = u"ppb"],
            CH2O(t) = 0.15, [unit = u"ppb"],
            CO(t) = 275.0, [unit = u"ppb"],
            CH3OOH(t) = 1.6, [unit = u"ppb"],
            #DMS(t) = 50.0, [unit = u"ppb"],
            #SO2(t) = 2.0, [unit = u"ppb"],
            ISOP(t) = 0.15, [unit = u"ppb"],
            H2O2(t) = 2.34, [unit = u"ppb"],
            HNO3(t) = 10.0, [unit = u"ppb"],
        )

        #Gas-phase reactions
        arr(T, P, 1.7e-12, 0.0, -940.0), O3 + OH --> HO2 + O2
        arr(T, P, 1.00e-14, 0.0, -490.0), O3 + HO2 --> OH + O2 + O2
        arr(T, P, 4.80e-11, 0.0, 250.0), OH + HO2 --> H2O + O2
        rHO2HO2(T, P, H2O, 3.00e-13, 460.0, 2.1e-33, 920.0), HO2 + HO2 --> H2O2 + O2
        arr(T, P, 1.80e-12, 0.0, 0.0), OH + H2O2 --> H2O + HO2
        arr(T, P, 3.00e-12, 0.0, -1500.0), O3 + NO --> NO2 + O2
        arr(T, P, 3.30e-12, 0.0, 270.0), HO2 + NO --> OH + NO2
        arr3(T, P, 1.80e-30, 3.0, 0.0, 2.8e-11, 0.0, 0.0, 0.6), NO2 + OH --> HNO3
        arr(T, P, 2.45e-12, 0.0, -1775.0), OH + CH4 --> CH3O2 + H2O
        rOHCO(T, P), OH + CO --> HO2
        arr(T, P, 5.50e-12, 0.0, 125.0), CH2O + OH --> CO + HO2 + H2O 
        arr(T, P, 4.10e-13, 0.0, 750.0), CH3O2+ HO2 --> CH3OOH + O2                       
        arr(T, P, 2.66e-12, 0.0, 200.0), CH3OOH + OH --> CH3O2+ H2O                        
        arr(T, P, 1.14e-12, 0.0, 200.0), CH3OOH + OH --> CH2O + OH + H2O                   
        arr(T, P, 2.80e-12, 0.0, 300.0), CH3O2+ NO --> CH2O + HO2 + NO2                
        arr(T, P, 9.50e-14, 0.0, 390.0), CH3O2+ CH3O2 --> 2.000CH2O + 0.800HO2  #{two different reactions simplified into one}
        cons(4.00e-24, T, P), H2O + NO2 --> 0.5HNO3                       
        arr(T, P, 2.7e-11, 0.0, 390.0), ISOP + OH --> CH3O2+ CH3O2 #{Isoprene chemistry parametrized from UCI for ISOP + OH}
        arr(T, P, 2.7e-11, 0.0, 390.0), ISOP + OH --> ISOP   #{Isoprene chemistry parametrized from UCI for ISOP + OH}
        arr(T, P, 2.7e-11, 0.0, 390.0), ISOP + OH --> ISOP + 0.500OH #{Isoprene chemistry parametrized from UCI for ISOP + OH}
        arr(T, P, 5.59e-15, 0.0, -1814.0), ISOP + O3 --> 0.870CH2O + 1.860CH3O2+ 0.060HO2 + 0.050CO   #{Isoprene chemistry parametrized from LLNL-IMPACT for ISOP + O3}

        #photolysis reactions
        jO32OH, O3 --> 2OH #simplified reaction of: O3 --> O2 + O1d; O1d + H2O --> 2OH
        jH2O2, H2O2 --> 2OH
        jNO2, NO2 --> NO + O3
        jCH2Oa, CH2O --> HO2 + HO2 + CO
        jCH2Ob, CH2O --> CO
        jCH3OOH, CH3OOH --> CH2O + HO2 + OH
    end
    rxns = compose(rx_sys, rate_systems)
    if rxn_sys
        return rxns
    end 
    # We set `combinatoric_ratelaws=false` because we are modeling macroscopic rather than microscopic behavior. 
    # See [here](https://docs.juliahub.com/ModelingToolkit/Qmdqu/3.14.0/systems/ReactionSystem/#ModelingToolkit.oderatelaw) 
    # and [here](https://github.com/SciML/Catalyst.jl/issues/311).
    convert(ODESystem, complete(rxns); combinatoric_ratelaws=false, name=name, metadata=Dict(:coupletype => SuperFastCoupler))
end