export superfast

# Add unit "ppb" to Unitful 
module MyUnits
using Unitful
@unit ppb "ppb" Number 1 / 1000000000 false
end
Unitful.register(MyUnits)

"""
    superfast()

This atmospheric chemical system model is built based on the Super Fast Chemical Mechanism, which is one of the simplest representations of atmospheric chemistry. It can efficiently simulate background tropheric ozone chemistry and perform well for those species included in the mechanism. The chemical equations we used is included in the supporting table S2 of the paper:

"Evaluating simplified chemical mechanisms within present-day simulations of the Community Earth System Model version 1.2 with CAM4 (CESM1.2 CAM-chem):
MOZART-4 vs. Reduced Hydrocarbon vs. Super-Fast chemistry" (2018), Benjamin Brown-Steiner, Noelle E. Selin, Ronald G. Prinn, Simone Tilmes, Louisa Emmons, Jean-François Lamarque, and Philip Cameron-Smith.

The input of the function is Temperature, concentrations of all chemicals, and reaction rates of photolysis reactions 

# Example
```
using OrdinaryDiffEq, Plots
rs = superfast()
sol = solve(ODEProblem(rs, [], (0,360), [], combinatoric_ratelaws=false), Tsit5())
plot(sol)
```

We set `combinatoric_ratelaws=false` because we are modeling macroscopic rather than microscopic behavior. See [here](https://docs.juliahub.com/ModelingToolkit/Qmdqu/3.14.0/systems/ReactionSystem/#ModelingToolkit.oderatelaw) and [here](https://github.com/SciML/Catalyst.jl/issues/311).
"""
function superfast(t)
    @variables jO31D(t) = 4.0 * 10.0^-3 [unit = u"s^-1"]
    @parameters j2OH = 2.2 * 10.0^-10 [unit = u"(s*ppb)^-1"]
    @variables jH2O2(t) = 1.0097 * 10.0^-5 [unit = u"s^-1"]
    @variables jNO2(t) = 0.0149 [unit = u"s^-1"]
    @variables jCH2Oa(t) = 0.00014 [unit = u"s^-1"]
    # TODO(JL): What's difference between the two photolysis reactions of CH2O, do we really need both? (@variables jCH20b(t) = 0.00014 [unit = u"s^-1"])
    @variables jCH3OOH(t) = 8.9573 * 10.0^-6 [unit = u"s^-1"]

    @parameters k1 = 1.7e-12 [unit = u"(s*ppb)^-1"] T1 = -940 [unit = u"K"]
    @parameters k2 = 1.0e-14 [unit = u"(s*ppb)^-1"] T2 = -490 [unit = u"K"]
    @parameters k3 = 4.8e-11 [unit = u"(s*ppb)^-1"] T3 = 250 [unit = u"K"]
    @parameters k4 = 3.0e-12 [unit = u"(s*ppb)^-1"] T4 = -1500 [unit = u"K"]
    @parameters k5 = 3.5e-12 [unit = u"(s*ppb)^-1"] T5 = 250 [unit = u"K"]
    @parameters k6 = 2.45e-12 [unit = u"(s*ppb)^-1"] T6 = -1775 [unit = u"K"]
    @parameters k7 = 5.5e-12 [unit = u"(s*ppb)^-1"] T7 = 125 [unit = u"K"]
    @parameters k8 = 4.1e-13 [unit = u"(s*ppb)^-1"] T8 = 750 [unit = u"K"]
    @parameters k9 = 2.7e-12 [unit = u"(s*ppb)^-1"] T9 = 200 [unit = u"K"]
    @parameters k10 = 1.1e-12 [unit = u"(s*ppb)^-1"] T10 = 200 [unit = u"K"]
    @parameters k11 = 2.8e-12 [unit = u"(s*ppb)^-1"] T11 = 300 [unit = u"K"]
    @parameters k12 = 9.5e-14 / 10 [unit = u"s^-1*ppb^-19"] T12 = 390 [unit = u"K"]
    @parameters k13 = 1.1e-11 [unit = u"(s*ppb)^-1"] T13 = -240 [unit = u"K"]
    @parameters k14 = 2.7e-11 [unit = u"(s*ppb)^-1"] T14 = 390 [unit = u"K"]
    @parameters k15 = 2.7e-11 / 2 [unit = u"s^-1*ppb^-3"] T15 = 390 [unit = u"K"]
    @parameters k16 = 5.59e-15 [unit = u"(s*ppb)^-1"] T16 = -1814 [unit = u"K"]
    @parameters k17 = 3.0e-13 [unit = u"(s*ppb)^-1"] T17 = 460 [unit = u"K"]
    @parameters k18 = 1.8e-12 [unit = u"(s*ppb)^-1"] k19 = 1.5e-13 [unit = u"(s*ppb)^-1"]
    @parameters t [unit = u"s"]
    @parameters T = 280.0 [unit = u"K"]

    @variables O3(t) = 10.0 [unit = u"ppb"]
    @variables O1d(t) = 0 [unit = u"ppb"]
    @variables OH(t) = 10.0 [unit = u"ppb"]
    @variables HO2(t) = 10.0 [unit = u"ppb"]
    @variables O2(t) = 2.1 * (10.0^8) [unit = u"ppb"]
    @variables H2O(t) = 450.0 [unit = u"ppb"]
    @variables NO(t) = 0.0 [unit = u"ppb"]
    @variables NO2(t) = 10.0 [unit = u"ppb"]
    @variables CH4(t) = 1700.0 [unit = u"ppb"]
    @variables CH3O2(t) = 0.01 [unit = u"ppb"]
    @variables CH2O(t) = 0.15 [unit = u"ppb"]
    @variables CO(t) = 275.0 [unit = u"ppb"]
    @variables CH3OOH(t) = 1.6 [unit = u"ppb"]
    @variables CH3O(t) = 0.0 [unit = u"ppb"]
    @variables DMS(t) = 50 [unit = u"ppb"]
    @variables SO2(t) = 2.0 [unit = u"ppb"]
    @variables ISOP(t) = 0.15 [unit = u"ppb"]
    @variables H2O2(t) = 2.34 [unit = u"ppb"]

    c = 2.46e10 # TODO(JL): What is this constant?
    rate(k, Tc) = k * exp(Tc / T) * c

    # Create reaction system, ignoring aqueous chemistry.
    rxs = [
        #O3 + OH --> HO2 + O2
        Reaction(rate(k1, T1), [O3, OH], [HO2, O2], [1, 1], [1, 1])
        #HO2 + O3 --> 2O2 + OH
        Reaction(rate(k2, T2), [HO2, O3], [O2, OH], [1, 1], [2, 1])
        #HO2 + OH --> H2O + O2
        Reaction(rate(k3, T3), [HO2, OH], [HO2, O2], [1, 1], [1, 1])
        #NO + O3 --> NO2 + O2
        Reaction(rate(k4, T4), [NO, O3], [NO2, O2], [1, 1], [1, 1])
        #HO2 + NO --> NO2 + OH 
        Reaction(rate(k5, T5), [HO2, NO], [NO2, OH], [1, 1], [1, 1])
        #CH4 + OH --> CH3O2 + H2O
        Reaction(rate(k6, T6), [CH4, OH], [CH3O2, H2O], [1, 1], [1, 1])
        #CH2O + OH --> CO + H2O + HO2
        Reaction(rate(k7, T7), [CH2O, OH], [CO, H2O, HO2], [1, 1], [1, 1, 1])
        #CH3O2 + HO2 --> CH3OOH + O2
        Reaction(rate(k8, T8), [CH3O2, HO2], [CH3OOH, O2], [1, 1], [1, 1])
        #CH3OOH + OH --> CH3O2 + H2O
        Reaction(rate(k9, T9), [CH3OOH, OH], [CH3O2, H2O], [1, 1], [1, 1])
        #CH3OOH + OH --> CH3O + H2O + OH
        Reaction(rate(k10, T10), [CH3OOH, OH], [CH3O, H2O, OH], [1, 1], [1, 1, 1])
        #CH3O2 + NO --> CH2O + HO2 + NO2    
        Reaction(rate(k11, T11), [CH3O2, NO], [CH2O, HO2, NO2], [1, 1], [1, 1, 1])
        #10CH3O2 + 10CH3O2 --> 20CH2O + 8HO2
        Reaction(rate(k12, T12), [CH3O2, CH3O2], [CH2O, H2O], [10, 10], [20, 8])
        #DMS + OH --> SO2
        Reaction(rate(k13, T13), [DMS, OH], [SO2], [1, 1], [1])
        #ISOP +OH --> 2CH3O2
        Reaction(rate(k14, T14), [ISOP, OH], [CH3O2], [1, 1], [2])
        #2ISOP + 2OH --> 2ISOP + OH
        Reaction(rate(k15, T15), [ISOP, OH], [ISOP, OH], [2, 2], [2, 1])
        #ISOP + O3 --> 0.87CH2O + 1.86CH3O2 + 0.06HO2 + 0.05CO
        Reaction(
            rate(k16, T16),
            [ISOP, O3],
            [CH2O, CH3O2, HO2, CO],
            [1, 1.0],
            [0.87, 1.86, 0.06, 0.05],
        )
        #O3 -> O2 + O(1D)
        Reaction(jO31D * 10^(-20), [O3], [O1d, O2], [1], [1, 1]) # TODO(JL): Is 10^(-20) a reasonable value?
        #O(1D) + H2O -> 2OH
        Reaction(j2OH, [O1d, H2O], [OH], [1, 1], [2])
        #H2O2 --> 2OH
        Reaction(jH2O2, [H2O2], [OH], [1], [2])
        #NO2 --> NO + O3
        Reaction(jNO2, [NO2], [NO, O3], [1], [1, 1])
        #CH2O --> CO + 2HO2
        Reaction(jCH2Oa, [CH2O], [CO, HO2], [1], [1, 2])
        #CH2O --> CO
        # TODO(JL): What's difference between the two photolysis reactions of CH2O, do we really need both? (Reaction(jCH20b, [CH2O], [CO], [1], [1]))
        #CH3OOH --> CH2O + HO2 + OH
        Reaction(jCH3OOH, [CH3OOH], [CH2O, HO2, OH], [1], [1, 1, 1])
        #HO2 + HO2 = H2O2 + O2
        Reaction(rate(k17, T17), [HO2], [H2O2, O2], [2], [1, 1])
        #OH + H2O2 = H2O + HO2
        Reaction(k18 * c, [OH, H2O2], [H2O, HO2], [1, 1], [1, 1])
        #OH + CO = HO2
        Reaction(k19 * c, [OH, CO], [HO2], [1, 1], [1])
    ]

    @named superfast = ReactionSystem(rxs, t)
    return superfast
end
