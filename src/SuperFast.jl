using Catalyst
using Unitful

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
sol = solve(ODEProblem(rs, [], (0,360), []), Tsit5())
plot(sol)
```
"""
function superfast()
    @parameters jH2O2 = 1.0097 * 10.0^-5 [unit = u"ppb/s"]
    @parameters jNO2 = 0.0149 [unit = u"ppb/s"]
    @parameters jCH20a = 0.00014 [unit = u"ppb/s"]
    @parameters jCH20b = 0.00014 [unit = u"ppb/s"]
    @parameters jCH3OOH = 8.9573 * 10.0^-6 [unit = u"ppb/s"]
    @parameters t [unit = u"s"]
    @parameters T = 280.0 # [unit = u"K"] # Adding units causes a problem when using in exp function.

    @variables O3(t) = 10.0 [unit = u"ppb"]
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

    # Create reaction system, ignoring aqueous chemistry.
    # It is in this form rather than the more compact form because it has non-integer coefficients.
    rxs = [
        #O3 + OH --> HO2 + O2
        Reaction(
            1.7 * 10^(-12) * exp(-940 / T) * 2.46 * 10^10,
            [O3, OH], [HO2, O2], [1, 1], [1, 1]
        )
        #HO2 + O3 --> 2O2 + OH
        Reaction(
            1.0 * 10^(-14) * exp(-490 / T) * 2.46 * 10^10,
            [HO2, O3], [O2, OH], [1, 1], [2, 1]
        )
        #HO2 + OH --> H2O + O2
        Reaction(
            4.8 * 10^(-11) * exp(250 / T) * 2.46 * 10^10,
            [HO2, OH], [HO2, O2], [1, 1], [1, 1]
        )
        #NO + O3 --> NO2 + O2
        Reaction(
            3.0 * 10^(-12) * exp(-1500 / T) * 2.46 * 10^10,
            [NO, O3], [NO2, O2], [1, 1], [1, 1]
        )
        #HO2 + NO --> NO2 + OH 
        Reaction(
            3.5 * 10^(-12) * exp(250 / T) * 2.46 * 10^10,
            [HO2, NO], [NO2, OH], [1, 1], [1, 1]
        )
        #CH4 + OH --> CH3O2 + H2O
        Reaction(
            2.45 * 10^(-12) * exp(-1775 / T) * 2.46 * 10^10,
            [CH4, OH], [CH3O2, H2O], [1, 1], [1, 1]
        )
        #CH2O + OH --> CO + H2O + HO2
        Reaction(
            5.50 * 10^(-12) * exp(125 / T) * 2.46 * 10^10,
            [CH2O, OH], [CO, H2O, HO2], [1, 1], [1, 1, 1]
        )
        #CH3O2 + HO2 --> CH3OOH + O2
        Reaction(
            4.10 * 10^(-13) * exp(750 / T) * 2.46 * 10^10,
            [CH3O2, HO2], [CH3OOH, O2], [1, 1], [1, 1]
        )
        #CH3OOH + OH --> CH3O2 + H2O
        Reaction(
            2.70 * 10^(-12) * exp(200 / T) * 2.46 * 10^10,
            [CH3OOH, OH], [CH3O2, H2O], [1, 1], [1, 1]
        )
        #CH3OOH + OH --> CH3O + H2O + OH
        Reaction(
            1.10 * 10^(-12) * exp(200 / T) * 2.46 * 10^10,
            [CH3OOH, OH], [CH3O, H2O, OH], [1, 1], [1, 1, 1]
        )
        #CH3O2 + NO --> CH2O + HO2 + NO2    
        Reaction(
            2.80 * 10^(-12) * exp(300 / T) * 2.46 * 10^10,
            [CH3O2, NO], [CH2O, HO2, NO2], [1, 1], [1, 1, 1]
        )
        #10CH3O2 + 10CH3O2 --> 20CH2O + 8HO2
        Reaction(
            9.50 * 10^(-14) * exp(390 / T) / 10 * 2.46 * 10^10,
            [CH3O2, CH3O2], [CH2O, H2O], [10, 10], [20, 8]
        )
        #DMS + OH --> SO2
        Reaction(
            1.10 * 10^(-11) * exp(-240 / T) * 2.46 * 10^10,
            [DMS, OH], [SO2], [1, 1], [1]
        )
        #ISOP +OH --> 2CH3O2
        Reaction(
            2.70 * 10^(-11) * exp(390 / T) * 2.46 * 10^10,
            [ISOP, OH], [CH3O2], [1, 1], [2]
        )
        #2ISOP + 2OH --> 2ISOP + OH
        Reaction(
            2.70 * 10^(-11) * exp(390 / T) / 2 * 2.46 * 10^10,
            [ISOP, OH], [ISOP, OH], [2, 2], [2, 1]
        )
        #ISOP + O3 --> 0.87CH2O + 1.86CH3O2 + 0.06HO2 + 0.05CO
        Reaction(
            5.59 * 10^(-15) * exp(-1814 / T) * 2.46 * 10^10,
            [ISOP, O3], [CH2O, CH3O2, HO2, CO], [1, 1.0], [0.87, 1.86, 0.06, 0.05]
        )
        #H2O2 --> 2OH
        Reaction(jH2O2, [H2O2], [OH], [1], [2])
        #NO2 --> NO + O3
        Reaction(jNO2, [NO2], [NO, O3], [1], [1, 1])
        #CH2O --> CO + 2HO2
        Reaction(jCH20a, [CH2O], [CO, HO2], [1], [1, 2])
        #CH2O --> CO
        Reaction(jCH20b, [CH2O], [CO], [1], [1])
        #CH3OOH --> CH2O + HO2 + OH
        Reaction(jCH3OOH, [CH3OOH], [CH2O, HO2, OH], [1], [1, 1, 1])
        #HO2 + HO2 = H2O2 + O2
        Reaction(
            3.0 * 10^(-13) * exp(460 / T) * 2.46 * 10^10,
            [HO2], [H2O2, O2], [2], [1, 1]
        )
        #OH + H2O2 = H2O + HO2
        Reaction(
            1.8 * 10^(-12) * 2.46 * 10^10,
            [OH, H2O2], [H2O, HO2], [1, 1], [1, 1]
        )
        #OH + CO = HO2
        Reaction(
            1.5 * 10^(-13) * 2.46 * 10^10,
            [OH, CO], [HO2], [1, 1], [1]
        )
    ]

    @named superfast = ReactionSystem(rxs, t)
end
