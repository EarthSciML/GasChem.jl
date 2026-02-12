export Pollu

struct PolluCoupler
    sys::Any
end

"""
    Pollu()

This atmospheric chemical system model is built based on this paper:

"GAUSS-SEIDEL ITERATION FOR STIFF ODES FROM CHEMICAL KINETICS" (1994), J. G. VERWER.

The input of the function is Temperature, concentrations of all chemicals, and reaction rates of photolysis reactions

If the keyword argument `rxn_sys` is set to `true`, the function will return a reaction system instead of an ODE system.

# Example

```
using GasChem, EarthSciMLBase, DifferentialEquations, Plots
rs = Pollu()
sol = solve(ODEProblem(structural_simplify(rs), [], (0,360), [], combinatoric_ratelaws=false), AutoTsit5(Rosenbrock23()), saveat=10.0)
plot(sol)
```
"""
function Pollu(; name = :Pollu, rxn_sys = false)
    rx_sys = @reaction_network Pollu begin
        @ivs t [unit = u"s"]
        @parameters(jNO2_O3P=0.35/60.0,
            [unit=u"s^-1"],
            k2=0.266e2/(60.0*1e3),
            [unit=u"1/(ppb*s)"],
            k3=0.123e5/(60.0*1e3), # Original value is 0.120e5/(60.0*1e3), but the value is updated to 0.123e5/(60.0*1e3) to match the Fontana codes "https://archimede.uniba.it/~testset/src/problems/pollu.f"
            [unit=u"1/(ppb*s)"],
            jH2COa=0.86e-3/60.0,
            [unit=u"s^-1"],
            jH2COb=0.82e-3/60.0,
            [unit=u"s^-1"],
            k6=0.15e5/(60.0*1e3),
            [unit=u"1/(ppb*s)"],
            jALD=0.13e-3/60.0,
            [unit=u"s^-1"],
            k8=0.24e5/(60.0*1e3),
            [unit=u"1/(ppb*s)"],
            k9=0.165e5/(60.0*1e3),
            [unit=u"1/(ppb*s)"],
            k10=0.9e4/(60.0*1e3),
            [unit=u"1/(ppb*s)"],
            jPAN=0.22e-1/60.0,
            [unit=u"s^-1"],
            k12=0.12e5/(60.0*1e3),
            [unit=u"1/(ppb*s)"],
            k13=1.88/60.0,
            [unit=u"1/s"],
            k14=0.163e5/(60.0*1e3),
            [unit=u"1/(ppb*s)"],
            k15=0.48e7/60.0,
            [unit=u"s^-1"],
            jO3_O1D=0.35e-3/60.0,
            [unit=u"s^-1"],
            jO3_O3P=0.175e-1/60.0,
            [unit=u"s^-1"],
            k18=0.1e9/(60.0),
            [unit=u"s^-1"],
            k19=0.444e12/60.0,
            [unit=u"s^-1"],
            k20=0.124e4/(60.0*1e3),
            [unit=u"1/(ppb*s)"],
            jNO3_NO=0.21e1/60.0,
            [unit=u"s^-1"],
            jNO3_NO2=0.578e1/60.0,
            [unit=u"s^-1"],
            k23=0.474e-1/(60.0*1e3),
            [unit=u"1/(ppb*s)"],
            k24=0.178e4/(60.0*1e3),
            [unit=u"1/(ppb*s)"],
            jN2O5=0.312e1/60.0,
            [unit=u"s^-1"],) # The original reaction rate units are in 1/min and 1/(ppm*min).

        @species(NO2(t)=4e-4,
            [unit=u"ppb"],
            NO(t)=4e-4,
            [unit=u"ppb"],
            O3P(t)=0,
            [unit=u"ppb", description="ground-state atomic oxygen"],
            O3(t)=40.0,
            [unit=u"ppb"],
            HO2(t)=4e-6,
            [unit=u"ppb"],
            OH(t)=4e-6,
            [unit=u"ppb"],
            CH2O(t)=4e-6,
            [unit=u"ppb"],
            CO(t)=100,
            [unit=u"ppb"],
            ALD(t)=1e-11,
            [unit=u"ppb", description="lumped non-formaldehyde aldehydes"],
            MEO2(t)=4e-6,
            [unit=u"ppb", description="CH3O2 (methyl peroxy radical)"],
            C2O3(t)=0,
            [unit=u"ppb", description="CH3C(O)O2 (acetyl peroxy radical)"],
            CO2(t)=3.55e5,
            [unit=u"ppb"],
            PAN(t)=1e-11,
            [unit=u"ppb", description="CH3C(O)O2NO2 (peroxy nitrate)"],
            CH3O(t)=0,
            [unit=u"ppb", description="CH3O(methoxy radical)"],
            HNO3(t)=4e-6,
            [unit=u"ppb"],
            O1D(t)=0,
            [unit=u"ppb", description="excited-state atomic oxygen"],
            SO2(t)=1e-11,
            [unit=u"ppb"],
            SO4(t)=1e-11,
            [unit=u"ppb", description="SO4 (sulfate)"],
            NO3(t)=4e-6,
            [unit=u"ppb"],
            N2O5(t)=4e-6,
            [unit=u"ppb"],
            )

        #Gas-phase reactions
        k2, NO + O3 --> NO2
        k3, HO2 + NO --> NO2 + OH
        k6, CH2O + OH --> HO2 + CO
        k8, ALD + OH --> C2O3
        k9, C2O3 + NO --> NO2 + MEO2 + CO2
        k10, C2O3 + NO2 --> PAN
        k12, MEO2 + NO --> CH3O + NO2
        k13, CH3O --> CH2O + HO2
        k14, NO2 + OH --> HNO3
        k15, O3P --> O3
        k18, O1D --> 2OH
        k19, O1D --> O3P
        k20, SO2 + OH --> SO4 + HO2
        k23, NO2 + O3 --> NO3
        k24, NO3 + NO2 --> N2O5

        #photolysis reactions
        jNO2_O3P, NO2 --> NO + O3P
        jH2COa, CH2O --> 2HO2 + CO
        jH2COb, CH2O --> CO
        jALD, ALD --> MEO2 + HO2 + CO
        jPAN, PAN --> C2O3 + NO2
        jO3_O1D, O3 --> O1D
        jO3_O3P, O3 --> O3P
        jNO3_NO, NO3 --> NO
        jNO3_NO2, NO3 --> NO2 + O3P
        jN2O5, N2O5 --> NO3 + NO2

    end
    rxns = rx_sys
    if rxn_sys
        return rxns
    end
    # We set `combinatoric_ratelaws=false` because we are modeling macroscopic rather than microscopic behavior. 
    # See [here](https://docs.juliahub.com/ModelingToolkit/Qmdqu/3.14.0/systems/ReactionSystem/#ModelingToolkit.oderatelaw) 
    # and [here](https://github.com/SciML/Catalyst.jl/issues/311).
    convert(
        Catalyst.ReactionRateSystem,
        complete(rxns);
        combinatoric_ratelaws = false,
        name = name,
        metadata = Dict(CoupleType => PolluCoupler)
    )
end
