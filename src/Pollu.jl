export Pollu

struct PolluCoupler
    sys::Any
end

"""
    Pollu()

This atmospheric chemical system model is built based on the POLLU benchmark problem from:

"GAUSS-SEIDEL ITERATION FOR STIFF ODES FROM CHEMICAL KINETICS" (1994), J. G. Verwer.

It models a 20-species, 25-reaction atmospheric chemistry system with units of ppb for
concentrations and s for time.

If the keyword argument `rxn_sys` is set to `true`, the function will return a reaction system instead of an ODE system.

# Example

```
using GasChem, EarthSciMLBase
rs = mtkcompile(Pollu())
```
"""
function Pollu(; name = :Pollu, rxn_sys = false)
    rx_sys = @reaction_network Pollu begin
        @ivs t [unit = u"s"]
        @parameters begin
            jNO2_O3P = 0.35 / 60.0,
            [unit = u"s^-1", description = "NO2 photolysis rate (R1)"]
            k2 = 0.266e2 / (60.0 * 1e3),
            [unit = u"1/(ppb*s)", description = "NO + O3 rate (R2)"]
            k3 = 0.123e5 / (60.0 * 1e3),
            [unit = u"1/(ppb*s)", description = "HO2 + NO rate (R3)"]
            jH2COa = 0.86e-3 / 60.0,
            [unit = u"s^-1", description = "CH2O photolysis radical channel (R4)"]
            jH2COb = 0.82e-3 / 60.0,
            [unit = u"s^-1", description = "CH2O photolysis molecular channel (R5)"]
            k6 = 0.15e5 / (60.0 * 1e3),
            [unit = u"1/(ppb*s)", description = "CH2O + OH rate (R6)"]
            jALD = 0.13e-3 / 60.0,
            [unit = u"s^-1", description = "ALD photolysis rate (R7)"]
            k8 = 0.24e5 / (60.0 * 1e3),
            [unit = u"1/(ppb*s)", description = "ALD + OH rate (R8)"]
            k9 = 0.165e5 / (60.0 * 1e3),
            [unit = u"1/(ppb*s)", description = "C2O3 + NO rate (R9)"]
            k10 = 0.9e4 / (60.0 * 1e3),
            [unit = u"1/(ppb*s)", description = "C2O3 + NO2 rate (R10)"]
            jPAN = 0.22e-1 / 60.0,
            [unit = u"s^-1", description = "PAN photolysis rate (R11)"]
            k12 = 0.12e5 / (60.0 * 1e3),
            [unit = u"1/(ppb*s)", description = "MEO2 + NO rate (R12)"]
            k13 = 1.88 / 60.0,
            [unit = u"1/s", description = "CH3O thermal decomposition rate (R13)"]
            k14 = 0.163e5 / (60.0 * 1e3),
            [unit = u"1/(ppb*s)", description = "NO2 + OH rate (R14)"]
            k15 = 0.48e7 / 60.0,
            [unit = u"s^-1", description = "O3P + O2 + M recombination rate (R15)"]
            jO3_O1D = 0.35e-3 / 60.0,
            [unit = u"s^-1", description = "O3 photolysis to O1D (R16)"]
            jO3_O3P = 0.175e-1 / 60.0,
            [unit = u"s^-1", description = "O3 photolysis to O3P (R17)"]
            k18 = 0.1e9 / (60.0), [unit = u"s^-1", description = "O1D + H2O rate (R18)"]
            k19 = 0.444e12 / 60.0,
            [unit = u"s^-1", description = "O1D quenching rate (R19)"]
            k20 = 0.124e4 / (60.0 * 1e3),
            [unit = u"1/(ppb*s)", description = "SO2 + OH rate (R20)"]
            jNO3_NO = 0.21e1 / 60.0,
            [unit = u"s^-1", description = "NO3 photolysis to NO (R21)"]
            jNO3_NO2 = 0.578e1 / 60.0,
            [unit = u"s^-1", description = "NO3 photolysis to NO2 (R22)"]
            k23 = 0.474e-1 / (60.0 * 1e3),
            [unit = u"1/(ppb*s)", description = "NO2 + O3 rate (R23)"]
            k24 = 0.178e4 / (60.0 * 1e3),
            [unit = u"1/(ppb*s)", description = "NO3 + NO2 rate (R24)"]
            jN2O5 = 0.312e1 / 60.0,
            [unit = u"s^-1", description = "N2O5 photolysis rate (R25)"]
        end # The original reaction rate units are in 1/min and 1/(ppm*min).

        @species begin
            NO2(t) = 4e-4, [unit = u"ppb", description = "nitrogen dioxide"]
            NO(t) = 4e-4, [unit = u"ppb", description = "nitric oxide"]
            O3P(t) = 0, [unit = u"ppb", description = "ground-state atomic oxygen"]
            O3(t) = 40.0, [unit = u"ppb", description = "ozone"]
            HO2(t) = 4e-6, [unit = u"ppb", description = "hydroperoxyl radical"]
            OH(t) = 4e-6, [unit = u"ppb", description = "hydroxyl radical"]
            CH2O(t) = 4e-6, [unit = u"ppb", description = "formaldehyde"]
            CO(t) = 100, [unit = u"ppb", description = "carbon monoxide"]
            ALD(t) = 1e-11,
            [unit = u"ppb", description = "lumped non-formaldehyde aldehydes"]
            MEO2(t) = 4e-6, [unit = u"ppb", description = "CH3O2 (methyl peroxy radical)"]
            C2O3(t) = 0, [unit = u"ppb", description = "CH3C(O)O2 (acetyl peroxy radical)"]
            CO2(t) = 3.55e5, [unit = u"ppb", description = "carbon dioxide"]
            PAN(t) = 1e-11, [unit = u"ppb", description = "CH3C(O)O2NO2 (peroxy nitrate)"]
            CH3O(t) = 0, [unit = u"ppb", description = "CH3O (methoxy radical)"]
            HNO3(t) = 4e-6, [unit = u"ppb", description = "nitric acid"]
            O1D(t) = 0, [unit = u"ppb", description = "excited-state atomic oxygen"]
            SO2(t) = 1e-11, [unit = u"ppb", description = "sulfur dioxide"]
            SO4(t) = 1e-11, [unit = u"ppb", description = "SO4 (sulfate)"]
            NO3(t) = 4e-6, [unit = u"ppb", description = "nitrate radical"]
            N2O5(t) = 4e-6, [unit = u"ppb", description = "dinitrogen pentoxide"]
        end

        #Gas-phase reactions
        k2, NO + O3 --> NO2                       # R2
        k3, HO2 + NO --> NO2 + OH                 # R3
        k6, CH2O + OH --> HO2 + CO                # R6
        k8, ALD + OH --> C2O3                      # R8
        k9, C2O3 + NO --> NO2 + MEO2 + CO2        # R9
        k10, C2O3 + NO2 --> PAN                    # R10
        k12, MEO2 + NO --> CH3O + NO2              # R12
        k13, CH3O --> CH2O + HO2                   # R13
        k14, NO2 + OH --> HNO3                     # R14
        k15, O3P --> O3                            # R15
        k18, O1D --> 2OH                           # R18
        k19, O1D --> O3P                           # R19
        k20, SO2 + OH --> SO4 + HO2               # R20
        k23, NO2 + O3 --> NO3                      # R23
        k24, NO3 + NO2 --> N2O5                    # R24

        #photolysis reactions
        jNO2_O3P, NO2 --> NO + O3P                # R1
        jH2COa, CH2O --> 2HO2 + CO                # R4
        jH2COb, CH2O --> CO                        # R5
        jALD, ALD --> MEO2 + HO2 + CO             # R7
        jPAN, PAN --> C2O3 + NO2                   # R11
        jO3_O1D, O3 --> O1D                        # R16
        jO3_O3P, O3 --> O3P                        # R17
        jNO3_NO, NO3 --> NO                        # R21
        jNO3_NO2, NO3 --> NO2 + O3P               # R22
        jN2O5, N2O5 --> NO3 + NO2                  # R25
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
