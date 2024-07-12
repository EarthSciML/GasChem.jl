module GasChem

using EarthSciMLBase
using ModelingToolkit
using Catalyst
using Dates
using Unitful
using StaticArrays
using Interpolations

include("SuperFast.jl")
include("Fast-JX.jl")
include("geoschem_ratelaws.jl")
include("geoschem_fullchem.jl")

export register_couplings

# Global flag to track initialization to ensure that the coupling is only initialized once
const couplings_registered = Ref(false)

function register_couplings()
    if couplings_registered[]
        println("Couplings have already been registered.")
        return
    end
    
    println("Registering couplings")
    @parameters t [unit = u"s"]

    register_coupling(GEOSChemGasPhase(t), FastJX(t)) do c, p
        println("Registering coupling for GEOSChemGasPhase and FastJX")
        c = param_to_var(c, :j_3, :j_9, :j_11, :j_7, :j_10)
        @constants uconv = 1 [unit = u"s"]
        @constants c_fixme1 = 10^(-21) [unit = u"s"] # FIXME: Suspicious constant
        ConnectorSystem([
            c.j_9 ~ uconv * p.j_h2o2
            c.j_7 ~ uconv * p.j_CH2Oa
            c.j_10 ~ uconv * p.j_CH3OOH
            c.j_11 ~ uconv * p.j_NO2
            c.j_3 ~ c_fixme1 * p.j_o31D
        ], c, p)
    end

    register_coupling(SuperFast(t), FastJX(t)) do c, p
        c = param_to_var(convert(ODESystem, c), :jO31D, :jH2O2, :jNO2, :jCH2Oa, :jCH3OOH)
        ConnectorSystem([
            c.jH2O2 ~ p.j_h2o2
            c.jCH2Oa ~ p.j_CH2Oa
            c.jCH3OOH ~ p.j_CH3OOH
            c.jNO2 ~ p.j_NO2
            c.jO31D ~ p.j_o31D
        ], c, p)
    end

    couplings_registered[] = true
    println("Coupling registry after registration: ", EarthSciMLBase.coupling_registry)
end

function __init__()
    println("Initializing GasChem module")
    register_couplings()
end

end
