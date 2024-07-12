module EmissionExt
using GasChem, EarthSciData, ModelingToolkit, Dates, EarthSciMLBase, Unitful

@parameters t [unit = u"s"]

@parameters lat 
@parameters lon 
@parameters lev

@parameters Δz  [unit = u"m"]
emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", t, lon, lat, lev, Δz; dtype=Float64)


P = 101325 # Pa
Tep = 298.15 # K
RR = 8.314 # Pa*m3/(mol*K)

emis_vars = Dict{String,Float64}(
            "NO" => 30.01,
            "NO2" => 46.0055,
            "FORM" => 30.0260,
            "CH4" => 16.0425,
            "CO" => 28.0101,
            "SO2" => 64.0638,
            "ISOP" => 68.12,
            )

export register_couplings_ext

# Use a global flag to track initialization and ensure that the coupling between SuperFast and NEI Emission is only initialized once
const couplings_registered_ext = Ref(false)

function register_couplings_ext()
    if couplings_registered_ext[]
        println("Couplings have already been registered.")
        return
    end
    
    println("Registering couplings in ext")

    register_coupling(SuperFast(t), emis) do c, e
        operator_compose(convert(ODESystem, c), e, Dict(
            c.NO2 => e.mrggrid_withbeis_withrwc₊NO2 => 2000 / 2240 * 1e12 * RR * Tep / (emis_vars["NO2"] * P),
            c.NO => e.mrggrid_withbeis_withrwc₊NO => 2000 / 2240 * 1e12 * RR * Tep / (emis_vars["NO"] * P),
            c.CH2O => e.mrggrid_withbeis_withrwc₊FORM => 2000 / 2240 * 1e12 * RR * Tep / (emis_vars["FORM"] * P),
            c.CH4 => e.mrggrid_withbeis_withrwc₊CH4 => 2000 / 2240 * 1e12 * RR * Tep / (emis_vars["CH4"] * P), 
            c.CO => e.mrggrid_withbeis_withrwc₊CO => 2000 / 2240 * 1e12 * RR * Tep / (emis_vars["CO"] * P), 
            c.SO2 => e.mrggrid_withbeis_withrwc₊SO2 => 2000 / 2240 * 1e12 * RR * Tep / (emis_vars["SO2"] * P), 
            c.ISOP => e.mrggrid_withbeis_withrwc₊ISOP => 2000 / 2240 * 1e12 * RR * Tep / (emis_vars["ISOP"] * P)
        ))
    end # The equation on the far right is used to convert the emission rate from "kg/s/m^3" to "ppb/s"
    
    couplings_registered_ext[] = true
    println("Coupling registry after registration: ", EarthSciMLBase.coupling_registry)
end

function __init__()
    println("Initializing EmissionExt module")
    register_couplings_ext()
end
end