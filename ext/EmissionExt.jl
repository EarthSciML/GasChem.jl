module EmissionExt
using GasChem, EarthSciData, ModelingToolkit, Dates, EarthSciMLBase, Unitful

@parameters t [unit = u"s"]

@parameters lat = 40 
@parameters lon = -97 
@parameters lev = 1

@parameters Δz = 60 [unit = u"m"]
emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", t, lon, lat, lev, Δz; dtype=Float64)


P = 101325 # Pa
Tep = 298.15 # K
RR = 8.314 # Pa*m3/(mol*K)
Δz2 = 60.0 # m

emis_vars = Dict{String,Float64}(
            "NO" => 30.01,
            "NO2" => 46.0055,
            "FORM" => 30.0260,
            "CH4" => 16.0425,
            "CO" => 28.0101,
            "SO2" => 64.0638,
            "ISOP" => 68.12,
            )

register_coupling(SuperFast(t), emis) do c, e
    operator_compose(convert(ODESystem, c), e, Dict(
        c.NO2 => e.mrggrid_withbeis_withrwc₊NO2 => 2000 / 2240 * 1e12 * RR * Tep / (Δz2 * emis_vars["NO2"] * P),
        c.NO => e.mrggrid_withbeis_withrwc₊NO => 2000 / 2240 * 1e12 * RR * Tep / (Δz2 * emis_vars["NO"] * P),
        c.CH2O => e.mrggrid_withbeis_withrwc₊FORM => 2000 / 2240 * 1e12 * RR * Tep / (Δz2 * emis_vars["FORM"] * P),
        c.CH4 => e.mrggrid_withbeis_withrwc₊CH4 => 2000 / 2240 * 1e12 * RR * Tep / (Δz2 * emis_vars["CH4"] * P), 
        c.CO => e.mrggrid_withbeis_withrwc₊CO => 2000 / 2240 * 1e12 * RR * Tep / (Δz2 * emis_vars["CO"] * P), 
        c.SO2 => e.mrggrid_withbeis_withrwc₊SO2 => 2000 / 2240 * 1e12 * RR * Tep / (Δz2 * emis_vars["SO2"] * P), 
        c.ISOP => e.mrggrid_withbeis_withrwc₊ISOP => 2000 / 2240 * 1e12 * RR * Tep / (Δz2 * emis_vars["ISOP"] * P)
    ))
end