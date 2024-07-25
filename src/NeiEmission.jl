export NeiEmissions

"""
This is a box model used to couple with the emis model (from the US National Emissions Inventory) to convert the emission rate to the ODE system of species shared in the SuperFast mechanism.
Build NeiEmissions model to couple with SuperFast model and emis model
# Example
``` julia
    using GasChem, EarthSciData    
    @parameters t [unit = u"s"]
    @parameters lat = 40 
    @parameters lon = -97 
    @parameters lev = 1
    @parameters Δz = 60 [unit = u"m"]
    emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", t, lon, lat, lev, Δz; dtype=Float64)

    model = couple(SuperFast(t), NeiEmissions(t), emis)
```
"""
function NeiEmissions(t)

    emis_vars = Dict{String,Float64}(
            "NO" => 30.01,
            "NO2" => 46.0055,
            "FORM" => 30.0260,
            "CH4" => 16.0425,
            "CO" => 28.0101,
            "SO2" => 64.0638,
            "ISOP" => 68.12,
            )

    @variables NO2(t) [unit = u"nmol/mol"]
    @variables NO(t) [unit = u"nmol/mol"]
    @variables FORM(t) [unit = u"nmol/mol"]
    @variables CH4(t) [unit = u"nmol/mol"]
    @variables CO(t) [unit = u"nmol/mol"]
    @variables SO2(t) [unit = u"nmol/mol"]
    @variables ISOP(t) [unit = u"nmol/mol"]

    @parameters u_conv = 1 [unit = u"nmol/mol*m^3/kg"]
    @parameters P = 101325 # Pa
    @parameters Tep = 298.15 # K
    @parameters RR = 8.314 # Pa*m3/(mol*K)
    @parameters Δz2 = 60.0 # m  
    @parameters r_NO2 [unit = u"kg/m^3/s"]
    @parameters r_NO [unit = u"kg/m^3/s"]
    @parameters r_FORM [unit = u"kg/m^3/s"]
    @parameters r_CH4 [unit = u"kg/m^3/s"]
    @parameters r_CO [unit = u"kg/m^3/s"]
    @parameters r_SO2 [unit = u"kg/m^3/s"]
    @parameters r_ISOP [unit = u"kg/m^3/s"]

    eqs = [
        Differential(t)(NO2) ~ r_NO2 * u_conv  * 1e12 * RR * Tep / (emis_vars["NO2"] * P),
        Differential(t)(NO) ~ r_NO * u_conv  * 1e12 * RR * Tep / (emis_vars["NO"] * P),
        Differential(t)(FORM) ~ r_FORM * u_conv  * 1e12 * RR * Tep / (emis_vars["FORM"] * P),
        Differential(t)(CH4) ~ r_CH4 * u_conv  * 1e12 * RR * Tep / (emis_vars["CH4"] * P),
        Differential(t)(CO) ~ r_CO * u_conv  * 1e12 * RR * Tep / (emis_vars["CO"] * P),
        Differential(t)(SO2) ~ r_SO2 * u_conv  * 1e12 * RR * Tep / (emis_vars["SO2"] * P),
        Differential(t)(ISOP) ~ r_ISOP * u_conv  * 1e12 * RR * Tep / (emis_vars["ISOP"] * P)
    ]

    ODESystem(eqs, t, [NO2, NO, FORM, CH4, CO, SO2, ISOP], [r_NO2, r_NO, r_FORM, r_CH4, r_CO, r_SO2, r_ISOP, P, Tep, RR, Δz2, u_conv]; name=:new_emis)
end