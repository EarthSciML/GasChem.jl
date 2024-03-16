module EmissionExt
using GasChem, EarthSciData, ModelingToolkit, Dates, EarthSciMLBase, Unitful

struct Emission <: EarthSciMLODESystem
    sys::ODESystem
    function GasChem.Emission(t)
        P = 101325 # Pa
        Tep = 298.15 # K
        R = 8.314 # Pa*m3/(mol*K)
        A = 1 # m2 asume unit area for the box model

        default_time = DateTime(2016, 5, 1)
		lon = -100.0
        lat = 30.0
        lev = 1.0
        @parameters Δz = 60 [unit = u"m"]
        @parameters uu = 1 [unit = u"nmol/mol/s"]
        @parameters t [unit = u"s"]
        Δz2 = 60.0
        emis = NEI2016MonthlyEmis{Float64}("mrggrid_withbeis_withrwc", t, lon, lat, lev, Δz)
        fs = emis.fileset

        D = Differential(t)
        @variables NO(t) [unit = u"nmol/mol"]
		@variables NO2(t) [unit = u"nmol/mol"]
        @variables CH2O(t) [unit = u"nmol/mol"]
		@variables CH4(t) [unit = u"nmol/mol"]
		@variables CO(t) [unit = u"nmol/mol"]
		@variables SO2(t) [unit = u"nmol/mol"]
		@variables ISOP(t) [unit = u"nmol/mol"]

        emis_vars = Dict{Num,Tuple{String,Float64}}(
        NO => ("NO", 30.01),
        NO2 => ("NO2", 46.0055),
        CH2O => ("FORM", 30.0260),
        CH4 => ("CH4", 16.0425),
        CO => ("CO", 28.0101),
        SO2 => ("SO2", 64.0638),
        ISOP => ("ISOP", 68.12),
        )

        eqs = [
            D(NO) ~  uu * EarthSciData.interp!(EarthSciData.DataSetInterpolator{Float64}(fs, "NO", default_time), date(t), lon, lat, lev)* 2000 / 2240 * 1e15 * R * Tep / (Δz2 * emis_vars[NO][2] * P * 3600 * 24)
            D(NO2) ~ uu * EarthSciData.interp!(EarthSciData.DataSetInterpolator{Float64}(fs, "NO2", default_time), date(t), lon, lat, lev)* 2000 / 2240 * 1e15 * R * Tep / (Δz2 * emis_vars[NO2][2] * P * 3600 * 24)
			D(CH2O) ~ uu * EarthSciData.interp!(EarthSciData.DataSetInterpolator{Float64}(fs, "FORM", default_time), date(t), lon, lat, lev)* 2000 / 2240 * 1e15 * R * Tep / (Δz2 * emis_vars[CH2O][2] * P * 3600 * 24)
			D(CH4) ~ uu * EarthSciData.interp!(EarthSciData.DataSetInterpolator{Float64}(fs, "CH4", default_time), date(t), lon, lat, lev)* 2000 / 2240 * 1e15 * R * Tep / (Δz2 * emis_vars[CH4][2] * P * 3600 * 24)
			D(CO) ~ uu * EarthSciData.interp!(EarthSciData.DataSetInterpolator{Float64}(fs, "CO", default_time), date(t), lon, lat, lev)* 2000 / 2240 * 1e15 * R * Tep / (Δz2 * emis_vars[CO][2] * P * 3600 * 24)
			D(SO2) ~ uu * EarthSciData.interp!(EarthSciData.DataSetInterpolator{Float64}(fs, "SO2", default_time), date(t), lon, lat, lev)* 2000 / 2240 * 1e15 * R * Tep / (Δz2 * emis_vars[SO2][2] * P * 3600 * 24)
			D(ISOP) ~ uu * EarthSciData.interp!(EarthSciData.DataSetInterpolator{Float64}(fs, "ISOP", default_time), date(t), lon, lat, lev)* 2000 / 2240 * 1e15 * R * Tep / (Δz2 * emis_vars[ISOP][2] * P * 3600 * 24)
        ]

        new(ODESystem(eqs, t, [NO, NO2, CH2O, CH4, CO, SO2, ISOP], [Δz, uu]; name=:Emission))
    end
end 

Base.:(+)(e::Emission, b::SuperFast) = operator_compose(b, e)
Base.:(+)(b::SuperFast, e::Emission) = e + b
end