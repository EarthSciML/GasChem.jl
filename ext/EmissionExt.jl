module EmissionExt
using GasChem, EarthSciData, ModelingToolkit, Dates, EarthSciMLBase, Unitful

P = 101325 # Pa
Tep = 298.15 # K
R = 8.314 # Pa*m3/(mol*K)
lat = 30.0
lev = 1.0
@constants uu = 1 [unit = u"nmol/mol/s"]
Δz2 = 60.0

Base.:(+)(emis::NEI2016MonthlyEmis, sf::SuperFast) = operator_compose(sf, emis, Dict(
    sf.sys.NO2 => emis.sys.mrggrid_withbeis_withrwc₊NO2 => uu * 2000 / 2240 * 1e15 * R * Tep / (Δz2 * 46.0055 * P * 3600 * 24),
    sf.sys.NO => emis.sys.mrggrid_withbeis_withrwc₊NO => uu * 2000 / 2240 * 1e15 * R * Tep / (Δz2 * 46.0055 * P * 3600 * 24),
    sf.sys.CH2O => emis.sys.mrggrid_withbeis_withrwc₊FORM => uu * 2000 / 2240 * 1e15 * R * Tep / (Δz2 * 46.0055 * P * 3600 * 24),
    sf.sys.CH4 => emis.sys.mrggrid_withbeis_withrwc₊CH4 => uu * 2000 / 2240 * 1e15 * R * Tep / (Δz2 * 46.0055 * P * 3600 * 24),
    sf.sys.CO => emis.sys.mrggrid_withbeis_withrwc₊CO => uu * 2000 / 2240 * 1e15 * R * Tep / (Δz2 * 46.0055 * P * 3600 * 24),
    sf.sys.SO2 => emis.sys.mrggrid_withbeis_withrwc₊SO2 => uu * 2000 / 2240 * 1e15 * R * Tep / (Δz2 * 46.0055 * P * 3600 * 24),
    sf.sys.ISOP => emis.sys.mrggrid_withbeis_withrwc₊ISOP => uu * 2000 / 2240 * 1e15 * R * Tep / (Δz2 * 46.0055 * P * 3600 * 24)
))

Base.:(+)(b::SuperFast, e::NEI2016MonthlyEmis) = e + b
end