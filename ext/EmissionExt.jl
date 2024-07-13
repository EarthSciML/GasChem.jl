module EmissionExt
using GasChem, EarthSciData, ModelingToolkit, Dates, EarthSciMLBase, Unitful

@parameters t [unit = u"s"]

@parameters lat 
@parameters lon 
@parameters lev

@parameters Δz  [unit = u"m"]
emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", t, lon, lat, lev, Δz; dtype=Float64)

export register_couplings_ext

# Use a global flag to track initialization and ensure that the coupling between SuperFast and NEI Emission is only initialized once
const couplings_registered_ext = Ref(false)

function register_couplings_ext()
    if couplings_registered_ext[]
        println("Couplings have already been registered.")
        return
    end
    
    println("Registering couplings in ext")

    register_coupling(SuperFast(t), NeiEmissions(t)) do c, e
        operator_compose(convert(ODESystem, c), e, Dict(
            c.NO2 => e.NO2,
            c.NO => e.NO,
            c.CH2O => e.FORM,
            c.CH4 => e.CH4, 
            c.CO => e.CO, 
            c.SO2 => e.SO2, 
            c.ISOP => e.ISOP
        ))
    end
    
    
    register_coupling(NeiEmissions(t), emis) do c, e
        c = param_to_var(convert(ODESystem, c), :r_NO2, :r_NO, :r_FORM, :r_CH4, :r_CO, :r_SO2, :r_ISOP)
        ConnectorSystem([
            c.r_NO2 ~ e.mrggrid_withbeis_withrwc₊NO2 
            c.r_NO ~ e.mrggrid_withbeis_withrwc₊NO 
            c.r_FORM ~ e.mrggrid_withbeis_withrwc₊FORM 
            c.r_CH4 ~ e.mrggrid_withbeis_withrwc₊CH4 
            c.r_CO ~ e.mrggrid_withbeis_withrwc₊CO 
            c.r_SO2 ~ e.mrggrid_withbeis_withrwc₊SO2 
            c.r_ISOP ~ e.mrggrid_withbeis_withrwc₊ISOP 
        ], c, e)
    end
    
    couplings_registered_ext[] = true
    println("Coupling registry after registration: ", EarthSciMLBase.coupling_registry)
end

function __init__()
    println("Initializing EmissionExt module")
    register_couplings_ext()
end
end