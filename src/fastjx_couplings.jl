
function EarthSciMLBase.couple2(c::GEOSChemGasPhaseCoupler, p::FastJXCoupler)
    c, p = c.sys, p.sys
    c = param_to_var(c, :j_3, :j_7, :j_9, :j_10, :j_11)
    ConnectorSystem(
        [c.j_3 ~ p.j_o31D
         c.j_7 ~ p.j_CH2Oa
         c.j_9 ~ p.j_h2o2
         c.j_10 ~ p.j_CH3OOH
         c.j_11 ~ p.j_NO2],
        c,
        p
    )
end

function EarthSciMLBase.couple2(c::SuperFastCoupler, p::FastJXCoupler)
    c, p = c.sys, p.sys
    c = param_to_var(
        convert(ODESystem, c),
        :jH2O2,
        :jNO2,
        :jCH2Oa,
        :jCH2Ob,
        :jCH3OOH,
        :jO32OH
    )
    ConnectorSystem(
        [c.jH2O2 ~ p.j_h2o2
         c.jCH2Oa ~ p.j_CH2Oa
         c.jCH2Ob ~ p.j_CH2Ob
         c.jCH3OOH ~ p.j_CH3OOH
         c.jNO2 ~ p.j_NO2
         c.jO32OH ~ p.j_o32OH],
        c,
        p
    )
end
