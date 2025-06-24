
function EarthSciMLBase.couple2(c::GEOSChemGasPhaseCoupler, p::FastJXCoupler)
    c, p = c.sys, p.sys
    c = param_to_var(c, 
    :j_1, :j_2, :j_3, :j_6, :j_7, :j_8, :j_9, :j_10, :j_11, :j_12, :j_13, :j_14, :j_15, :j_16, :j_18, 
    :j_19, :j_20, :j_22, :j_24, :j_25, :j_26, :j_27, :j_28, :j_30, :j_32, :j_33, :j_34, :j_37, 
    :j_38, :j_39, :j_40, :j_41, :j_42, :j_43, :j_44, :j_45, :j_46, :j_47, :j_48, :j_49, :j_50, :j_51, 
    :j_53, :j_54, :j_55, :j_56, :j_59, :j_61, :j_63, :j_64, :j_65, :j_66, :j_68, :j_69, :j_70, :j_71, 
    :j_72, :j_73, :j_74, :j_76, :j_77, :j_122, :j_134 )
    ConnectorSystem(
        [
        c.j_1 ~ p.j_O2#
        c.j_2 ~ p.j_O3#
        c.j_3 ~ p.j_O31D#
        c.j_6 ~ p.j_NO#
        c.j_7 ~ p.j_H2COa
        c.j_8 ~ p.j_H2COb#
        c.j_9 ~ p.j_H2O2#
        c.j_10 ~ p.j_CH3OOH
        c.j_11 ~ p.j_NO2
        c.j_12 ~ p.j_NO3a#
        c.j_13 ~ p.j_NO3b#
        c.j_14 ~ p.j_N2O5#
        c.j_15 ~ p.j_HNO2#
        c.j_16 ~ p.j_HNO3#
        c.j_18 ~ p.j_HNO4#
        c.j_19 ~ p.j_ClNO3a#
        c.j_20 ~ p.j_ClNO3b#
        c.j_22 ~ p.j_Cl2#
        c.j_24 ~ p.j_HOCl#
        c.j_25 ~ p.j_OClO#
        c.j_26 ~ p.j_Cl2O2#
        c.j_27 ~ p.j_ClO#
        c.j_28 ~ p.j_BrO#
        c.j_30 ~ p.j_BrNO3#
        c.j_32 ~ p.j_HOBr#
        c.j_33 ~ p.j_BrCl#
        c.j_34 ~ p.j_OCS#
        #c.j_35 ~ p.j_N2O#IS NOT IN geos-chem/KPP/fullchem/fullchem.eqn
        c.j_37 ~ p.j_CFCl3#
        c.j_38 ~ p.j_CF2Cl2#
        c.j_39 ~ p.j_F113#
        c.j_40 ~ p.j_F114#
        c.j_41 ~ p.j_F115#
        c.j_42 ~ p.j_CCl4#
        c.j_43 ~ p.j_CH3Cl#
        c.j_44 ~ p.j_MeCCl3#
        c.j_45 ~ p.j_CH2Cl2#
        c.j_46 ~ p.j_CHF2Cl#
        c.j_47 ~ p.j_F123#
        c.j_48 ~ p.j_F141b#
        c.j_49 ~ p.j_F142b#
        c.j_50 ~ p.j_CH3Br#
        c.j_51 ~ p.j_H1211#
        c.j_53 ~ p.j_H1301#
        c.j_54 ~ p.j_H2402#
        c.j_55 ~ p.j_CH2Br2#
        c.j_56 ~ p.j_CHBr3#
        c.j_59 ~ p.j_PAN#
        c.j_61 ~ p.j_ActAld#
        c.j_63 ~ p.j_MeVKa#
        c.j_64 ~ p.j_MeVKb#
        c.j_65 ~ p.j_MeVKc#
        c.j_66 ~ p.j_MeAcr#
        c.j_68 ~ p.j_GlyAld#
        c.j_69 ~ p.j_MEKeto#
        c.j_70 ~ p.j_PrAld#
        c.j_71 ~ p.j_MGlyxl#
        c.j_72 ~ p.j_Glyxla#
        c.j_73 ~ p.j_Glyxlb#
        c.j_74 ~ p.j_Glyxlc#
        c.j_76 ~ p.j_Aceta#
        c.j_77 ~ p.j_Acetb#
        c.j_122 ~ p.j_CH3I#
        c.j_134 ~ p.j_CH3NO3#
        ],
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
