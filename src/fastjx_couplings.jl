function EarthSciMLBase.couple2(c::GEOSChemGasPhaseCoupler, p::FastJXCoupler)
    c, p = c.sys, p.sys
    c = param_to_var(
        c,
        :j_1, :j_2, :j_3, :j_6, :j_7, :j_8, :j_9, :j_10,
        :j_11, :j_12, :j_13, :j_14, :j_15, :j_16, :j_18,
        :j_19, :j_20, :j_22, :j_24, :j_25, :j_26, :j_27, :j_28, :j_30, :j_32, :j_33, :j_34, :j_37,
        :j_38, :j_39, :j_40, :j_41, :j_42, :j_43, :j_44,
        :j_45, :j_46, :j_47, :j_48, :j_49, :j_50, :j_51,
        :j_53, :j_54, :j_55, :j_56, :j_59, :j_61, :j_63,
        :j_64, :j_65, :j_66, :j_68, :j_69, :j_70, :j_71,
        :j_72, :j_73, :j_74, :j_76, :j_77, :j_122, :j_134,
        # organic surrogate channels (photolysis completion)
        :j_78, :j_79, :j_80, :j_81, :j_82, :j_83, :j_84, :j_85,
        :j_86, :j_87, :j_88, :j_89, :j_90, :j_91, :j_92, :j_93,
        :j_94, :j_95, :j_96, :j_97, :j_98, :j_99, :j_105, :j_106,
        :j_107, :j_108, :j_109, :j_110, :j_111, :j_112, :j_113, :j_135,
        :j_136, :j_137, :j_138, :j_139, :j_140, :j_141, :j_142, :j_143,
        :j_144, :j_145, :j_146, :j_147, :j_148, :j_149, :j_150, :j_151,
        :j_152, :j_153, :j_154, :j_155, :j_156, :j_157, :j_158, :j_159,
        :j_160, :j_161, :j_162, :j_164, :j_165, :j_166
    )
    return ConnectorSystem(
        [
            c.j_1 ~ p.j_O2 #
            c.j_2 ~ p.j_O3 #
            c.j_3 ~ p.j_O31D #
            c.j_6 ~ p.j_NO #
            c.j_7 ~ p.j_H2COa
            c.j_8 ~ p.j_H2COb #
            c.j_9 ~ p.j_H2O2 #
            c.j_10 ~ p.j_CH3OOH
            c.j_11 ~ p.j_NO2
            c.j_12 ~ p.j_NO3a #
            c.j_13 ~ p.j_NO3b #
            c.j_14 ~ p.j_N2O5 #
            c.j_15 ~ p.j_HNO2 #
            c.j_16 ~ p.j_HNO3 #
            c.j_18 ~ p.j_HNO4 #
            c.j_19 ~ p.j_ClNO3a #
            c.j_20 ~ p.j_ClNO3b #
            c.j_22 ~ p.j_Cl2 #
            c.j_24 ~ p.j_HOCl #
            c.j_25 ~ p.j_OClO #
            c.j_26 ~ p.j_Cl2O2 #
            c.j_27 ~ p.j_ClO #
            c.j_28 ~ p.j_BrO #
            c.j_30 ~ p.j_BrNO3 #
            c.j_32 ~ p.j_HOBr #
            c.j_33 ~ p.j_BrCl #
            c.j_34 ~ p.j_OCS #
            #c.j_35 ~ p.j_N2O#IS NOT IN geos-chem/KPP/fullchem/fullchem.eqn
            c.j_37 ~ p.j_CFCl3 #
            c.j_38 ~ p.j_CF2Cl2 #
            c.j_39 ~ p.j_F113 #
            c.j_40 ~ p.j_F114 #
            c.j_41 ~ p.j_F115 #
            c.j_42 ~ p.j_CCl4 #
            c.j_43 ~ p.j_CH3Cl #
            c.j_44 ~ p.j_MeCCl3 #
            c.j_45 ~ p.j_CH2Cl2 #
            c.j_46 ~ p.j_CHF2Cl #
            c.j_47 ~ p.j_F123 #
            c.j_48 ~ p.j_F141b #
            c.j_49 ~ p.j_F142b #
            c.j_50 ~ p.j_CH3Br #
            c.j_51 ~ p.j_H1211 #
            c.j_53 ~ p.j_H1301 #
            c.j_54 ~ p.j_H2402 #
            c.j_55 ~ p.j_CH2Br2 #
            c.j_56 ~ p.j_CHBr3 #
            c.j_59 ~ p.j_PAN #
            c.j_61 ~ p.j_ActAld #
            c.j_63 ~ p.j_MeVKa #
            c.j_64 ~ p.j_MeVKb #
            c.j_65 ~ p.j_MeVKc #
            c.j_66 ~ p.j_MeAcr #
            c.j_68 ~ p.j_GlyAld #
            c.j_69 ~ p.j_MEKeto #
            c.j_70 ~ p.j_PrAld #
            c.j_71 ~ p.j_MGlyxl #
            c.j_72 ~ p.j_Glyxla #
            c.j_73 ~ p.j_Glyxlb #
            c.j_74 ~ p.j_Glyxlc #
            c.j_76 ~ p.j_Aceta #
            c.j_77 ~ p.j_Acetb #
            c.j_122 ~ p.j_CH3I #
            c.j_134 ~ p.j_CH3NO3
            # === organic surrogate channels (photolysis completion) ===
            # Hydroperoxides/peroxides photolyze via the /CH3OOH/ surrogate, organic & alkyl
            # nitrates via /CH3NO3/, plus 13 dedicated Cloud-J v7.3e cross-sections. Surrogate
            # assignments + j-factors follow GEOS-Chem's FJX_j2j.dat (v7.3e, index-aligned in the
            # organic block, species-verified). ETP carries GEOS-Chem's 0.5 factor.
            c.j_78 ~ p.j_CH3NO3 # IDN (14.1.1 generic /ONIT2/)
            c.j_79 ~ p.j_CH3OOH # PRPN
            c.j_80 ~ 0.5 * p.j_CH3OOH # ETP
            c.j_81 ~ p.j_CH3OOH # RA3P
            c.j_82 ~ p.j_CH3OOH # RB3P
            c.j_83 ~ p.j_CH3OOH # R4P
            c.j_84 ~ p.j_CH3OOH # PP
            c.j_85 ~ p.j_CH3OOH # RP
            c.j_86 ~ p.j_HMHP # HMHP (dedicated v7.3e σ)
            c.j_87 ~ p.j_CH3OOH # HPETHNL (14.1.1 generic /PrAldP/)
            c.j_88 ~ p.j_MGlyxl # PYAC
            c.j_89 ~ p.j_CH3NO3 # PROPNN (14.1.1 generic /PROPNN/)
            c.j_90 ~ p.j_MGlyxl # MVKHC
            c.j_91 ~ p.j_PrAld # MVKHCB
            c.j_92 ~ p.j_CH3OOH # MVKHP
            c.j_93 ~ p.j_CH3OOH # MVKPC (14.1.1 generic /PrAldP/)
            c.j_94 ~ p.j_ENOL # MCRENOL (dedicated v7.3e σ)
            c.j_95 ~ p.j_CH3OOH # MCRHP (14.1.1 generic /PrAldP/)
            c.j_96 ~ p.j_CH3OOH # MACR1OOH (14.1.1 generic /PrAldP/)
            c.j_97 ~ p.j_CH3OOH # ATOOH
            c.j_98 ~ p.j_CH3NO3 # R4N2
            c.j_99 ~ p.j_CH3OOH # MAP
            c.j_105 ~ p.j_H2O2 # PIP
            c.j_106 ~ p.j_ICN # ICN (dedicated v7.3e σ)
            c.j_107 ~ p.j_ETHLN # ETHLN (dedicated v7.3e σ)
            c.j_108 ~ p.j_MVKN # MVKN (dedicated v7.3e σ)
            c.j_109 ~ p.j_MACRN # MCRHN (dedicated v7.3e σ)
            c.j_110 ~ p.j_MACRNP # MCRHNB (dedicated v7.3e σ)
            c.j_111 ~ p.j_ONIT1 # MONITS (dedicated v7.3e σ)
            c.j_112 ~ p.j_ONIT1 # MONITU (dedicated v7.3e σ)
            c.j_113 ~ p.j_ONIT1 # HONIT (dedicated v7.3e σ)
            c.j_135 ~ p.j_ETNO3 # ETNO3 (dedicated v7.3e σ)
            c.j_136 ~ p.j_IPRNO3 # IPRNO3 (dedicated v7.3e σ)
            c.j_137 ~ p.j_NPRNO3 # NPRNO3 (dedicated v7.3e σ)
            c.j_138 ~ p.j_CH3OOH # RIPA
            c.j_139 ~ p.j_CH3OOH # RIPB
            c.j_140 ~ p.j_CH3OOH # RIPC
            c.j_141 ~ p.j_CH3OOH # RIPD
            c.j_142 ~ p.j_CH3OOH # HPALD1 (14.1.1 generic /HPALD1/)
            c.j_143 ~ p.j_CH3OOH # HPALD2 (14.1.1 generic /HPALD2/)
            c.j_144 ~ p.j_CH3OOH # HPALD3 (14.1.1 generic /PrAldP/)
            c.j_145 ~ p.j_CH3OOH # HPALD4 (14.1.1 generic /PrAldP/)
            c.j_146 ~ p.j_ONIT1 # IHN1 (dedicated v7.3e σ)
            c.j_147 ~ p.j_ONIT1 # IHN2 (dedicated v7.3e σ)
            c.j_148 ~ p.j_ONIT1 # IHN3 (dedicated v7.3e σ)
            c.j_149 ~ p.j_ONIT1 # IHN4 (dedicated v7.3e σ)
            c.j_150 ~ p.j_NITP # INPB (dedicated v7.3e σ)
            c.j_151 ~ p.j_CH3OOH # INPD
            c.j_152 ~ p.j_ONIT1 # INPD (dedicated v7.3e σ)
            c.j_153 ~ p.j_PrAld # ICPDH
            c.j_154 ~ p.j_CH3OOH # ICPDH
            c.j_155 ~ p.j_HP2 # IDHDP (dedicated v7.3e σ)
            c.j_156 ~ p.j_CH3OOH # IDHPE
            c.j_157 ~ p.j_CH3OOH # IDCHP (14.1.1 generic /PrAldP/)
            c.j_158 ~ p.j_CH3OOH # ITHN
            c.j_159 ~ p.j_ONIT1 # ITHN (dedicated v7.3e σ)
            c.j_160 ~ p.j_MACRNP # ITCN (dedicated v7.3e σ)
            c.j_161 ~ p.j_PrAld # ITCN
            c.j_162 ~ p.j_CH3OOH # ETHP
            c.j_164 ~ p.j_CH3OOH # BZCO3H
            c.j_165 ~ p.j_CH3OOH # BENZP
            c.j_166 ~ p.j_CH3NO3 # NPHEN (14.1.1 generic /PROPNN/)
        ],
        c,
        p
    )
end

function EarthSciMLBase.couple2(c::SuperFastCoupler, p::FastJXCoupler)
    c, p = c.sys, p.sys
    c = param_to_var(
        c,
        :jH2O2,
        :jNO2,
        :jH2COa,
        :jH2COb,
        :jCH3OOH,
        :jO32OH
    )
    return ConnectorSystem(
        [
            c.jH2O2 ~ p.j_H2O2
            c.jH2COa ~ p.j_H2COa
            c.jH2COb ~ p.j_H2COb
            c.jCH3OOH ~ p.j_CH3OOH
            c.jNO2 ~ p.j_NO2
            c.jO32OH ~ p.j_o32OH
        ],
        c,
        p
    )
end

function EarthSciMLBase.couple2(c::PolluCoupler, p::FastJXCoupler)
    c, p = c.sys, p.sys
    c = param_to_var(
        c,
        :jNO2_O3P,
        :jH2COa,
        :jH2COb,
        :jALD,
        :jPAN,
        :jO3_O1D,
        :jO3_O3P,
        :jNO3_NO,
        :jNO3_NO2,
        :jN2O5
    )
    ConnectorSystem(
        [c.jNO2_O3P ~ p.j_NO2
         c.jH2COa ~ p.j_H2COa
         c.jH2COb ~ p.j_H2COb
         c.jALD ~ p.j_ActAld
         c.jPAN ~ p.j_PAN
         c.jO3_O1D ~ p.j_O31D
         c.jO3_O3P ~ p.j_O3
         c.jNO3_NO ~ p.j_NO3b
         c.jNO3_NO2 ~ p.j_NO3a
         c.jN2O5 ~ p.j_N2O5],
        c,
        p
    )
end