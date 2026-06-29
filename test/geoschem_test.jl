@testsnippet GEOSChemGasPhaseSetup begin
    using GasChem, EarthSciMLBase
    using OrdinaryDiffEqRosenbrock
    using ModelingToolkit

    tspan = (0.0, 360.0)
    sys = GEOSChemGasPhase()
    sys = mtkcompile(sys)
end

# Unit Test 0: Base case
@testitem "Base case" setup = [GEOSChemGasPhaseSetup] begin
    prob = ODEProblem(sys, [], tspan)
    sol = solve(prob, Rosenbrock23())
    @test sol[sys.O3][end] > 0
end

# Unit Test 1: O1D sensitivity to O3
@testitem "O1D sensitivity to O3" setup = [GEOSChemGasPhaseSetup] begin
    u_1 = 8.160593694128693e-7

    @unpack O3, O1D = sys
    o1 = solve(ODEProblem(sys, [O3 => 20, O1D => 1.0e-6], tspan), Rosenbrock23())
    o2 = solve(ODEProblem(sys, [O3 => 20, O1D => 1.0e-6 * 1.1], tspan), Rosenbrock23())
    test1 = o1[O3][end] - o2[O3][end]

    @test test1 ≈ u_1 rtol = 0.001
end

# Unit Test 2: OH sensitivity to O3
@testitem "OH sensitivity to O3" setup = [GEOSChemGasPhaseSetup] begin
    # Finite-difference O3 sensitivity to initial OH. This is a near-cancellation
    # quantity (a small difference of two O3 ~ 20 endpoints), so reproducibility
    # needs two things: (1) a tight solver tolerance (reltol=1e-10) so each O3
    # endpoint is accurate well below the difference, and (2) a perturbation large
    # enough that the resulting O3 difference sits far above the solver accuracy
    # floor (reltol * O3 ~ 2e-9). A 50% OH perturbation puts the signal ~80x above
    # the floor (vs ~8x for the old 5%, where it was near-floor and drifted ~2x with
    # dependency updates); the response is linear, so this is still a clean
    # sensitivity. Converged to <0.06% between reltol 1e-10 and 1e-11, so rtol=0.01
    # has ample margin while still catching any real chemistry regression.
    u_2 = 1.6170403327464555e-7

    @unpack O3, OH = sys
    o1 = solve(ODEProblem(sys, [O3 => 20, OH => 4.0e-6], tspan), Rosenbrock23(), abstol = 1.0e-10, reltol = 1.0e-10)
    o2 = solve(ODEProblem(sys, [O3 => 20, OH => 4.0e-6 * 1.5], tspan), Rosenbrock23(), abstol = 1.0e-10, reltol = 1.0e-10)
    test2 = o1[O3][end] - o2[O3][end]

    @test test2 ≈ u_2 rtol = 0.01
end

# Unit Test 3: NO2 sensitivity to O3
@testitem "NO2 sensitivity to O3" setup = [GEOSChemGasPhaseSetup] begin
    u_3 = 1.1784680253867919e-7

    @unpack O3, NO2 = sys
    o1 = solve(ODEProblem(sys, [O3 => 20, NO2 => 4.0e-4], tspan), Rosenbrock23(), abstol = 1.0e-6, reltol = 1.0e-6)
    o2 = solve(ODEProblem(sys, [O3 => 20, NO2 => 4.0e-4 * 1.05], tspan), Rosenbrock23(), abstol = 1.0e-6, reltol = 1.0e-6)
    test3 = o1[O3][end] - o2[O3][end]

    @test test3 ≈ u_3 rtol = 0.01
end

# Unit Test 4: HO2 sensitivity to O3
@testitem "HO2 sensitivity to O3" setup = [GEOSChemGasPhaseSetup] begin
    # Same near-cancellation construction as the OH sensitivity test above: tight
    # solver tolerance plus a 50% perturbation to put the O3 difference ~80x above
    # the solver accuracy floor (converged to <0.001% between reltol 1e-10 and 1e-11).
    u_4 = 1.5975674827473085e-7

    @unpack O3, HO2 = sys
    o1 = solve(ODEProblem(sys, [O3 => 20, HO2 => 4.0e-6], tspan), Rosenbrock23(), abstol = 1.0e-10, reltol = 1.0e-10)
    o2 = solve(ODEProblem(sys, [O3 => 20, HO2 => 4.0e-6 * 1.5], tspan), Rosenbrock23(), abstol = 1.0e-10, reltol = 1.0e-10)
    test4 = o1[O3][end] - o2[O3][end]

    @test test4 ≈ u_4 rtol = 0.01
end

@testitem "Compose GEOSChem FastJX" begin
    using GasChem, EarthSciMLBase
    using ModelingToolkit
    gc = GEOSChemGasPhase()
    fjx = FastJX(0.0)
    gf_coupled = couple(gc, fjx)
    gf = convert(System, gf_coupled, compile = false)

    eqs = string.(equations(gf))

    j_eqs = filter(eq -> contains(eq, r"^GEOSChemGasPhase₊j_"), eqs)

    wanteqs = [
        "GEOSChemGasPhase₊j_1(t) ~ FastJX₊j_O2(t)", "GEOSChemGasPhase₊j_2(t) ~ FastJX₊j_O3(t)",
        "GEOSChemGasPhase₊j_3(t) ~ FastJX₊j_O31D(t)", "GEOSChemGasPhase₊j_6(t) ~ FastJX₊j_NO(t)",
        "GEOSChemGasPhase₊j_7(t) ~ FastJX₊j_H2COa(t)", "GEOSChemGasPhase₊j_8(t) ~ FastJX₊j_H2COb(t)",
        "GEOSChemGasPhase₊j_9(t) ~ FastJX₊j_H2O2(t)", "GEOSChemGasPhase₊j_10(t) ~ FastJX₊j_CH3OOH(t)",
        "GEOSChemGasPhase₊j_11(t) ~ FastJX₊j_NO2(t)", "GEOSChemGasPhase₊j_12(t) ~ FastJX₊j_NO3a(t)",
        "GEOSChemGasPhase₊j_13(t) ~ FastJX₊j_NO3b(t)", "GEOSChemGasPhase₊j_14(t) ~ FastJX₊j_N2O5(t)",
        "GEOSChemGasPhase₊j_15(t) ~ FastJX₊j_HNO2(t)", "GEOSChemGasPhase₊j_16(t) ~ FastJX₊j_HNO3(t)",
        "GEOSChemGasPhase₊j_18(t) ~ FastJX₊j_HNO4(t)", "GEOSChemGasPhase₊j_19(t) ~ FastJX₊j_ClNO3a(t)",
        "GEOSChemGasPhase₊j_20(t) ~ FastJX₊j_ClNO3b(t)", "GEOSChemGasPhase₊j_22(t) ~ FastJX₊j_Cl2(t)",
        "GEOSChemGasPhase₊j_24(t) ~ FastJX₊j_HOCl(t)", "GEOSChemGasPhase₊j_25(t) ~ FastJX₊j_OClO(t)",
        "GEOSChemGasPhase₊j_26(t) ~ FastJX₊j_Cl2O2(t)", "GEOSChemGasPhase₊j_27(t) ~ FastJX₊j_ClO(t)",
        "GEOSChemGasPhase₊j_28(t) ~ FastJX₊j_BrO(t)", "GEOSChemGasPhase₊j_30(t) ~ FastJX₊j_BrNO3(t)",
        "GEOSChemGasPhase₊j_32(t) ~ FastJX₊j_HOBr(t)", "GEOSChemGasPhase₊j_33(t) ~ FastJX₊j_BrCl(t)",
        "GEOSChemGasPhase₊j_34(t) ~ FastJX₊j_OCS(t)", "GEOSChemGasPhase₊j_37(t) ~ FastJX₊j_CFCl3(t)",
        "GEOSChemGasPhase₊j_38(t) ~ FastJX₊j_CF2Cl2(t)", "GEOSChemGasPhase₊j_39(t) ~ FastJX₊j_F113(t)",
        "GEOSChemGasPhase₊j_40(t) ~ FastJX₊j_F114(t)", "GEOSChemGasPhase₊j_41(t) ~ FastJX₊j_F115(t)",
        "GEOSChemGasPhase₊j_42(t) ~ FastJX₊j_CCl4(t)", "GEOSChemGasPhase₊j_43(t) ~ FastJX₊j_CH3Cl(t)",
        "GEOSChemGasPhase₊j_44(t) ~ FastJX₊j_MeCCl3(t)", "GEOSChemGasPhase₊j_45(t) ~ FastJX₊j_CH2Cl2(t)",
        "GEOSChemGasPhase₊j_46(t) ~ FastJX₊j_CHF2Cl(t)", "GEOSChemGasPhase₊j_47(t) ~ FastJX₊j_F123(t)",
        "GEOSChemGasPhase₊j_48(t) ~ FastJX₊j_F141b(t)", "GEOSChemGasPhase₊j_49(t) ~ FastJX₊j_F142b(t)",
        "GEOSChemGasPhase₊j_50(t) ~ FastJX₊j_CH3Br(t)", "GEOSChemGasPhase₊j_51(t) ~ FastJX₊j_H1211(t)",
        "GEOSChemGasPhase₊j_53(t) ~ FastJX₊j_H1301(t)", "GEOSChemGasPhase₊j_54(t) ~ FastJX₊j_H2402(t)",
        "GEOSChemGasPhase₊j_55(t) ~ FastJX₊j_CH2Br2(t)", "GEOSChemGasPhase₊j_56(t) ~ FastJX₊j_CHBr3(t)",
        "GEOSChemGasPhase₊j_59(t) ~ FastJX₊j_PAN(t)", "GEOSChemGasPhase₊j_61(t) ~ FastJX₊j_ActAld(t)",
        "GEOSChemGasPhase₊j_63(t) ~ FastJX₊j_MeVKa(t)", "GEOSChemGasPhase₊j_64(t) ~ FastJX₊j_MeVKb(t)",
        "GEOSChemGasPhase₊j_65(t) ~ FastJX₊j_MeVKc(t)", "GEOSChemGasPhase₊j_66(t) ~ FastJX₊j_MeAcr(t)",
        "GEOSChemGasPhase₊j_68(t) ~ FastJX₊j_GlyAld(t)", "GEOSChemGasPhase₊j_69(t) ~ FastJX₊j_MEKeto(t)",
        "GEOSChemGasPhase₊j_70(t) ~ FastJX₊j_PrAld(t)", "GEOSChemGasPhase₊j_71(t) ~ FastJX₊j_MGlyxl(t)",
        "GEOSChemGasPhase₊j_72(t) ~ FastJX₊j_Glyxla(t)", "GEOSChemGasPhase₊j_73(t) ~ FastJX₊j_Glyxlb(t)",
        "GEOSChemGasPhase₊j_74(t) ~ FastJX₊j_Glyxlc(t)", "GEOSChemGasPhase₊j_76(t) ~ FastJX₊j_Aceta(t)",
        "GEOSChemGasPhase₊j_77(t) ~ FastJX₊j_Acetb(t)", "GEOSChemGasPhase₊j_122(t) ~ FastJX₊j_CH3I(t)",
        "GEOSChemGasPhase₊j_134(t) ~ FastJX₊j_CH3NO3(t)",
    ]
    for eq in wanteqs
        @test contains(string(j_eqs), eq)
    end

    # Photolysis-completion: 62 organic surrogate channels added on top of the
    # original 63. Hydroperoxides -> /CH3OOH/, organic & alkyl nitrates -> /CH3NO3/, plus 13
    # dedicated Cloud-J v7.3e cross-sections (ONIT1/ETNO3/.../ENOL). j_80 (ETP) keeps the 0.5
    # j-factor and is checked separately below.
    new_parent = Dict(
        78 => "CH3NO3", 79 => "CH3OOH", 81 => "CH3OOH", 82 => "CH3OOH", 83 => "CH3OOH",
        84 => "CH3OOH", 85 => "CH3OOH", 86 => "HMHP", 87 => "CH3OOH", 88 => "MGlyxl",
        89 => "CH3NO3", 90 => "MGlyxl", 91 => "PrAld", 92 => "CH3OOH", 93 => "CH3OOH",
        94 => "ENOL", 95 => "CH3OOH", 96 => "CH3OOH", 97 => "CH3OOH", 98 => "CH3NO3",
        99 => "CH3OOH", 105 => "H2O2", 106 => "ICN", 107 => "ETHLN", 108 => "MVKN",
        109 => "MACRN", 110 => "MACRNP", 111 => "ONIT1", 112 => "ONIT1", 113 => "ONIT1",
        135 => "ETNO3", 136 => "IPRNO3", 137 => "NPRNO3", 138 => "CH3OOH", 139 => "CH3OOH",
        140 => "CH3OOH", 141 => "CH3OOH", 142 => "CH3OOH", 143 => "CH3OOH", 144 => "CH3OOH",
        145 => "CH3OOH", 146 => "ONIT1", 147 => "ONIT1", 148 => "ONIT1", 149 => "ONIT1",
        150 => "NITP", 151 => "CH3OOH", 152 => "ONIT1", 153 => "PrAld", 154 => "CH3OOH",
        155 => "HP2", 156 => "CH3OOH", 157 => "CH3OOH", 158 => "CH3OOH", 159 => "ONIT1",
        160 => "MACRNP", 161 => "PrAld", 162 => "CH3OOH", 164 => "CH3OOH", 165 => "CH3OOH",
        166 => "CH3NO3",
    )
    for (n, par) in new_parent
        @test contains(string(j_eqs), "GEOSChemGasPhase₊j_$(n)(t) ~ FastJX₊j_$(par)(t)")
    end
    # ETP (j_80) carries GEOS-Chem's 0.5 j-factor on the CH3OOH surrogate.
    @test any(
        e -> contains(e, "GEOSChemGasPhase₊j_80(t) ~") && contains(e, "0.5") &&
            contains(e, "FastJX₊j_CH3OOH(t)"),
        j_eqs
    )

    @test_nowarn convert(System, gf_coupled)
end

@testitem "Compose GEOSChem FastJX_interpolation_troposphere" begin
    using GasChem, EarthSciMLBase
    using ModelingToolkit

    # Photolysis connection equations produced when coupling GEOSChemGasPhase to a
    # given Fast-JX constructor.
    function jconnections(fjx)
        gf = convert(System, couple(GEOSChemGasPhase(), fjx), compile = false)
        sort(filter(eq -> contains(eq, r"^GEOSChemGasPhase₊j_"), string.(equations(gf))))
    end

    online = jconnections(FastJX(0.0))
    interp = jconnections(FastJX_interpolation_troposphere(0.0))   # mech=:all (full set)

    # With `mech=:all` the interpolated Fast-JX must wire the *same* GEOS-Chem
    # photolysis connections as the online `FastJX`; the `mech=:superfast` subset
    # (the historical default) only exposes the SuperFast species and cannot.
    @test !isempty(interp)
    @test interp == online

    # And the fully-coupled, compiled system must build.
    @test_nowarn convert(System, couple(GEOSChemGasPhase(), FastJX_interpolation_troposphere(0.0)))
end
