@testsnippet GEOSChemGasPhaseSetup begin
    using EarthSciMLBase
    using OrdinaryDiffEqRosenbrock
    using ModelingToolkit

    tspan = (0.0, 360.0)
    sys = GEOSChemGasPhase()
    sys = mtkcompile(sys)

    # Helper to replace removed ModelingToolkit.get_defaults in MTK v11
    function _get_defaults(sys)
        d = Dict{Any,Any}()
        for v in unknowns(sys)
            if ModelingToolkit.hasdefault(v)
                d[v] = ModelingToolkit.getdefault(v)
            end
        end
        for p in parameters(sys)
            if ModelingToolkit.hasdefault(p)
                d[p] = ModelingToolkit.getdefault(p)
            end
        end
        return d
    end
end

# Unit Test 0: Base case
@testitem "Base case" setup=[GEOSChemGasPhaseSetup] begin
    u_0 = [19.995176711847932, 0.002380578760821661, 0.0, 9.428652108253236e-251, 
    5.879149823894928e-6, 0.0, 1.840000000235526e7]

    vals = _get_defaults(sys)
    for k in setdiff(unknowns(sys), keys(vals))
        vals[k] = 0 # Set variables with no default to zero.
    end
    prob = ODEProblem(sys, vals, tspan)
    sol = solve(prob, Rosenbrock23())
    test0 = [sol[v][end]
             for v in [sys.O3, sys.NO2, sys.ISOP, sys.O1D, sys.OH, sys.DMS, sys.H2O]]

    test0 ≈ u_0
end

# Unit Test 1: O1D sensitivity to O3
@testitem "O1D sensitivity to O3" setup=[GEOSChemGasPhaseSetup] begin
    u_1 = 8.160593694128693e-7

    vals = _get_defaults(sys)
    @unpack O3, O1D = sys
    vals[O3] = 20
    vals[O1D] = 1E-6
    o1 = solve(ODEProblem(sys, vals, tspan), Rosenbrock23())
    vals[O1D] = 1E-6 * 1.1E0
    o2 = solve(ODEProblem(sys, vals, tspan), Rosenbrock23())
    test1 = o1[O3][end] - o2[O3][end]

    @test test1≈u_1 rtol=0.001
end

# Unit Test 2: OH sensitivity to O3
@testitem "OH sensitivity to O3" setup=[GEOSChemGasPhaseSetup] begin
    u_2 = 8.209088520061414e-9

    vals = _get_defaults(sys)
    @unpack O3, OH = sys
    vals[O3] = 20
    vals[OH] = 4E-6
    o1 = solve(ODEProblem(sys, vals, tspan), Rosenbrock23(), abstol = 1e-6, reltol = 1e-6)
    vals[OH] = 4E-6 * 1.05E0
    o2 = solve(ODEProblem(sys, vals, tspan), Rosenbrock23(), abstol = 1e-6, reltol = 1e-6)
    test2 = o1[O3][end] - o2[O3][end]

    @test test2≈u_2 rtol=0.001
end

# Unit Test 3: NO2 sensitivity to O3
@testitem "NO2 sensitivity to O3" setup=[GEOSChemGasPhaseSetup] begin
    u_3 =  1.1784680253867919e-7

    vals = _get_defaults(sys)
    @unpack O3, NO2 = sys
    vals[O3] = 20
    vals[NO2] = 4E-4
    o1 = solve(ODEProblem(sys, vals, tspan), Rosenbrock23(), abstol = 1e-6, reltol = 1e-6)
    vals[NO2] = 4E-4 * 1.05E0
    o2 = solve(ODEProblem(sys, vals, tspan), Rosenbrock23(), abstol = 1e-6, reltol = 1e-6)
    test3 = o1[O3][end] - o2[O3][end]

    @test test3≈u_3 rtol=0.01
end

# Unit Test 4: HO2 sensitivity to O3
@testitem "HO2 sensitivity to O3" setup=[GEOSChemGasPhaseSetup] begin
    u_4 = 1.6049252593575147e-8

    vals = _get_defaults(sys)
    @unpack O3, HO2 = sys
    vals[O3] = 20
    vals[HO2] = 4E-6
    o1 = solve(ODEProblem(sys, vals, tspan), Rosenbrock23(), abstol = 1e-6, reltol = 1e-6)
    vals[HO2] = 4E-6 * 1.05E0
    o2 = solve(ODEProblem(sys, vals, tspan), Rosenbrock23(), abstol = 1e-6, reltol = 1e-6)
    test4 = o1[O3][end] - o2[O3][end]

    @test test4≈u_4 rtol=0.001
end

@testitem "Compose GEOSChem FastJX" begin
    using EarthSciMLBase
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
        "GEOSChemGasPhase₊j_134(t) ~ FastJX₊j_CH3NO3(t)"
    ]
    for eq in wanteqs
        @test contains(string(j_eqs), eq)
    end

    @test_nowarn convert(System, gf_coupled)
end
