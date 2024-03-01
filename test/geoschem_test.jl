using GasChem
using Test
using DifferentialEquations, ModelingToolkit

tspan = (0.0, 360.0)
@variables t
sys = GEOSChemGasPhase(t).sys
structural_simplify(sys)

# Unit Test 0: Base case
@testset "Base case" begin
    u_0 = [
        2.1327512182118245
        211.8938586170381
        83.2370879010155
        17.646210820369
        606.1557601209427
        5.8614918579823065
        56.63543840034316
    ]

    defaults = ModelingToolkit.get_defaults(sys)
    prob = ODEProblem(sys, defaults, tspan, defaults)
    sol = solve(prob, AutoTsit5(Rosenbrock23()))
    test0 = [sol[v][end] for v in [sys.O3, sys.NO2, sys.ISOP, sys.O1D, sys.OH, sys.DMS, sys.H2O]]

    test0 ≈ u_0
end


# Unit Test 1: O1D sensitivity to O3
@testset "O1D sensitivity to O3" begin
    u_1 = 1.7279825730298626e-5

    vals = ModelingToolkit.get_defaults(sys)
    @unpack O3, O1D = sys
    vals[O3] = 20
    vals[O1D] = 0
    o1 = solve(ODEProblem(sys, vals, tspan, vals), AutoTsit5(Rosenbrock23()))
    vals[O1D] = 10
    o2 = solve(ODEProblem(sys, vals, tspan, vals), AutoTsit5(Rosenbrock23()))
    test1 = o1[O3][end] - o2[O3][end]

    @test test1 ≈ u_1 rtol = 0.001
end

# Unit Test 2: OH sensitivity to O3
@testset "OH sensitivity to O3" begin
    u_2 = 5.209463225241961e-7

    vals = ModelingToolkit.get_defaults(sys)
    @unpack O3, OH = sys
    vals[O3] = 20
    vals[OH] = 0
    o1 = solve(ODEProblem(sys, vals, tspan, vals), AutoTsit5(Rosenbrock23()))
    vals[OH] = 1000
    o2 = solve(ODEProblem(sys, vals, tspan, vals), AutoTsit5(Rosenbrock23()))
    test2 = o1[O3][end] - o2[O3][end]

    @test test2 ≈ u_2 rtol = 0.001
end

# Unit Test 3: NO2 sensitivity to O3
@testset "NO2 sensitivity to O3" begin
    u_3 = 1.8946337831948767e-9

    vals = ModelingToolkit.get_defaults(sys)
    @unpack O3, NO2 = sys
    vals[O3] = 20
    vals[NO2] = 20
    o1 = solve(ODEProblem(sys, vals, tspan, vals), AutoTsit5(Rosenbrock23()))
    vals[NO2] = 4000
    o2 = solve(ODEProblem(sys, vals, tspan, vals), AutoTsit5(Rosenbrock23()))
    test3 = o1[O3][end] - o2[O3][end]

    @test test3 ≈ u_3 rtol = 0.01
end

# Unit Test 4: HO2 sensitivity to O3
@testset "HO2 sensitivity to O3" begin
    u_4 = -6.469934550779044e-6

    vals = ModelingToolkit.get_defaults(sys)
    @unpack O3, HO2 = sys
    vals[O3] = 20
    vals[HO2] = 0
    o1 = solve(ODEProblem(sys, vals, tspan, vals), AutoTsit5(Rosenbrock23()))
    vals[HO2] = 20
    o2 = solve(ODEProblem(sys, vals, tspan, vals), AutoTsit5(Rosenbrock23()))
    test4 = o1[O3][end] - o2[O3][end]

    @test test4 ≈ u_4 rtol = 0.001
end

@testset "Compose GEOSChem FastJX" begin
    gc = GEOSChemGasPhase(t)
    fjx = FastJX(t)
    gf = gc + fjx
    gf = get_mtk(gf)

    structural_simplify(gf)

    eqs = string(equations(gf))

    wanteqs = ["GEOSChemGasPhase₊j_9(t) ~ uconv*fastjx₊j_h2o2(t)",
        "GEOSChemGasPhase₊j_7(t) ~ uconv*fastjx₊j_CH2Oa(t)",
        "GEOSChemGasPhase₊j_10(t) ~ uconv*fastjx₊j_CH3OOH(t)",
        "GEOSChemGasPhase₊j_11(t) ~ uconv*fastjx₊j_NO2(t)",
        "GEOSChemGasPhase₊j_3(t) ~ c_fixme1*fastjx₊j_o31D(t)"]

    for eq in wanteqs
        @test contains(string(eqs), wanteqs[1])
    end
end