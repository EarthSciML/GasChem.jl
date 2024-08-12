using GasChem
using DifferentialEquations, ModelingToolkit, Unitful

tspan = (0.0, 360.0)

@testset "Base case" begin
    answer = 18.861830827565885

    @parameters t [unit = u"s"]
    rs = structural_simplify(SuperFast(t))
    sol = solve(
        ODEProblem(rs, [], tspan, []),
        Tsit5(),
        saveat=10.0,
        abstol=1e-12,
        reltol=1e-12,
    )

    @test sol[rs.O3][end] ≈ answer
end

@testset "DMS sensitivity" begin
    u_dms = 0.8842096169345286

    @parameters t [unit = u"s"]
    rs1 = structural_simplify(SuperFast(t))
    o1 = solve(
        ODEProblem(rs1, [rs1.DMS => 76], tspan, []),
        Tsit5(),
        saveat=10.0,
        abstol=1e-12,
        reltol=1e-12,
    )
    rs2 = structural_simplify(SuperFast(t))
    o2 = solve(
        ODEProblem(rs2, [rs2.DMS => 46], tspan, []),
        Tsit5(),
        saveat=10.0,
        abstol=1e-12,
        reltol=1e-12,
    )
    test1 = o1[rs1.SO2][end] - o2[rs2.SO2][end]

    @test test1 ≈ u_dms
end

@testset "ISOP sensitivity" begin
    u_isop = 0.19386790460198

    @parameters t [unit = u"s"]
    rs1 = structural_simplify(SuperFast(t))
    o1 = solve(
        ODEProblem(rs1, [rs1.ISOP => 0.54], tspan, []),
        Tsit5(),
        saveat=10.0,
        abstol=1e-12,
        reltol=1e-12,
    )
    rs2 = structural_simplify(SuperFast(t))
    o2 = solve(
        ODEProblem(rs2, [rs2.ISOP => 0.13], tspan, []),
        Tsit5(),
        saveat=10.0,
        abstol=1e-12,
        reltol=1e-12,
    )
    test2 = o1[rs1.O3][end] - o2[rs2.O3][end]

    @test test2 ≈ u_isop
end

@testset "NO2 sensitivity" begin
    u_no2 = 45.85359224356945

    @parameters t [unit = u"s"]
    rs1 = structural_simplify(SuperFast(t))
    o1 = solve(
        ODEProblem(
            rs1,
            [rs1.NO2 => 100.0, rs1.DMS => 0.1],
            tspan,
            [],
            combinatoric_ratelaws=false,
        ),
        Tsit5(),
        saveat=10.0,
        abstol=1e-12,
        reltol=1e-12,
    )
    rs2 = structural_simplify(SuperFast(t))
    o2 = solve(
        ODEProblem(rs2, [rs2.DMS => 0.1], tspan, []),
        Tsit5(),
        saveat=10.0,
        abstol=1e-12,
        reltol=1e-12,
    )
    test3 = o1[rs1.O3][end] - o2[rs2.O3][end]

    @test test3 ≈ u_no2
end

@testset "CO sensitivity" begin
    u_co = -0.1938631791778107

    @parameters t [unit = u"s"]
    rs1 = structural_simplify(SuperFast(t))
    o1 = solve(
        ODEProblem(rs1, [rs1.CO => 50.0], tspan, []),
        Tsit5(),
        saveat=10.0,
        abstol=1e-12,
        reltol=1e-12,
    )
    rs2 = structural_simplify(SuperFast(t))
    o2 = solve(
        ODEProblem(rs2, [rs2.CO => 500.0], tspan, []),
        Tsit5(),
        saveat=10.0,
        abstol=1e-12,
        reltol=1e-12,
    )
    test4 = o1[rs1.O3][end] - o2[rs2.O3][end]

    @test test4 ≈ u_co
end

@testset "CH4 sensitivity" begin
    u_ch4 = 0.006664852234028018

    @parameters t [unit = u"s"]
    rs1 = structural_simplify(SuperFast(t))
    o1 = solve(
        ODEProblem(rs1, [rs1.CH4 => 1900.0], tspan, []),
        Tsit5(),
        saveat=10.0,
        abstol=1e-12,
        reltol=1e-12,
    )
    rs2 = structural_simplify(SuperFast(t))
    o2 = solve(
        ODEProblem(rs2, [rs2.CH4 => 1600.0], tspan, []),
        Tsit5(),
        saveat=10.0,
        abstol=1e-12,
        reltol=1e-12,
    )
    test5 = o1[rs1.O3][end] - o2[rs1.O3][end]

    @test test5 ≈ u_ch4
end
