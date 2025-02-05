using GasChem
using DifferentialEquations, ModelingToolkit, DynamicQuantities

tspan = (0.0, 360.0)

@testset "Base case" begin
    answer = 23.400804629407062

    rs = structural_simplify(SuperFast())
    sol = solve(
        ODEProblem(rs, [], tspan, []),
        Tsit5(),
        saveat=10.0,
        abstol=1e-12,
        reltol=1e-12,
    )

    @test sol[rs.O3][end] ≈ answer
end

@testset "ISOP sensitivity" begin
    u_isop = 1.3266279542331567

    rs1 = structural_simplify(SuperFast())
    o1 = solve(
        ODEProblem(rs1, [rs1.ISOP => 0.54], tspan, []),
        Tsit5(),
        saveat=10.0,
        abstol=1e-12,
        reltol=1e-12,
    )
    rs2 = structural_simplify(SuperFast())
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
    u_no2 = 28.430159226086346

    rs1 = structural_simplify(SuperFast())
    o1 = solve(
        ODEProblem(
            rs1,
            [rs1.NO2 => 100.0],
            tspan,
            [],
            combinatoric_ratelaws=false,
        ),
        Tsit5(),
        saveat=10.0,
        abstol=1e-12,
        reltol=1e-12,
    )
    rs2 = structural_simplify(SuperFast())
    o2 = solve(
        ODEProblem(rs2, [], tspan, []),
        Tsit5(),
        saveat=10.0,
        abstol=1e-12,
        reltol=1e-12,
    )
    test3 = o1[rs1.O3][end] - o2[rs2.O3][end]

    @test test3 ≈ u_no2
end

@testset "CO sensitivity" begin
    u_co = -3.2610697264494846

    rs1 = structural_simplify(SuperFast())
    o1 = solve(
        ODEProblem(rs1, [rs1.CO => 50.0], tspan, []),
        Tsit5(),
        saveat=10.0,
        abstol=1e-12,
        reltol=1e-12,
    )
    rs2 = structural_simplify(SuperFast())
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
    u_ch4 = 0.08634138773470212

    rs1 = structural_simplify(SuperFast())
    o1 = solve(
        ODEProblem(rs1, [rs1.CH4 => 1900.0], tspan, []),
        Tsit5(),
        saveat=10.0,
        abstol=1e-12,
        reltol=1e-12,
    )
    rs2 = structural_simplify(SuperFast())
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