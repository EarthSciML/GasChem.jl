@testitem "Base case" begin
    using OrdinaryDiffEqRosenbrock, ModelingToolkit
    tspan = (0.0, 360.0)
    answer = 23.400804629407062

    rs = mtkcompile(SuperFast())
    sol = solve(
        ODEProblem(rs, [], tspan),
        Rosenbrock23(),
        saveat = 10.0,
        abstol = 1.0e-12,
        reltol = 1.0e-12
    )

    @test sol[rs.O3][end] ≈ answer
end

@testitem "ISOP sensitivity" begin
    using OrdinaryDiffEqRosenbrock, ModelingToolkit
    tspan = (0.0, 360.0)
    u_isop = 1.3266279542331567

    rs1 = mtkcompile(SuperFast())
    o1 = solve(
        ODEProblem(rs1, [rs1.ISOP => 0.54], tspan),
        Rosenbrock23(),
        saveat = 10.0,
        abstol = 1.0e-12,
        reltol = 1.0e-12
    )
    rs2 = mtkcompile(SuperFast())
    o2 = solve(
        ODEProblem(rs2, [rs2.ISOP => 0.13], tspan),
        Rosenbrock23(),
        saveat = 10.0,
        abstol = 1.0e-12,
        reltol = 1.0e-12
    )
    test2 = o1[rs1.O3][end] - o2[rs2.O3][end]

    @test test2 ≈ u_isop
end

@testitem "NO2 sensitivity" begin
    using OrdinaryDiffEqRosenbrock, ModelingToolkit
    tspan = (0.0, 360.0)
    u_no2 = 28.430159226086346

    rs1 = mtkcompile(SuperFast())
    o1 = solve(
        ODEProblem(rs1, [rs1.NO2 => 100.0], tspan),
        Rosenbrock23(),
        saveat = 10.0,
        abstol = 1.0e-12,
        reltol = 1.0e-12
    )
    rs2 = mtkcompile(SuperFast())
    o2 = solve(
        ODEProblem(rs2, [], tspan),
        Rosenbrock23(),
        saveat = 10.0,
        abstol = 1.0e-12,
        reltol = 1.0e-12
    )
    test3 = o1[rs1.O3][end] - o2[rs2.O3][end]

    @test test3 ≈ u_no2
end

@testitem "CO sensitivity" begin
    using OrdinaryDiffEqRosenbrock, ModelingToolkit
    tspan = (0.0, 360.0)
    u_co = -3.2610697264494846

    rs1 = mtkcompile(SuperFast())
    o1 = solve(
        ODEProblem(rs1, [rs1.CO => 50.0], tspan),
        Rosenbrock23(),
        saveat = 10.0,
        abstol = 1.0e-12,
        reltol = 1.0e-12
    )
    rs2 = mtkcompile(SuperFast())
    o2 = solve(
        ODEProblem(rs2, [rs2.CO => 500.0], tspan),
        Rosenbrock23(),
        saveat = 10.0,
        abstol = 1.0e-12,
        reltol = 1.0e-12
    )
    test4 = o1[rs1.O3][end] - o2[rs2.O3][end]

    @test test4 ≈ u_co
end

@testitem "CH4 sensitivity" begin
    using OrdinaryDiffEqRosenbrock, ModelingToolkit
    tspan = (0.0, 360.0)
    u_ch4 = 0.08634138773470212

    rs1 = mtkcompile(SuperFast())
    o1 = solve(
        ODEProblem(rs1, [rs1.CH4 => 1900.0], tspan),
        Rosenbrock23(),
        saveat = 10.0,
        abstol = 1.0e-12,
        reltol = 1.0e-12
    )
    rs2 = mtkcompile(SuperFast())
    o2 = solve(
        ODEProblem(rs2, [rs2.CH4 => 1600.0], tspan),
        Rosenbrock23(),
        saveat = 10.0,
        abstol = 1.0e-12,
        reltol = 1.0e-12
    )
    test5 = o1[rs1.O3][end] - o2[rs1.O3][end]

    @test test5 ≈ u_ch4
end
