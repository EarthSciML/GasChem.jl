using GasChem
using EarthSciMLBase, EarthSciData
using OrdinaryDiffEqRosenbrock, ModelingToolkit, DynamicQuantities
using Dates

tspan1 = (0.0, 60.0)
tspan2 = (0.0, 3600.0)

rs = structural_simplify(Pollu())
initial_conditions_test = [rs.NO2 => 0, 
            rs.NO => 0.2*1e3, 
            rs.O3P => 0, 
            rs.O3 => 0.04*1e3, 
            rs.HO2 => 0, 
            rs.OH => 0, 
            rs.CH2O => 0.1*1e3, 
            rs.CO => 0.3*1e3, 
            rs.ALD => 0.01*1e3, 
            rs.MEO2 => 0, 
            rs.C2O3 => 0, 
            rs.CO2 => 0, 
            rs.PAN => 0, 
            rs.CH3O => 0, 
            rs.HNO3 => 0, 
            rs.O1D => 0, 
            rs.SO2 => 0.007*1e3, 
            rs.SO4 => 0, 
            rs.NO3 => 0, 
            rs.N2O5 => 0] # Convert to ppb

@testset "Case1 - 60s" begin
    sol = solve(
        ODEProblem(rs, initial_conditions_test, tspan1, []),
        Rosenbrock23(),
        saveat = 10.0,
        abstol = 1e-12,
        reltol = 1e-12)
    @test sol[rs.O3][end] ≈ 0.32994065756620e-2*1e3 rtol = 0.001 # Convert to ppb
end

@testset "Case2 - 3600s" begin
    sol = solve(
        ODEProblem(rs, initial_conditions_test, tspan2, []),
        Rosenbrock23(),
        saveat = 10.0,
        abstol = 1e-12,
        reltol = 1e-12)
    @test sol[rs.O3][end] ≈ 0.552314020747798e-2*1e3 rtol = 0.001 # Convert to ppb
end

@testset "Couple Pollu and FastJX" begin
    sol_middle = 39.99906107410295

    sf1 = couple(Pollu(), FastJX_interpolation_troposphere(0.0))
    sf2 = couple(Pollu(), FastJX(0.0))
    sys1 = convert(ODESystem, sf1)
    sys2 = convert(ODESystem, sf2)
    tspan = (0.0, 3600 * 24)
    prob1 = ODEProblem(sys1, [], tspan, [])
    prob2 = ODEProblem(sys2, [], tspan, [])
    sol1 = solve(prob1, Rosenbrock23(), saveat = 10.0)
    sol2 = solve(prob2, Rosenbrock23(), saveat = 10.0)
    @test sol1[sys1.Pollu₊O3][4320] ≈ sol_middle rtol=1e-4
    @test sol2[sys2.Pollu₊O3][4320] ≈ sol_middle rtol=1e-4
end

@testset "Couple Pollu and NEI2016MonthlyEmis" begin
    domain = DomainInfo(
        DateTime(2016, 5, 1),
        DateTime(2016, 5, 4);
        lonrange = deg2rad(-115):deg2rad(2.5):deg2rad(-68.75),
        latrange = deg2rad(25):deg2rad(2):deg2rad(53.7),
        levrange = 1:15
    )

    model_3way = couple(
        FastJX(get_tref(domain)),
        Pollu(),
        NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain)
    )

    sys = convert(ODESystem, model_3way)
    @test length(unknowns(sys)) ≈ 20

    eqs = string(equations(sys))

    wanteq = "Differential(t)(Pollu₊ALD(t)) ~ Pollu₊NEI2016MonthlyEmis_ALD2(t)"
    @test contains(string(eqs), wanteq)
end