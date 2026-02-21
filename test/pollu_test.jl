@testitem "Pollu O3 Reference @60s and @3600s" tags = [:pollu] begin
    using GasChem
    using OrdinaryDiffEqRosenbrock, ModelingToolkit

    rs = mtkcompile(Pollu())

    # Initial conditions (ppm → ppb)
    ic = [
        rs.NO2 => 0.0,
        rs.NO => 0.2 * 1.0e3,
        rs.O3P => 0.0,
        rs.O3 => 0.04 * 1.0e3,
        rs.HO2 => 0.0,
        rs.OH => 0.0,
        rs.CH2O => 0.1 * 1.0e3,
        rs.CO => 0.3 * 1.0e3,
        rs.ALD => 0.01 * 1.0e3,
        rs.MEO2 => 0.0,
        rs.C2O3 => 0.0,
        rs.CO2 => 0.0,
        rs.PAN => 0.0,
        rs.CH3O => 0.0,
        rs.HNO3 => 0.0,
        rs.O1D => 0.0,
        rs.SO2 => 0.007 * 1.0e3,
        rs.SO4 => 0.0,
        rs.NO3 => 0.0,
        rs.N2O5 => 0.0,
    ]

    sol_60 = solve(
        ODEProblem(rs, ic, (0.0, 60.0)),
        Rosenbrock23();
        saveat = 10.0,
        abstol = 1.0e-12,
        reltol = 1.0e-12
    )
    @test sol_60[rs.O3][end] ≈ 0.3299406575662e-2 * 1.0e3 rtol = 1.0e-3

    sol_3600 = solve(
        ODEProblem(rs, ic, (0.0, 3600.0)),
        Rosenbrock23();
        saveat = 10.0,
        abstol = 1.0e-12,
        reltol = 1.0e-12
    )
    @test sol_3600[rs.O3][end] ≈ 0.552314020747798e-2 * 1.0e3 rtol = 1.0e-3
end

@testitem "Couple Pollu and FastJX" tags = [:pollu] begin
    using GasChem, EarthSciMLBase
    using OrdinaryDiffEqRosenbrock
    using ModelingToolkit

    sol_middle = 39.99906107410295

    sf1 = couple(Pollu(), FastJX_interpolation_troposphere(0.0))
    sf2 = couple(Pollu(), FastJX(0.0))
    sys1 = convert(System, sf1)
    sys2 = convert(System, sf2)

    tspan = (0.0, 3600 * 24)
    prob1 = ODEProblem(sys1, [], tspan; build_initializeprob = false)
    prob2 = ODEProblem(sys2, [], tspan; build_initializeprob = false)

    sol1 = solve(prob1, Rosenbrock23(); saveat = 10.0)
    sol2 = solve(prob2, Rosenbrock23(); saveat = 10.0)

    @test sol1[sys1.Pollu₊O3][4320] ≈ sol_middle rtol = 1.0e-4
    @test sol2[sys2.Pollu₊O3][4320] ≈ sol_middle rtol = 1.0e-4
end

@testitem "Couple Pollu and NEI2016MonthlyEmis" tags = [:pollu] begin
    using GasChem, EarthSciMLBase, EarthSciData
    using Dates, ModelingToolkit

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

    sys = convert(System, model_3way)
    @test length(unknowns(sys)) ≈ 20

    eqs = string(equations(sys))
    wanteq = "Differential(t)(Pollu₊ALD(t)) ~ Pollu₊NEI2016MonthlyEmis_ALD2(t)"
    @test contains(eqs, wanteq)
end
