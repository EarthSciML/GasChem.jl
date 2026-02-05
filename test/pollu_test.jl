@testitem "Pollu Structural Verification" tags = [:pollu] begin
    using ModelingToolkit
    rs = mtkcompile(Pollu())
    @test length(unknowns(rs)) == 20
end

@testitem "Pollu Reference Solution at t=60min" tags = [:pollu] begin
    using OrdinaryDiffEqRosenbrock, ModelingToolkit

    # Reference solution at t=60 minutes (3600 seconds) from Verwer (1994),
    # computed with RADAU5 at high precision, converted from ppm to ppb.
    rs = mtkcompile(Pollu())
    ic = [
        rs.NO2 => 0.0,         # y[1] = 0 ppm
        rs.NO => 0.2 * 1e3,    # y[2] = 0.2 ppm = 200 ppb
        rs.O3P => 0.0,         # y[3] = 0 ppm
        rs.O3 => 0.04 * 1e3,   # y[4] = 0.04 ppm = 40 ppb
        rs.HO2 => 0.0,         # y[5] = 0 ppm
        rs.OH => 0.0,          # y[6] = 0 ppm
        rs.CH2O => 0.1 * 1e3,  # y[7] = 0.1 ppm = 100 ppb
        rs.CO => 0.3 * 1e3,    # y[8] = 0.3 ppm = 300 ppb
        rs.ALD => 0.01 * 1e3,  # y[9] = 0.01 ppm = 10 ppb
        rs.MEO2 => 0.0,        # y[10] = 0 ppm
        rs.C2O3 => 0.0,        # y[11] = 0 ppm
        rs.CO2 => 0.0,         # y[12] = 0 ppm
        rs.PAN => 0.0,         # y[13] = 0 ppm
        rs.CH3O => 0.0,        # y[14] = 0 ppm
        rs.HNO3 => 0.0,        # y[15] = 0 ppm
        rs.O1D => 0.0,         # y[16] = 0 ppm
        rs.SO2 => 0.007 * 1e3, # y[17] = 0.007 ppm = 7 ppb
        rs.SO4 => 0.0,         # y[18] = 0 ppm
        rs.NO3 => 0.0,         # y[19] = 0 ppm
        rs.N2O5 => 0.0,        # y[20] = 0 ppm
    ]
    tspan = (0.0, 3600.0)  # 60 minutes in seconds

    sol = solve(
        ODEProblem(rs, ic, tspan),
        Rosenbrock23(),
        saveat = 10.0,
        abstol = 1e-12,
        reltol = 1e-12,
    )

    # Reference values from Verwer (1994) Fortran benchmark, converted ppm → ppb
    @test sol[rs.NO2][end] ≈ 0.5646255480022769e-01 * 1e3 rtol = 1e-3  # y[1]
    @test sol[rs.NO][end] ≈ 0.1342484130422339e+00 * 1e3 rtol = 1e-3   # y[2]
    @test sol[rs.O3][end] ≈ 0.5523140207484359e-02 * 1e3 rtol = 1e-3   # y[4]
    @test sol[rs.CH2O][end] ≈ 0.7784249118997964e-01 * 1e3 rtol = 1e-3 # y[7]
    @test sol[rs.CO][end] ≈ 0.3245075353396018e+00 * 1e3 rtol = 1e-3   # y[8]
    @test sol[rs.ALD][end] ≈ 0.7494013383880406e-02 * 1e3 rtol = 1e-3  # y[9]
    @test sol[rs.CO2][end] ≈ 0.2230505975721359e-02 * 1e3 rtol = 1e-3  # y[12]
    @test sol[rs.PAN][end] ≈ 0.2087162882798630e-03 * 1e3 rtol = 1e-3  # y[13]
    @test sol[rs.HNO3][end] ≈ 0.8964884856898295e-02 * 1e3 rtol = 1e-3 # y[15]
    @test sol[rs.SO2][end] ≈ 0.6899219696263405e-02 * 1e3 rtol = 1e-3  # y[17]
    @test sol[rs.SO4][end] ≈ 0.1007803037365946e-03 * 1e3 rtol = 1e-3  # y[18]
end

@testitem "Pollu Short Integration" tags = [:pollu] begin
    using OrdinaryDiffEqRosenbrock, ModelingToolkit
    rs = mtkcompile(Pollu())
    ic = [
        rs.NO2 => 0.0, rs.NO => 200.0, rs.O3P => 0.0, rs.O3 => 40.0,
        rs.HO2 => 0.0, rs.OH => 0.0, rs.CH2O => 100.0, rs.CO => 300.0,
        rs.ALD => 10.0, rs.MEO2 => 0.0, rs.C2O3 => 0.0, rs.CO2 => 0.0,
        rs.PAN => 0.0, rs.CH3O => 0.0, rs.HNO3 => 0.0, rs.O1D => 0.0,
        rs.SO2 => 7.0, rs.SO4 => 0.0, rs.NO3 => 0.0, rs.N2O5 => 0.0,
    ]
    tspan = (0.0, 60.0)  # 1 minute in seconds

    sol = solve(
        ODEProblem(rs, ic, tspan),
        Rosenbrock23(),
        saveat = 10.0,
        abstol = 1e-12,
        reltol = 1e-12,
    )

    # Basic qualitative checks
    @test sol[rs.O3][end] < sol[rs.O3][1]  # O3 should decrease
    @test sol[rs.CO][end] > sol[rs.CO][1]  # CO should increase
end

@testitem "Pollu Conservation Properties" tags = [:pollu] begin
    using OrdinaryDiffEqRosenbrock, ModelingToolkit
    rs = mtkcompile(Pollu())
    ic = [
        rs.NO2 => 0.0, rs.NO => 200.0, rs.O3P => 0.0, rs.O3 => 40.0,
        rs.HO2 => 0.0, rs.OH => 0.0, rs.CH2O => 100.0, rs.CO => 300.0,
        rs.ALD => 10.0, rs.MEO2 => 0.0, rs.C2O3 => 0.0, rs.CO2 => 0.0,
        rs.PAN => 0.0, rs.CH3O => 0.0, rs.HNO3 => 0.0, rs.O1D => 0.0,
        rs.SO2 => 7.0, rs.SO4 => 0.0, rs.NO3 => 0.0, rs.N2O5 => 0.0,
    ]
    tspan = (0.0, 3600.0)

    sol = solve(
        ODEProblem(rs, ic, tspan),
        Rosenbrock23(),
        saveat = 10.0,
        abstol = 1e-12,
        reltol = 1e-12,
    )

    # Sulfur conservation: SO2 + SO4 should be constant
    sulfur_init = sol[rs.SO2][1] + sol[rs.SO4][1]
    sulfur_final = sol[rs.SO2][end] + sol[rs.SO4][end]
    @test sulfur_init ≈ sulfur_final rtol = 1e-6
end

@testitem "Couple Pollu and FastJX" tags = [:pollu] begin
    using EarthSciMLBase
    using OrdinaryDiffEqRosenbrock
    using ModelingToolkit

    sf1 = couple(Pollu(), FastJX_interpolation_troposphere(0.0))
    sf2 = couple(Pollu(), FastJX(0.0))
    sys1 = convert(System, sf1)
    sys2 = convert(System, sf2)
    tspan = (0.0, 3600 * 24)
    prob1 = ODEProblem(sys1, [], tspan)
    prob2 = ODEProblem(sys2, [], tspan)
    sol1 = solve(prob1, Rosenbrock23(), saveat = 10.0)
    sol2 = solve(prob2, Rosenbrock23(), saveat = 10.0)
    @test sol1[sys1.Pollu₊O3][4320] ≈ sol2[sys2.Pollu₊O3][4320] rtol = 1e-4
end

@testitem "Couple Pollu and NEI2016MonthlyEmis" tags = [:pollu] begin
    using EarthSciMLBase, EarthSciData
    using Dates, ModelingToolkit

    domain = DomainInfo(
        DateTime(2016, 5, 1),
        DateTime(2016, 5, 4);
        lonrange = deg2rad(-115):deg2rad(2.5):deg2rad(-68.75),
        latrange = deg2rad(25):deg2rad(2):deg2rad(53.7),
        levrange = 1:15,
    )

    model_3way = couple(
        FastJX(get_tref(domain)),
        Pollu(),
        NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain),
    )

    sys = convert(System, model_3way)
    @test length(unknowns(sys)) ≈ 20

    eqs = string(equations(sys))
    wanteq = "Differential(t)(Pollu₊ALD(t)) ~ Pollu₊NEI2016MonthlyEmis_ALD2(t)"
    @test contains(string(eqs), wanteq)
end
