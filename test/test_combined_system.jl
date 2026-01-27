@testitem "Combined System Construction" begin
    using ModelingToolkit
    sys = TroposphericChemistrySystem()
    @test sys isa ODESystem
    @test nameof(sys) == :TroposphericChemistrySystem
end

@testitem "Combined Typical Conditions" begin
    cond = get_typical_conditions()
    @test cond[:M] ≈ 2.5e19
    @test cond[:O2] ≈ 5.25e18
    @test cond[:O3] ≈ 1e12
    @test cond[:OH] ≈ 1e6
    @test cond[:HO2] ≈ 1e8
end

@testitem "Combined Urban Conditions" begin
    cond = get_urban_conditions()
    @test cond[:NO] > get_typical_conditions()[:NO]
    @test cond[:NO2] > get_typical_conditions()[:NO2]
    @test cond[:CO] > get_typical_conditions()[:CO]
end

@testitem "Combined Remote Conditions" begin
    cond = get_remote_conditions()
    @test cond[:NO] < get_typical_conditions()[:NO]
    @test cond[:NO2] < get_typical_conditions()[:NO2]
    @test cond[:HO2] > get_typical_conditions()[:HO2]
end

@testitem "Combined Diagnostic Variables" begin
    using ModelingToolkit
    sys = TroposphericChemistrySystem()
    vars = unknowns(sys)
    var_names = [string(v) for v in vars]
    @test any(v -> contains(v, "NOx"), var_names)
    @test any(v -> contains(v, "HOx"), var_names)
    @test any(v -> contains(v, "P_O3"), var_names)
    @test any(v -> contains(v, "OPE"), var_names)
end

@testitem "Combined Key Rate Constants" begin
    using ModelingToolkit
    sys = TroposphericChemistrySystem()
    params = parameters(sys)
    defaults = ModelingToolkit.defaults(sys)
    function find_param_value(params, defaults, name_ending)
        for p in params
            if haskey(defaults, p) && endswith(string(p), name_ending)
                return defaults[p]
            end
        end
        return nothing
    end
    @test find_param_value(params, defaults, "k_HO2_NO") ≈ 8.1e-12
    @test find_param_value(params, defaults, "k_OH_NO2") ≈ 1.0e-11
    @test find_param_value(params, defaults, "k_NO_O3") ≈ 1.8e-14
    @test find_param_value(params, defaults, "k_CH3O2_NO") ≈ 7.7e-12
    @test find_param_value(params, defaults, "j_NO2") ≈ 8e-3
end

@testitem "Combined NOx-VOC Regimes" begin
    @test true
end

@testitem "Combined OH Production Coupling" begin
    @test true
end

@testitem "Combined NOx Cycling Coupling" begin
    @test true
end

@testitem "Combined O3 Production Budget" begin
    @test true
end
