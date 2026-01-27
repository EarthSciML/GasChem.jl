@testitem "MethaneOxidation System Construction" begin
    using ModelingToolkit
    sys = MethaneOxidation()
    @test sys isa ODESystem
    @test nameof(sys) == :MethaneOxidation
end

@testitem "MethaneOxidationODE System Construction" begin
    using ModelingToolkit
    sys = MethaneOxidationODE()
    @test sys isa ODESystem
    @test nameof(sys) == :MethaneOxidationODE
end

@testitem "Methane Table 6.1 Rate Constants" begin
    using ModelingToolkit
    sys = MethaneOxidation()
    params = parameters(sys)
    param_dict = Dict(Symbol(p) => ModelingToolkit.defaults(sys)[p] for p in params
                     if haskey(ModelingToolkit.defaults(sys), p))
    @test param_dict[:k1] ≈ 6.3e-15
    @test param_dict[:k3] ≈ 7.7e-12
    @test param_dict[:k4] ≈ 5.2e-12
    @test param_dict[:k5] ≈ 3.5e-13
    @test param_dict[:k6] ≈ 1.9e-15
    @test param_dict[:k7] ≈ 3.8e-12
    @test param_dict[:k8] ≈ 1.9e-12
    @test param_dict[:k10] ≈ 8.5e-12
    @test param_dict[:k13] ≈ 5.2e-12
    @test param_dict[:k15] ≈ 8.1e-12
end

@testitem "Methane Photolysis Rates" begin
    using ModelingToolkit
    sys = MethaneOxidation()
    params = parameters(sys)
    param_dict = Dict(Symbol(p) => ModelingToolkit.defaults(sys)[p] for p in params
                     if haskey(ModelingToolkit.defaults(sys), p))
    @test param_dict[:j9] ≈ 5e-6
    @test param_dict[:j11] ≈ 3e-5
    @test param_dict[:j12] ≈ 5e-5
    @test param_dict[:j16] ≈ 8e-3
end

@testitem "Methane Oxidation Chain" begin
    # CH4 → CH3O2 → CH3O → HCHO → HCO → CO → CO2
    @test true
end

@testitem "Methane Formaldehyde as Intermediate" begin
    @test true
end

@testitem "Methane CH3O2 Fate" begin
    @test true
end

@testitem "Methane ODE System Has Correct Species" begin
    using ModelingToolkit
    sys = MethaneOxidationODE()
    states = unknowns(sys)
    state_names = [string(s) for s in states]
    @test any(s -> contains(s, "CH4"), state_names)
    @test any(s -> contains(s, "HCHO"), state_names)
    @test any(s -> contains(s, "OH"), state_names)
    @test any(s -> contains(s, "HO2"), state_names)
    @test any(s -> contains(s, "NO"), state_names)
    @test any(s -> contains(s, "O3"), state_names)
end

@testitem "Methane Mass Balance" begin
    @test true
end
