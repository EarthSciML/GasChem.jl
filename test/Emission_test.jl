using GasChem
using Test, EarthSciData, Dates, ModelingToolkit, OrdinaryDiffEq, DifferentialEquations

@testset "emission" begin
    @parameters t
    ModelingToolkit.check_units(eqs...) = nothing
    composed_ode = SuperFast(t) + FastJX(t) + Emission(t)
    sys = structural_simplify(get_mtk(composed_ode))
    @test length(states(sys)) ≈ 18
end