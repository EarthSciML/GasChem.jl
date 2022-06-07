using GasChem
using Test, SafeTestsets

@testset "GasChem.jl" begin
    @safetestset "SuperFast" begin include("superfast_test.jl") end
    @safetestset "FastJX" begin include("fastjx_test.jl") end
end
