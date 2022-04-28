using GasChem
using Test, SafeTestsets

@testset "GasChem.jl" begin
    @safetestset "SuperFast" begin include("superfast_test.jl") end
end
