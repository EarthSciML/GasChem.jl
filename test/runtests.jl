using GasChemMTK
using Test, SafeTestsets

@testset "GasChemMTK.jl" begin
    @safetestset "SuperFast" begin include("superfast_test.jl") end
end
