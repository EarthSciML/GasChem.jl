using GasChem
using Test, SafeTestsets

@testset "GasChem.jl" begin
    # @safetestset "SuperFast" begin
    #      include("superfast_test.jl")
    # end
    # @safetestset "FastJX" begin
    #     include("fastjx_test.jl")
    # end
    # @safetestset "compose_fastjx_superfast" begin
    #     include("compose_fastjx_superfast_test.jl")
    # end
    # @safetestset "geoschem_test" begin
    #      include("geoschem_test.jl")
    # end
    @safetestset "Emission_test.jl" begin
         include("Emission_test.jl")
    end
end
