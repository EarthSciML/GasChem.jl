using GasChem, EarthSciData
using Test, Dates, ModelingToolkit, DifferentialEquations, EarthSciMLBase, Unitful

@testset "NEI2016Extension3way"
    @parameters t [unit = u"s"]
    model_3way = couple(FastJX(t), SuperFast(t), emis)

    sys = structural_simplify(get_mtk(sys))
    @test length(states(sys)) ≈ 18
end