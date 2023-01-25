module GasChem

using ModelingToolkit
using Catalyst
using DifferentialEquations
using Symbolics
using Dates
using StaticArrays
using Unitful
include("SuperFast.jl")
include("Fast-JX.jl")
include("compose_fastjx_superfast.jl")

end
