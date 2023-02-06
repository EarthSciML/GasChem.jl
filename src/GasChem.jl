module GasChem

using EarthSciMLBase
using ModelingToolkit
using Catalyst
using Dates
using Unitful
include("SuperFast.jl")
include("Fast-JX.jl")
include("compose_fastjx_superfast.jl")

end
