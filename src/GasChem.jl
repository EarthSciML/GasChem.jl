module GasChem

using EarthSciMLBase
using ModelingToolkit
using Catalyst
using Dates
using Unitful
using StaticArrays
using Interpolations

include("SuperFast.jl")
include("geoschem_ratelaws.jl")
include("geoschem_fullchem.jl")
include("Fast-JX.jl")
include("fastjx_couplings.jl")

end
