module GasChem

using EarthSciMLBase
using ModelingToolkit
using Catalyst
using Dates
using DynamicQuantities
using StaticArrays
using Interpolations
using ModelingToolkit:t

@register_unit ppb 1

include("SuperFast.jl")
include("geoschem_ratelaws.jl")
include("geoschem_fullchem.jl")
include("Fast-JX.jl")
include("fastjx_couplings.jl")

end
