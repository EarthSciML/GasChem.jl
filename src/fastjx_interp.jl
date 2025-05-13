function interp2_func(x1, x2, T1, T2)
    if x1 ≈ x2
        return (T) -> x1
    end
    function interp2(T)
        T = clamp(T, T1, T2)
        x1 + (T - T1) * (x2 - x1) / (T2 - T1)
    end
end

function interp3_func(x1, x2, x3, T1, T2, T3)
    if x1 ≈ x2 && x2 ≈ x3
        return (T) -> x1
    end
    function interp3(T)
        T = clamp(T, T1, T3)
        ifelse(
            T < T2,
            x1 + (T - T1) * (x2 - x1) / (T2 - T1),
            x2 + (T - T2) * (x3 - x2) / (T3 - T2)
        )
    end
end

function interp_func(temperatures, cross_sections)
    if length(temperatures) == 2
        return interp2_func(
            cross_sections[1],
            cross_sections[2],
            temperatures[1],
            temperatures[2]
        )
    elseif length(temperatures) == 3
        return interp3_func(
            cross_sections[1],
            cross_sections[2],
            cross_sections[3],
            temperatures[1],
            temperatures[2],
            temperatures[3]
        )
    end
    LinearInterpolation(temperatures, cross_sections, extrapolation_bc = Flat())
end

"""
Create a vector of interpolators to interpolate the cross sections σ (TODO: What are the units?) for different wavelengths (in nm) and temperatures (in K).

We use use linear interpolation with flat extrapolation.
"""
function create_fjx_interp(
        temperatures::Vector{Float32},
        cross_sections::Vector{SVector{18, Float32}}
)
    [interp_func(temperatures, [x[i] for x in cross_sections])
     for
     i in 1:length(cross_sections[1])]
    #[linear_interpolation(temperatures, [x[i] for x ∈ cross_sections], extrapolation_bc=Flat()) for i ∈ 1:length(WL)]
end
