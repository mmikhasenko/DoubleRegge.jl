@recipe function f(x, mass_PWs::Vector{TwoBodyPartialWaveAs{N, T}} where {N, T},
    what::Symbol, i::Integer; iref = 2)
    mass_PWsIϕ = changerepresentation.(mass_PWs; iref = iref)
    (x, mass_PWsIϕ, what, i)
end

@recipe function f(x, mass_PWs::Vector{TwoBodyPartialWaveIϕs{N, V}} where {N, V},
    what::Symbol, i::Integer)
    y = []
    label --> ""
    L, M = mass_PWs[1].LMs[i]
    title --> "LM=$L$M"
    if what == :I
        intensities = intensity(mass_PWs, i)
        y = intensities
    end
    if what == :ϕ
        phases_adj = phase(mass_PWs, i)
        y = phases_adj .* (180 / π)
    end
    (x, y)
end