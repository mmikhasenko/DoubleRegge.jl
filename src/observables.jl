"""
    phi_asymmetry(f; start = -π / 2)

Compute left right asymmetry of angular distribution.
The asymmetry is sensitive to M=2 component via interference of odd and even waves.

For example, with M=1,2 waves, the asymmetry is given by expression
```math
asymmetry = 8 / 3π * Re(w1ˣw2) / (|w1|^2 + |w2|^2)
```
where w1, and w2 are weights in front of sin(ϕ)/√π, and sin(2ϕ)/√π, respectively.
"""
function phi_asymmetry(f; start = -π / 2)
    L = quadgk(f, start + π, start + 2π)[1]
    R = quadgk(f, start, start + π)[1]
    @show L, R
    return (L - R) / (L + R)
end

function phi_asymmetry_2d(f2d; start = -π / 2)
    L = cuhre((x, f) -> f[1] = f2d(2x[1] - 1, start + π * x[2]), 2, 1)[1][1]
    R = cuhre((x, f) -> f[1] = f2d(2x[1] - 1, start + π + π * x[2]), 2, 1)[1][1]
    return (L - R) / (L + R)
end


function intensity(mass_PWs::Vector{TwoBodyPartialWaveIϕs{N, V}} where {N, V}, i::Integer)
    (((mass_PWs) .. :PWs) .. i) .. :I
end
function phase(mass_PWs::Vector{TwoBodyPartialWaveIϕs{N, V}} where {N, V}, i::Integer)
    phases = alignperiodicsequence(((mass_PWs .. :PWs) .. i) .. :ϕ)
    phases_adj = meanshiftbyperiod(phases)
    if mean(phases_adj) < π / 2
        phases_adj .+= 2π
    end
    return phases_adj
end
