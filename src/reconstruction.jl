arg(z) = atan(imag(z), real(z))
#
amplitude(I,ϕ) = sqrt(I) * cis(ϕ) 
amplitude(; I,ϕ) = amplitude(I,ϕ)
Iϕ(A; ref=1.0) = NamedTuple{(:I,:ϕ)}((abs2(A),arg(A*ref')))
# 
recamp(cosθ,ϕ,amps,LMs) =
    sum(a*Psi(L,M,cosθ,ϕ) for
    (a, (L, M)) in zip(amps,LMs))
#
recamp(cosθ,ϕ, expansion::TwoBodyPartialWaves{N,Complex{V}} where N where V) =
    recamp(cosθ, ϕ, expansion.PWs, expansion.LMs)
# 
function dNdcosθ(cosθ, expansion)
    list_of_all = collect(zip(expansion.PWs,expansion.LMs))
    # M=1
    l1 = filter(x->x[2][2]==1, list_of_all)
    l2 = filter(x->x[2][2]==2, list_of_all)
    v = abs2(sum(a*Psi(L,M,cosθ,π/2) for (a, (L, M)) in l1)) +
        abs2(sum(a*Psi(L,M,cosθ,π/4) for (a, (L, M)) in l2))
    v *= π
    return v
end
#
dNdϕ(ϕ, expansion; cosθlims::Tuple{Real,Real}=(-1,1), LMs) =
    quadgk(cosθ->abs2(recamp(cosθ,ϕ,expansion)), cosθlims...)[1]

function integrate_dcosθdϕ(g, cosθlims=(-1,1), ϕlims=(-π+0.31,π+0.31); dims::Int=1)
    integration = cuhre((x,f)->f[:] .= g(
        cosθlims[1]+x[1]*(cosθlims[2]-cosθlims[1]),
           ϕlims[1]+x[2]*(   ϕlims[2]-   ϕlims[1])),2, dims)
    return integration[1] .* ((cosθlims[2]-cosθlims[1])*(ϕlims[2]-ϕlims[1]))
end

