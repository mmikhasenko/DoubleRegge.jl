
recamp(cosθ,ϕ,amps) =
    sum(a*Psi(L,M,cosθ,ϕ) for
    (a, (L, M)) in zip(amps,LMs))

"""
    constructamps(intensities, phases)

Important: the phases are expected in degrees.
"""
constructamps(intensities, phases) = [@. sqrt(iv) * cis(ϕv/180*π)
    for (iv, ϕv) in zip(intensities, phases)];
#
#
function dNdcosθ(cosθ; amps=amps)
    list_of_all = collect(zip(amps,LMs))
    # M=1
    l1 = filter(x->x[2].M==1, list_of_all)
    l2 = filter(x->x[2].M==2, list_of_all)
    v = abs2(sum(a*Psi(L,M,cosθ,π/2) for (a, (L, M)) in l1)) +
        abs2(sum(a*Psi(L,M,cosθ,π/4) for (a, (L, M)) in l2))
    v *= π
    return v
end
#
dNdϕ(ϕ; amps=amps, cosθlims::Tuple{Real,Real}=(-1,1)) = quadgk(cosθ->abs2(recamp(cosθ,ϕ,amps)), cosθlims...)[1]
