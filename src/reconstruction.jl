
recamp(cosθ,ϕ,amps,LMs) =
    sum(a*Psi(L,M,cosθ,ϕ) for
    (a, (L, M)) in zip(amps,LMs))
#
function dNdcosθ(cosθ; amps, LMs)
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
dNdϕ(ϕ; amps, cosθlims::Tuple{Real,Real}=(-1,1), LMs) = quadgk(cosθ->abs2(recamp(cosθ,ϕ,amps,LMs)), cosθlims...)[1]

integrate_dcosθdϕ(g) = cuhre((x,f)->f[1]=g(2*x[1]-1, π*(2*x[2]-1)),2,1)[1][1]*(4π)
