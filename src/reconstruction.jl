
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

function integrate_dcosθdϕ(g, cosθlims=(-1,1), ϕlims=(-π,π))
    integration = cuhre((x,f)->f[1]=g(
        cosθlims[1]+x[1]*(cosθlims[2]-cosθlims[1]),
           ϕlims[1]+x[2]*(   ϕlims[2]-   ϕlims[1])),2,1)
    return integration[1][1]*(cosθlims[2]-cosθlims[1])*(ϕlims[2]-ϕlims[1])
end

