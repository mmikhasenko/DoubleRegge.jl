using DrWatson
@quickactivate "DoubleRegge"

using DoubleRegge
using QuadGK
using Plots
using Test

theme(:wong; size=(500,350))
#
# fit
const labs=["a2/ℙ",  "a2/f2", "f2/ℙ", "f2/f2", "ℙ/ℙ", "ℙ/f2"]
const exchanges = [(α_a2, α_ℙ, true),  (α_a2, α_f2, true),
                   (α_f2, α_ℙ, false), (α_f2, α_f2, false),
                   (α_ℙ,  α_ℙ, false), (α_ℙ,  α_f2, false)];
const Np = length(exchanges)

function model(mηπ,cosθ,ϕ; pars)
    vars = (s = DoubleRegge.s0, s1 = mηπ^2,
        cosθ = cosθ, ϕ = ϕ, t2 = -0.2)
    return sum(p*modelDR(t[1], t[2], vars; η_forward=t[3])
        for (p,t) in zip(pars, exchanges))
end
# 
phi_asymmetry_model(mηπ, cosθ; pars) = phi_asymmetry(ϕ->abs2(model(mηπ, cosθ, ϕ; pars=pars)))

let
    plot()
    
    for i in 1:6
        plot!(cosθ->phi_asymmetry_model(2.8, cosθ; pars=[k==i for k in 1:Np]), -1:0.03:1,
            seriescolor=div(i+1,2), lab = labs[i],
            ls=(isodd(i) ? :solid : :dash))
    end
    plot!(xlab="cosθ", ylab="Assymetry in position of ϕ peak")
end

phi_asymmetry_2d_model(mηπ; pars) = phi_asymmetry_2d((cosθ,ϕ)->abs2(model(mηπ, cosθ, ϕ; pars=pars)))
phi_asymmetries = [phi_asymmetry_2d_model(2.6; pars=[k==i for k in 1:Np])
    for i in 1:6]
#
[println(l,": ",round(Σϕ,digits=3)) for (l, Σϕ) in zip(labs, phi_asymmetries)]
let
    plot()
    labs=["a2/ℙ",  "a2/f2", "f2/ℙ", "f2/f2", "ℙ/ℙ", "ℙ/f2"]
    for i in 1:6
        plot!(mηπ->phi_asymmetry_2d_model(mηπ; pars=[k==i for k in 1:Np]), 2.0:0.1:3.0,
            seriescolor=div(i+1,2), lab = labs[i],
            ls=(isodd(i) ? :solid : :dash))
    end
    plot!(xlab="m(ηπ)", ylab="Assymetry in position of ϕ peak")
end