using DrWatson
@quickactivate "DoubleRegge"

using DoubleRegge
using QuadGK
using Plots
using Test

theme(:wong2;
    frame = :box, grid = false, lab = "",
    xlims = (:auto, :auto), ylims = (:auto, :auto))
#
#
# fit
const exchanges = six_exchanges
const labs = getproperty.(exchanges, :label)
const Np = length(exchanges)
const reaction_system = compass_ηπ
const model = DoubleReggeModel(exchanges, -0.2, 0.9, reaction_system, ones(Np))
# 
phi_asymmetry_model(mηπ, cosθ; pars) = phi_asymmetry(ϕ -> abs2(amplitude(with_parameters(model, pars), mηπ, cosθ, ϕ)))

let
    plot()

    for i in 1:6
        plot!(cosθ -> phi_asymmetry_model(2.8, cosθ; pars = [k == i for k in 1:Np]), -1:0.03:1,
            seriescolor = div(i + 1, 2), lab = labs[i],
            ls = (isodd(i) ? :solid : :dash))
    end
    plot!(xlab = "cosθ", ylab = "Assymetry in position of ϕ peak")
end

phi_asymmetry_2d_model(mηπ; pars) = phi_asymmetry_2d((cosθ, ϕ) -> abs2(amplitude(with_parameters(model, pars), mηπ, cosθ, ϕ)))
phi_asymmetries = [phi_asymmetry_2d_model(2.6; pars = [k == i for k in 1:Np])
                   for i in 1:6]
#
[println(l, ": ", round(Σϕ, digits = 3)) for (l, Σϕ) in zip(labs, phi_asymmetries)]
let
    plot()
    for i in 1:6
        plot!(mηπ -> phi_asymmetry_2d_model(mηπ; pars = [k == i for k in 1:Np]), 2.0:0.1:3.0,
            seriescolor = div(i + 1, 2), lab = labs[i],
            ls = (isodd(i) ? :solid : :dash))
    end
    plot!(xlab = "m(ηπ)", ylab = "Assymetry in position of ϕ peak")
end
