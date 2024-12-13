using DrWatson
@quickactivate "DoubleRegge"

using DoubleRegge
setsystem!(:compass_ηπ)

using UnROOT
using StaticArrays
using Plots
using Parameters
using DelimitedFiles
using LinearAlgebra

theme(:wong)

vof = let
    file_name = datadir("exp_raw", "dimas_toys", "pi03pic_mc_pwa_model_2.root")
    !(isfile(file_name)) && error("no file $file_name")
    # 
    file = ROOTFile(file_name)
    # 
    tree_name = "pi03pic"
    tree = LazyTree(file, tree_name)
    #
    branchnames = ["pBeam", "pRecoil", "pPim", "pEta"]
    components = ['X', 'Y', 'Z', 'E']
    varnames = (:pb, :pr, :pπ, :pη)
    # 
    n = size(tree, 1)
    v = map(1:n) do i
        row = tree[i, :]
        map(branchnames) do p
            SVector([getproperty(row, Symbol(p * c)) for c in components]...)
        end
    end
    NamedTuple{varnames}.(v)
end;

invariants(pb, pr, pπ, pη) =
    (s0 = invmasssq(pr + pπ + pη),
        s1 = invmasssq(pη + pπ),
        s2 = invmasssq(pπ + pr),
        t1 = invmasssq(pb - pη),
        t2 = invmasssq(pb - pη - pπ))

voi = [invariants(v.pb, v.pr, v.pπ, v.pη) for v in vof];

stephist(ϕTY.(vof))
stephist(cosθ1.(voi))
histogram2d(sqrt.(getproperty.(voi, :s1)), cosθ1.(voi))

vok = [(ϕ = ϕ, cosθ = cosθ, mηπ = mηπ) for
       (ϕ, cosθ, mηπ) in zip(ϕTY.(vof), cosθ1.(voi), sqrt.(getproperty.(voi, :s1)))];

function selected_distr(vars; mηπ_bin_range = error("range"))
    sample = filter(v -> mηπ_bin_range[1] < sqrt(v.s1) < mηπ_bin_range[2], vars)
    return cosθ1.(sample)
end

const Δx = 0.04;
const x0 = 0.74;
bin_range(i) = (x0 + (i - 1) * Δx, x0 + i * Δx)

const main = readdlm(datadir("exp_pro", "main_point.txt"));
main

let
    ps = [
        let
            norm = sum(main[i, :]) / 100 * (1 - (-1))
            plot(range(-1, 1, length = 100), main[i, :] / norm, lab = (i != 56 ? "" : "corr. PLB PWA rec."), lw = 3) # , α=0.5
            stephist!(selected_distr(voi; mηπ_bin_range = bin_range(i)), bins = range(-1, 1, length = 100), st = :stephist, norm = true,
                lab = (i != 56 ? "" : "Dima's MC"), c = :black, title = "bin: $(round(x0 + (2i-1)*Δx/2, digits=2))")
        end for i in 41:56
    ]
    plot(ps..., size = (2000, 1700))
end
savefig(plotsdir("dimas_toys_etapi_costheta.pdf"))

# let bin = 15
#     s = filter(x->(x0+(bin-1)*Δx < x.mηπ < x0+bin*Δx) && (x.cosθ<0), vok)
#     stephist(getproperty.(s,:ϕ), bins=50, norm=true, lab="Dima's old", xlab="phi_TY", c=:black)
#     plot!(ϕ->((sin(ϕ)+0.0*sin(2ϕ))^2)/π, -π, π, lab="M=0 is set to zero", leg=:bottom)
#
#     dϕ = ϕ_data(bin)
#     plot!( dϕ[:,1], dϕ[:,2] ./ sum( dϕ[:,2]) ./ (2π/length(dϕ[:,2])), lab = "M=2 as in the paper")
# end
# savefig(joinpath("plots","dimas_toys_etapi_phi.pdf"))
