using DrWatson
@quickactivate "DoubleRegge"

using DoubleRegge

using UpROOT
using StaticArrays
using Plots
using Parameters
using DelimitedFiles
using LinearAlgebra

theme(:wong)
##########################################

# let f = TFile(datadir("exp_raw","prod4_data","Slot4EtaPi_MC.root"))
#     keys(f)
#     t = f["EtaPi"];
#     propertynames(t)
# end
#
##########################################
function list_of_vectors(tree_entry;
        varnames=error("list of four-vector names"),
        branchnames=error("list of branch names"),
        components=['X', 'Y', 'Z', 'E'],
        between="",
        before="")
    tpl = NamedTuple{varnames}(
        [SVector([getproperty(tree_entry,Symbol(before*p*between*c)) for c in ['X', 'Y', 'Z', 'E']]...)
            for p in branchnames])
    return tpl
end
invariants(pb,pr,pπ,pη) =
    (s0 = invmasssq(pr+pπ+pη),
     s1 = invmasssq(pη+pπ),
     s2 = invmasssq(pπ+pr),
     t1 = invmasssq(pb-pη),
     t2 = invmasssq(pb-pη-pπ))

##########################################

vof = broadcast_over_tree(
    row->list_of_vectors(row;
        varnames=(:pb, :pr, :pπ, :pη),
        branchnames = ["pBeam", "pRecoil", "pPim", "pEta"],
        components=['X', 'Y', 'Z', 'E']);
    filename=datadir("exp_raw","prod4_data","Slot4EtaPi.root"),
    treename="EtaPi");

# constraint_mass

# DoubleRegge.mπ, sqrt(invmasssq(vof[2].pb))
# DoubleRegge.mπ, sqrt(invmasssq(vof[2].pπ))
# DoubleRegge.mη, sqrt(invmasssq(vof[2].pη))
# DoubleRegge.mp, sqrt(invmasssq(vof[2].pr))

vofc =
    [(pb = SVector(v.pb[1:3]...,sqrt(DoubleRegge.mπ2+sum(abs2,v.pb[1:3]))),
      pr = SVector(v.pr[1:3]...,sqrt(DoubleRegge.mp2+sum(abs2,v.pr[1:3]))),
      pπ = SVector(v.pπ[1:3]...,sqrt(DoubleRegge.mπ2+sum(abs2,v.pπ[1:3]))),
      pη = SVector(v.pη[1:3]...,sqrt(DoubleRegge.mη2+sum(abs2,v.pη[1:3])))) for v in vof];

voi = [invariants(v.pb,v.pr,v.pπ,v.pη) for v in vofc];

##########################################
# MC

vof_mc = broadcast_over_tree(
    row->list_of_vectors(row;
        varnames=(:pb, :pr, :pπ, :pη),
        branchnames = ["pBeam", "pRecoil", "pPim", "pEta"],
        components=['X', 'Y', 'Z', 'E']);
    filename=datadir("exp_raw","prod4_data","Slot4EtaPi_MC.root"),
    treename="EtaPi");

# constraint_mass

# DoubleRegge.mπ, sqrt(invmasssq(vof_mc[2].pb))
# DoubleRegge.mπ, sqrt(invmasssq(vof_mc[2].pπ))
# DoubleRegge.mη, sqrt(invmasssq(vof_mc[2].pη))
# DoubleRegge.mp, sqrt(invmasssq(vof_mc[2].pr))

vofc_mc =
    [(pb = SVector(v.pb[1:3]...,sqrt(DoubleRegge.mπ2+sum(abs2,v.pb[1:3]))),
      pr = SVector(v.pr[1:3]...,sqrt(DoubleRegge.mp2+sum(abs2,v.pr[1:3]))),
      pπ = SVector(v.pπ[1:3]...,sqrt(DoubleRegge.mπ2+sum(abs2,v.pπ[1:3]))),
      pη = SVector(v.pη[1:3]...,sqrt(DoubleRegge.mη2+sum(abs2,v.pη[1:3])))) for v in vof_mc];

voi_mc = [invariants(v.pb,v.pr,v.pπ,v.pη) for v in vofc_mc];

##########################################

function selected_distr(vars; mηπ_bin_range=error("range"))
    sample = filter(v->mηπ_bin_range[1] < sqrt(v.s1) < mηπ_bin_range[2], vars)
    return cosθ1.(sample)
end

const Δx = 0.04;
const x0 = 0.74;
bin_range(i) = (x0 + (i-1)*Δx, x0 + i*Δx)

const main = readdlm(datadir("exp_pro","main_point.txt"));

##########################################

let
    ps = [let
        norm = sum(main[i,:]) / 100 * (1-(-1))
        plot(range(-1,1,length=100), main[i,:] / norm, lab=(i!=56 ? "" : "corr. PLB PWA rec."), lw=3) # , α=0.5
        stephist!(selected_distr(voi; mηπ_bin_range=bin_range(i)), bins=range(-1,1,length=100), norm=true,
            lab=(i!=56 ? "" : "Henri's data"), c=:black, title="bin: $(round(x0 + (2i-1)*Δx/2, digits=2))")
        #
        d = selected_distr(voi_mc; mηπ_bin_range=bin_range(i))
        stephist!(d, weights=fill(1.5*30/length(d), length(d)), bins=range(-1,1,length=30),
            lab=(i!=56 ? "" : "Henri's MC"), c=:red)
    end for i in 41:56]
    plot(ps..., size=(2000,1700))
end
savefig(plotsdir("prod4_etapi_costheta.pdf"))

# let bin = 15
#     s = filter(x->(x0+(bin-1)*Δx < x.mηπ < x0+bin*Δx) && (x.cosθ<0), vok)
#     stephist(getproperty.(s,:ϕ), bins=50, norm=true, lab="Dima's old", xlab="phi_TY", c=:black)
#     plot!(ϕ->((sin(ϕ)+0.0*sin(2ϕ))^2)/π, -π, π, lab="M=0 is set to zero", leg=:bottom)
#
#     dϕ = ϕ_data(bin)
#     plot!( dϕ[:,1], dϕ[:,2] ./ sum( dϕ[:,2]) ./ (2π/length(dϕ[:,2])), lab = "M=2 as in the paper")
# end
# savefig(joinpath("plots","dimas_toys_etapi_phi.pdf"))
