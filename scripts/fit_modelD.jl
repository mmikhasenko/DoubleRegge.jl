using DrWatson
@quickactivate "DoubleRegge"

using DoubleRegge
using QuadGK
using Plots
using Cuba # HCubature
using Optim
theme(:wong; size=(500,350))
#
#
#
# data
const data = [(inten = get_intesity(L,M),
         phase = get_phase(L,M)) for (L,M) in LMs];
const xdata = data[1].inten[:,1];
const Nbins = length(xdata)
#
const amplitudes = constructamps(
    slice(data, :inten, 1:Nbins),
    slice(data, :phase, 1:Nbins)
);

# error band
const lower_bands = readdlm(datadir("exp_pro","bootstrap_lower.txt"));
const upper_bands = readdlm(datadir("exp_pro","bootstrap_upper.txt"));

# range
fit_range = (2.5, 3.0)
filt_fr = fit_range[1] .< xdata .< fit_range[2]
const amplitudes_fr = map(as->as[filt_fr], amplitudes)
const xdata_fr = xdata[filt_fr];
const lower_bands_fr = lower_bands[filt_fr,:];
const upper_bands_fr = upper_bands[filt_fr,:];

# fit
function model(mηπ,cosθ,ϕ; pars)
    vars = (s = DoubleRegge.s0, s1 = mηπ^2,
        cosθ = cosθ, ϕ = ϕ, t2 = -0.2)
    return pars[1]*modelDR(α_a2, α_ℙ, vars; η_forward=true)+
           -pars[2]*modelDR(α_f2, α_ℙ, vars; η_forward=false) +
           -pars[3]*modelDR(α_ℙ, α_ℙ, vars; η_forward=false)
end
#
# model(2.3, cos(π/3),π/4; pars=[1.0, 1.0])
#
function χ2(cosθ,ϕ,c)
    # @show cosθ,ϕ
    Id = abs2.(recamp(cosθ,ϕ, amplitudes_fr))
    Am = model.(xdata_fr,cosθ,ϕ; pars=c)
    Im = abs2.(Am)
    return sum(abs2, Id - Im)
end
#
@time χ2(0.3,0.2,[1,1.1,1e-3]) # test
#
integrate_dcosθdϕ(g) = cuhre((x,f)->f[1]=g(x),2,1).integral[1]
# integrate_dcosθdϕ(g) = quadgk(x->g([x,1/2+1/8]),0,1)[1]
# integrate_dcosθdϕ(g) =  hcubature(g, [0.0, 0.0], [1.0, 1.0])[1]
χ2(pars) = integrate_dcosθdϕ(x->χ2(2x[1]-1, π*(2x[2]-1),pars))*(4π)
#
@time χ2([1,1.1,1e-3]) # test
# @time ForwardDiff.gradient(χ2, [1,1.1,1e-3])
#
@time ft = Optim.minimizer(Optim.optimize(χ2, [1.1,1.1,1e-3], BFGS(),
               Optim.Options(show_trace = true))) # ; autodiff = :forwarddiff

# plotting
let pars = ft# [1,1.2,1e-3] # ft,
    cosθv = range(-1,1, length=100)
    length(cosθv) != size(lower_bands_fr,2) && error("mismatch")
    ps = [let
        calv = dNdcosθ.(cosθv; amps=slice(amplitudes_fr, bin))
        l = lower_bands_fr[bin,:]; u = upper_bands_fr[bin,:]
        # calv = log10.(calv); l = log10.(l); u = log10.(u);
        plot(cosθv, calv,
            ribbon=(calv-l,u-calv), lab=(bin==1 ? "data" : ""),
            lw=2, title="$(round(xdata_fr[bin], digits=2))")
        #
        calv = [quadgk(ϕ->abs2(model(xdata_fr[bin], cosθ, ϕ; pars=pars)), -π, π)[1] for cosθ in cosθv]
        # calv = log10.(calv)
        plot!(cosθv, calv, lab=(bin==1 ? "model" : ""), lc=:black, lw=2)
    end for bin in 1:length(xdata_fr)]
    plot(ps..., size=(1500,1000))
end
# savefig(plotsdir("P+a2.pdf"))
# savefig(plotsdir("f2+a2.pdf"))
savefig(plotsdir("f2+a2+P.pdf"))
# savefig(plotsdir("f2+a2.log10.pdf"))

let bin = 1
    cosθv = range(-1,1,length=50)
    ϕv = range(-π,π,length=50)
    calv_d = [abs2.(recamp(cosθ,ϕ,slice(amplitudes_fr, bin))) for (ϕ, cosθ) in Iterators.product(ϕv, cosθv)]
    calv_m = [abs2.(model(xdata_fr[bin],cosθ,ϕ; pars=ft)) for (ϕ, cosθ) in Iterators.product(ϕv, cosθv)]
    plot(heatmap(cosθv, ϕv, calv_d),
         heatmap(cosθv, ϕv, calv_m),
         heatmap(cosθv, ϕv, (calv_d-calv_m).^2),
         size=(1400,400), layout=grid(1,3), title=["data" "model" "diff^2"])
end
