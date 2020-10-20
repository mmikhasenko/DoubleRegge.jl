using DrWatson
@quickactivate "DoubleRegge"
using TOML
using TypedTables
using QuadGK
using Cuba
# 
using Plots
import Plots.PlotMeasures.mm
theme(:wong2; size=(500,350), bottom_margin=5mm)
# 
using DoubleRegge
setsystem!(:compass_ηπ)

#                            _|            
#    _|_|_|    _|_|      _|_|_|    _|_|    
#  _|        _|    _|  _|    _|  _|_|_|_|  
#  _|        _|    _|  _|    _|  _|        
#    _|_|_|    _|_|      _|_|_|    _|_|_|  

# settings_file = joinpath("data", "exp_pro","fit-results_bottom-Po_Np=3.toml")
settings_file = joinpath("data", "exp_pro","fit-results_a2Po-f2f2-PoPo_Np=3.toml")
! isfile(settings_file) && error("no file")
# 
parsed = TOML.parsefile(settings_file)
settings = parsed["settings"]
const fit_results = parsed["fit_results"]

# 
setsystem!(Symbol(settings["system"]))

# fit
const exchanges = sixexchages[settings["exchanges"]]
const model = build_model(exchanges, settings["t2"], settings["scale_α"])
fixed_model(m,cosθ,ϕ; pars=fit_results["fit_minimizer"]) = model(m,cosθ,ϕ; pars=pars)
intensity(m, cosθ, ϕ; pars=fit_results["fit_minimizer"]) = abs2(fixed_model(m, cosθ, ϕ; pars=pars))*q(m)

# data
LMs = compass_ηπ_LMs
data = Table(x_IδI_ϕδϕ_compass_ηπ(settings["pathtodata"]))
amplitudes = [sqrt.(is) .* cis.(ϕs) for (is,ϕs) in zip(data.I, data.ϕ)]
data = Table(data, amps=amplitudes)
# fit range
fitrangemap = map(x->inlims(x.x, settings["fitrange"]), data)
fitdata = data[fitrangemap]
# plot 
plotmap = map(x->inlims(x.x, (2.4,3.0)), data)
plotdata = data[plotmap]


function pw_iϕ(m,L,M)
    amplitude(cosθ,ϕ) = fixed_model(m,cosθ,ϕ)*sqrt(q(m))
    return pw_project(amplitude,L,M)
end
#
model_integral(m) = integrate_dcosθdϕ((cosθ,ϕ)->abs2(fixed_model(m,cosθ,ϕ)))*q(m)
model_integral_forward(m) = integrate_dcosθdϕ((cosθ,ϕ)->abs2(fixed_model(m,cosθ,ϕ)),(0,1))*q(m)
model_integral_backward(m) = integrate_dcosθdϕ((cosθ,ϕ)->abs2(fixed_model(m,cosθ,ϕ)),(-1,0))*q(m)

# ellh
fit_results["fit_minimum"]

# constraint
sum(sum, fitdata.I), sum(model_integral.(fitdata.x))

# intensity
intensity_in_bins = model_integral.(plotdata.x)

# asymmetry
intensity_forward_in_bins = model_integral_forward.(plotdata.x)
intensity_backward_in_bins = model_integral_backward.(plotdata.x)
data_forward_in_bins = [quadgk(cosθ->dNdcosθ(cosθ; amps=a, LMs=LMs),0,1)[1] for a in plotdata.amps]
data_backward_in_bins = [quadgk(cosθ->dNdcosθ(cosθ; amps=a, LMs=LMs),-1,0)[1] for a in plotdata.amps]
# 
asymmetry(f,b) = (f-b)/(f+b)
δasymmetry(f,b,δ) = 2δ*sqrt(f^2+b^2)/(f+b)^2

asymm_data = [asymmetry(f,b) for (f,b) in zip(data_forward_in_bins,data_backward_in_bins)]
δasymm_data = [δasymmetry(f,b,δ) for (f,b,δ) in zip(
    data_forward_in_bins,
    data_backward_in_bins,
    sqrt.(sum.(x->x^2, plotdata.δI)/2))]
# 
asymm_model = [asymmetry(f,b) for (f,b) in zip(intensity_forward_in_bins, intensity_backward_in_bins)]
# 
let
    plot(size=(900,700), layout=grid(2,2),
        xlab="m(ηπ) (GeV)",
        ylab=["intensity" "(F-B) / (F+B)" "intensity" "intensity"],
        title=["total" "asymmetry" "forward" "backward"])
    #
    common_options = (xerr=(plotdata.x[2]-plotdata.x[1])/2, m=(3,:black))
    # 
    scatter!(sp=3, plotdata.x, data_forward_in_bins; lab="data",
        common_options..., yerr = sqrt.(sum.(x->x^2, plotdata.δI)/2))
    scatter!(sp=4, plotdata.x, data_backward_in_bins; lab=""   ,
        common_options..., yerr = sqrt.(sum.(x->x^2, plotdata.δI)/2))
    plot!(sp=3, plotdata.x, intensity_forward_in_bins, lab="model", lw=2)
    plot!(sp=4, plotdata.x, intensity_backward_in_bins, lab="", lw=2)
    # 
    scatter!(sp=2, plotdata.x, asymm_data; lab="", common_options..., yerr=δasymm_data)
    plot!(sp=2, plotdata.x, asymm_model, lab="", ylim=(-1,1), lw=2)
    #
    scatter!(sp=1, plotdata.x, sum.(plotdata.I);
        common_options..., lab="",
        yerr = sqrt.(sum.(x->x^2, plotdata.δI)))
    plot!(sp=1, plotdata.x, intensity_in_bins, lw=2, lab="")
    # 
    vspan!(sp=1, fitdata.x[[1,end]], lab="", α=0.1, seriescolor=7)
    vspan!(sp=2, fitdata.x[[1,end]], lab="", α=0.1, seriescolor=7)
    vspan!(sp=3, fitdata.x[[1,end]], lab="", α=0.1, seriescolor=7)
    vspan!(sp=4, fitdata.x[[1,end]], lab="", α=0.1, seriescolor=7)
end

# cosθ distributions
using Statistics
using LinearAlgebra

let
    cosθv = range(-1,1, length=101)
    function make_plot(bin)
        calv = dNdcosθ.(cosθv; amps=fitdata.amps[bin], LMs=LMs)
        #
        stat = hcat([dNdcosθ.(cosθv; amps=randA(fitdata)[bin], LMs=LMs) for _ in 1:1000]...)
        err = sqrt.(diag(cov(stat; dims=2)))
        # 
        plot(xlab="cosθ", title="$(round(fitdata.x[bin]; digits=2)) GeV")
        plot!(cosθv, calv, ribbon=err, lab="")
        # 
        projection(cosθ) = quadgk(ϕ->intensity(fitdata.x[bin], cosθ, ϕ), -π, π)[1]
        plot!(cosθv, projection.(cosθv), lab="", lw=2)
    end
    ps = make_plot.(1:length(fitdata.x))
    plot(ps..., size=(900,600))
end
# savefig()


# projections
pw_projections = [pw_iϕ.(data.x[plotmap],L,M) for (L,M) in LMs]
let
    plot(layout=grid(3,3), size=(900,900))
    for (i,(m,(L,M))) in enumerate(zip(pw_projections,LMs))
        scatter!(sp=i, plotdata.x, getindex.(plotdata.I, i),
            yerr=getindex.(plotdata.δI, i), xerr=(plotdata.x[2]-plotdata.x[1])/2,
            c=:black, title="LM=$L$M", ms=3,
            lab=i!=1 ? "" : "data",)
        plot!(sp=i, plotdata.x, getindex.(m,1), lab=i!=1 ? "" : "a2Po-f2Po-PoPo", l=(2))
        vspan!(sp=i, fitdata.x[[1,end]], lab="", α=0.1, seriescolor=7)
    end
    plot!(xlab="m(ηπ) (GeV)")
end

savefig(
    joinpath("data", "exp_pro",
        "pw-projections_$(settings["tag"])_Np=$(length(settings["exchanges"])).pdf"))
#

