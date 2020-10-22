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
# 
using Statistics
using LinearAlgebra

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
const pfr = fit_results["fit_minimizer"]

# 
setsystem!(Symbol(settings["system"]))

# fit
const exchanges = sixexchages[settings["exchanges"]]
const model = build_model(exchanges, settings["t2"], settings["scale_α"])
fixed_model(m,cosθ,ϕ; pars=pfr) = model(m,cosθ,ϕ; pars=pars)
fixed_model_sqrtq(m,cosθ,ϕ; pars=pfr) = fixed_model(m,cosθ,ϕ; pars=pars)*sqrt(q(m))
intensity(m, cosθ, ϕ; pars=pfr) = abs2(fixed_model_sqrtq(m,cosθ,ϕ; pars=pars))
# 
function pw_project_fixed(m::Float64,L,M)
    amplitude(cosθ,ϕ) = fixed_model_sqrtq(m,cosθ,ϕ)
    return pw_project(amplitude,L,M)
end
# 
model_integral(m)          = integrate_dcosθdϕ((cosθ,ϕ)->abs2(fixed_model(m,cosθ,ϕ)))[1]*q(m)
model_integral_forward(m)  = integrate_dcosθdϕ((cosθ,ϕ)->abs2(fixed_model(m,cosθ,ϕ)),(0,1))[1]*q(m)
model_integral_backward(m) = integrate_dcosθdϕ((cosθ,ϕ)->abs2(fixed_model(m,cosθ,ϕ)),(-1,0))[1]*q(m)


# data
const LMs = compass_ηπ_LMs
data = Table(x_IδI_ϕδϕ_compass_ηπ(settings["pathtodata"]))
amplitudes = [sqrt.(is) .* cis.(ϕs) for (is,ϕs) in zip(data.I, data.ϕ)]
data = Table(data, amps=amplitudes)
# fit range
fitrangemap = map(x->inlims(x.x, settings["fitrange"]), data)
fitdata = data[fitrangemap]
# plot 
plotmap = map(x->inlims(x.x, (2.4,3.0)), data)
plotdata = data[plotmap]


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
    plot(size=(900,350), layout=grid(1,2),
        xlab="m(ηπ) (GeV)",
        ylab=["intensity" "(F-B) / (F+B)" "intensity" "intensity"],
        title=["number of events" "asymmetry"])
    #
    common_options = (xerr=(plotdata.x[2]-plotdata.x[1])/2, m=(3,))
    # 
    scatter!(sp=1, plotdata.x, data_forward_in_bins; lab="forward",
        common_options..., yerr = sqrt.(sum.(x->x^2, plotdata.δI)/2), seriescolor=3)
    scatter!(sp=1, plotdata.x, data_backward_in_bins; lab="backward",
        common_options..., yerr = sqrt.(sum.(x->x^2, plotdata.δI)/2), seriescolor=4)
    plot!(sp=1, plotdata.x, intensity_forward_in_bins, lab="", lw=2, seriescolor=3)
    plot!(sp=1, plotdata.x, intensity_backward_in_bins, lab="", lw=2, seriescolor=4)
    # 
    scatter!(sp=2, plotdata.x, asymm_data; lab="", common_options..., yerr=δasymm_data)
    plot!(sp=2, plotdata.x, asymm_model, lab="", ylim=(-1,1), lw=2)
    #
    scatter!(sp=1, plotdata.x, sum.(plotdata.I);
        common_options..., lab="intensity",
        yerr = sqrt.(sum.(x->x^2, plotdata.δI)), seriescolor=2)
    plot!(sp=1, plotdata.x, intensity_in_bins, lw=2, lab="", seriescolor=2)
    # 
    vspan!(sp=1, fitdata.x[[1,end]], lab="", α=0.1, seriescolor=7)
    vspan!(sp=2, fitdata.x[[1,end]], lab="", α=0.1, seriescolor=7)
end
savefig(
    joinpath("data", "exp_pro",
        "intensity-assymetry_$(settings["tag"])_Np=$(length(settings["exchanges"])).pdf"))

# cosθ distributions
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
savefig(
    joinpath("data", "exp_pro",
        "cos-distributions_$(settings["tag"])_Np=$(length(settings["exchanges"])).pdf"))


# # projections
pw_projections = [map(LM->pw_project_fixed(m,LM...), LMs) for m in data.x[plotmap]]
pw_intensities = map(x->abs2.(x), pw_projections)
let
    plot(layout=grid(3,3), size=(900,900))
    for (i,(L,M)) in enumerate(LMs)
        scatter!(sp=i, plotdata.x, getindex.(plotdata.I, i),
            yerr=getindex.(plotdata.δI, i), xerr=(plotdata.x[2]-plotdata.x[1])/2,
            c=:black, title="LM=$L$M", ms=3,
            lab=i!=1 ? "" : "data",)
        #
        plot!(sp=i, plotdata.x, getindex.(pw_intensities,i), lab=i!=1 ? "" : "PW projection", l=(2,))
        vspan!(sp=i, fitdata.x[[1,end]], lab="", α=0.1, seriescolor=7)
    end
    plot!(xlab="m(ηπ) (GeV)")
end
savefig(
    joinpath("data", "exp_pro",
        "pw-projections_$(settings["tag"])_Np=$(length(settings["exchanges"])).pdf"))
#
# Odd and Even waves
fHeigher = (intensity_in_bins .- sum.(pw_intensities)) ./ intensity_in_bins
filtodd = isodd.(getproperty.(LMs,:L))
fOdd = map(x->sum(x .* filtodd), pw_intensities) ./ intensity_in_bins
fEven = map(x->sum(x .* iszero.(filtodd)), pw_intensities) ./ intensity_in_bins
let
    plot(ylab="fraction", xlab="m(ηπ) (GeV)", size=(500,350), title=settings["tag"])
    plot!(plotdata.x, fOdd, lab="Odd waves L ≤ 6")
    plot!(plotdata.x, fEven, lab="Even waves L ≤ 6")
    plot!(plotdata.x, fHeigher, lab="Higher waves L > 6")
    plot!(ylims=(0,1))
end
savefig(
    joinpath("data", "exp_pro",
        "odd-and-even_$(settings["tag"])_Np=$(length(settings["exchanges"])).pdf"))
