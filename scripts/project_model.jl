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
using LaTeXStrings
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

# settings_file = joinpath("data", "exp_pro", "a2Po-f2Po-PoPo_same-sign", "fit-results_a2Po-f2Po-PoPo_same-sign_Np=3_alpha=0.8.toml")
# settings_file = joinpath("data", "exp_pro", "a2Po-f2Po-PoPo_opposite-sign", "fit-results_a2Po-f2Po-PoPo_opposite-sign_Np=3_alpha=0.8.toml")
settings_file = joinpath("data", "exp_pro", "a2Po-f2f2-PoPo_opposite-sign", "fit-results_a2Po-f2f2-PoPo_opposite-sign_Np=3_alpha=0.8.toml")
! isfile(settings_file) && error("no file")
# 
parsed = TOML.parsefile(settings_file)
settings = parsed["settings"]
fit_results = parsed["fit_results"]
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
model_integral(m; pars=pfr)          = integrate_dcosθdϕ((cosθ,ϕ)->abs2(fixed_model_sqrtq(m,cosθ,ϕ; pars=pars)))[1]
model_integral_forward(m; pars=pfr)  = integrate_dcosθdϕ((cosθ,ϕ)->abs2(fixed_model_sqrtq(m,cosθ,ϕ; pars=pars)),(0,1))[1]
model_integral_backward(m; pars=pfr) = integrate_dcosθdϕ((cosθ,ϕ)->abs2(fixed_model_sqrtq(m,cosθ,ϕ; pars=pars)),(-1,0))[1]

δ(i; n=3) = (1:n .== i)
integral_interf(m,i,j) = model_integral(m; pars=pfr.*δ(i)+pfr.*δ(j)) - 
    model_integral(m; pars=pfr.*δ(i)) -
    model_integral(m; pars=pfr.*δ(j));
# 
contribution_matrix(m) = [(i>j ? 0 : integral_interf(m,i,j)/(i==j ? 2 : 1)) for i in 1:3, j in 1:3]
let
    m = sum(contribution_matrix.(fitdata.x[1:2]))
    heatmap([mi==0.0 ? NaN : mi for mi in m],
        xticks = (1:3, getindex.(exchanges,4)),
        yticks = (1:3, getindex.(exchanges,4)), colorbar=false)
    # 
    mn = m ./ sum(m)
    for i in Iterators.CartesianIndices(m)
        (i[1]>i[2]) && continue
        annotate!([(i[2], i[1],
            text("$(Int(round(m[i], digits=0))) / $(Int(round(100*mn[i], digits=0)))%", 10, :red))])
    end
    plot!(size=(400,350), title="contributions of different diagrams")
    ellh = fit_results["fit_minimum"]
    s = prod("$l: $v,\n" for (l,v) in zip(getindex.(exchanges,4), round.(pfr, digits=2))) 
    s *= "extLLH: "*string(round(ellh, digits=1))
    annotate!([(0.6,3,text(s, 10, :left))])
end
savefig(
    joinpath("data", "exp_pro", settings["tag"],
        "contributions_$(settings["tag"])_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf"))
# 

# data
const LMs = compass_ηπ_LMs
data = Table(x_IδI_ϕδϕ_compass_ηπ(settings["pathtodata"]))
amplitudes = [sqrt.(is) .* cis.(ϕs) for (is,ϕs) in zip(data.I, data.ϕ)]
data = Table(data, amps=amplitudes)
# fit range
fitrangemap = map(x->inlims(x.x, settings["fitrange"]), data)
fitdata = data[fitrangemap]
# plot 
plotmap = map(x->inlims(x.x, (2.2,3.0)), data)
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
    joinpath("data", "exp_pro", settings["tag"],
        "intensity-assymetry_$(settings["tag"])_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf"))
#
# let bin = 1
#     mηπ = fitdata.x[bin]
#     cosθv = range(-1,1, length=100)
#     ϕv = range(-π, π, length=95)
#     calv = intensity.(mηπ,cosθv',ϕv)
#     heatmap(cosθv, ϕv, calv, colorbar=false,
#         xlab=L"\cos\theta", ylab=L"\phi", title="$(round(mηπ, digits=2)) GeV")
# end

setfirstargument(x1,f) = (x2,x3)->f(x1,x2,x3)
setthirdfourtharguments(x3,x4,f) = (x1,x2)->f(x1,x2,x3,x4)
# 
function phi_moment(f, l; forward=true)
    cosrange = forward ? (0,1) : (-1,0)
    return integrate_dcosθdϕ((cosθ,ϕ)->cos(l*ϕ)*f(cosθ,ϕ), cosrange)[1] /
           integrate_dcosθdϕ(f, cosrange)[1]
end

phi_moment_data(amps; l, forward) = phi_moment.(
    setthirdfourtharguments.(amps, Ref(LMs),(x...)->abs2(recamp(x...))),l;
        forward=forward)
#

function take_stdiv_of_vector(vector_of_vector)
    mat = hcat(vector_of_vector...)
    err = sqrt.(diag(cov(mat; dims=2)))
    return err
end

# phi asymmetry
function phiasymmplot(l)
    plot(xlab=L"m_{\eta\pi}\,(\textrm{GeV})", ylab=(l==1 ? L"<\cos\,\phi>" : L"<\cos\,2\phi>"), size=(500,350))
    plot!(plotdata.x, phi_moment.(setfirstargument.(plotdata.x, intensity),l; forward=true ), lw=2, c=2, lab="forward")
    plot!(plotdata.x, phi_moment.(setfirstargument.(plotdata.x, intensity),l; forward=false), lw=2, c=3, lab="backward")
    # let
    vals = phi_moment_data(plotdata.amps; l=l, forward=true)
    err = take_stdiv_of_vector([phi_moment_data(randA(plotdata); l=l, forward=true) for _ in 1:100])
    scatter!(plotdata.x, vals, yerr=err, c=2, lab="")
    #
    vals = phi_moment_data(plotdata.amps; l=l, forward=false)
    err = take_stdiv_of_vector([phi_moment_data(randA(plotdata); l=l, forward=false) for _ in 1:100])
    scatter!(plotdata.x, vals, yerr=err, c=3, lab="")
end
plot(phiasymmplot(1), phiasymmplot(2), size=(900,350), layout=grid(1,2))
vspan!(sp=1, fitdata.x[[1,end]], lab="", α=0.1, seriescolor=7)
vspan!(sp=2, fitdata.x[[1,end]], lab="", α=0.1, seriescolor=7)
savefig(
    joinpath("data", "exp_pro", settings["tag"],
        "intensity-cosphi_$(settings["tag"])_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf"))

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
        plot!(cosθv, projection.(cosθv), lab="", lw=2, yaxis=nothing)
    end
    ps = make_plot.(1:length(fitdata.x))
    plot(ps..., size=(1100,500), layout=grid(3,5))
end
savefig(
    joinpath("data", "exp_pro", settings["tag"],
        "cos-distributions_$(settings["tag"])_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf"))
#


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
    joinpath("data", "exp_pro", settings["tag"],
        "pw-projections_$(settings["tag"])_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf"))
#
tolab(LM) = "$(LM.L)$(LM.M)"

bar([sum(getindex.(fitdata.I, i)) for (i,LM) in enumerate(LMs)],
    xticks=(1:7, tolab.(LMs)), xlab="LM", ylab="intensity", lab="")

bar([sum(getindex.(pw_intensities[6:end], i)) for (i,LM) in enumerate(LMs)],
    xticks=(1:7, tolab.(LMs)), xlab="LM", ylab="intensity", lab="")
#

# Odd and Even waves
fHeigher = (intensity_in_bins .- sum.(pw_intensities)) ./ intensity_in_bins
filtodd = isodd.(getproperty.(LMs,:L))
fOdd = map(x->sum(x .* filtodd), pw_intensities) ./ intensity_in_bins
fEven = map(x->sum(x .* iszero.(filtodd)), pw_intensities) ./ intensity_in_bins
# 
fOdd_compass = map(x->sum(x .* filtodd), plotdata.I) ./ sum.(plotdata.I)
fEven_compass = map(x->sum(x .* iszero.(filtodd)), plotdata.I) ./ sum.(plotdata.I)

let
    plot(ylab="fraction", xlab="m(ηπ) (GeV)", size=(500,350), title=settings["tag"])
    plot!(plotdata.x, fHeigher, lab="Higher waves L > 6", lw=2)
    plot!(plotdata.x, fEven, lab="Even waves L ≤ 6", lw=2)
    plot!(plotdata.x, fOdd, lab="Odd waves L ≤ 6", lw=2)
    scatter!(plotdata.x, xerr=(plotdata.x[2]-plotdata.x[1])/2, fEven_compass, lab="", mc=2, ms=3)
    scatter!(plotdata.x, xerr=(plotdata.x[2]-plotdata.x[1])/2, fOdd_compass, lab="", mc=3, ms=3)
    plot!(ylims=(0,1), leg=:left)
    vspan!(fitdata.x[[1,end]], lab="", α=0.1, seriescolor=7)
end
savefig(
    joinpath("data", "exp_pro", settings["tag"],
        "odd-and-even_$(settings["tag"])_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf"))
#


let 
    pathtofolder = joinpath("data", "exp_pro", settings["tag"])
    inputfiles = readdir(pathtofolder)
    outputfile = joinpath(pathtofolder, "_combined.pdf")
    inputfiles = filter(f->splitext(f)[2]==".pdf" && f!=outputfile, inputfiles)
end


produced_files = [
    "intensity-assymetry_$(settings["tag"])_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf",
    "cos-distributions_$(settings["tag"])_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf",
    "contributions_$(settings["tag"])_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf",
    "odd-and-even_$(settings["tag"])_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf",
    "intensity-cosphi_$(settings["tag"])_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf",
    "pw-projections_$(settings["tag"])_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf"
    ]

let 
    pathtofolder = joinpath("data", "exp_pro", settings["tag"])
    # inputfiles = readdir(pathtofolder, join=true)
    outputfile = joinpath(pathtofolder, "combined_$(settings["tag"])_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf")
    # inputfiles = filter(f->splitext(f)[2]==".pdf" && f!=outputfile, inputfiles)
    inputfiles = joinpath.(Ref(pathtofolder), produced_files)
    #
    run(`pdftk $inputfiles cat output $outputfile`)
end
