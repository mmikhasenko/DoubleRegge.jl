using DrWatson
@quickactivate "DoubleRegge"
using TOML
using TypedTables
using QuadGK
using Cuba
using ForwardDiff
using Optim
using DelimitedFiles
using BenchmarkTools
using Cuba

# 
using LaTeXStrings
using Plots
import Plots.PlotMeasures.mm
theme(:wong2; size=(500,350), bottom_margin=5mm)
# 

using DoubleRegge
setsystem!(:compass_ηπ)

# settings_file = joinpath("data", "exp_pro","fit-results_bottom-Po_Np=3.toml")
settings_file = joinpath("data", "exp_pro", "a2Po-f2f2-PoPo_opposite-sign", "fit-results_a2Po-f2f2-PoPo_opposite-sign_Np=3_alpha=0.8.toml")
# settings_file = joinpath("data", "exp_pro", "a2Po-f2Po-PoPo_opposite-sign", "fit-results_a2Po-f2Po-PoPo_opposite-sign_Np=3_alpha=0.8.toml")
! isfile(settings_file) && error("no file")
# 
parsed = TOML.parsefile(settings_file)
settings = parsed["settings"]
fit_results = parsed["fit_results"]
# 
constrained_pw_projection_fixed_model(m, init_pars) =
    constrained_pw_projection((cosθ,ϕ)->intensity(m,cosθ,ϕ), init_pars, compass_ηπ_LMs)

# fit
const exchanges = sixexchages[settings["exchanges"]]
const model = build_model(exchanges, settings["t2"], settings["scale_α"])
const fixed_pars = fit_results["fit_minimizer"]
fixed_model(m,cosθ,ϕ) = model(m,cosθ,ϕ; pars=fixed_pars)
intensity(m, cosθ, ϕ) = abs2(fixed_model(m, cosθ, ϕ))*q(m)

function pw_project_fixed_model(m)
    amplitude(cosθ,ϕ) = fixed_model(m,cosθ,ϕ)*sqrt(q(m))
    pws = [pw_project(amplitude,L,M) for (L,M) in LMs]
    return pws
end
#
# forward
function pull_through_forward(mpoints)
    x1 = mpoints[1]
    first_init_pars = pw_project_fixed_model(x1)
    cpw = [constrained_pw_projection_fixed_model(x1, first_init_pars)]
    for x in mpoints[2:end]
        cpwi = constrained_pw_projection_fixed_model(x, cpw[end].pars)
        push!(cpw, cpwi)
    end
    return cpw
end
# backward
function pull_through_backward(mpoints)
    x1 = mpoints[end]
    first_init_pars = pw_project_fixed_model(x1)
    cpw = [constrained_pw_projection_fixed_model(x1, first_init_pars)]
    for x in reverse(mpoints[1:end-1])
        cpwi = constrained_pw_projection_fixed_model(x, cpw[end].pars)
        push!(cpw, cpwi)
    end
    return reverse(cpw)
end

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
# 

function morefrequentrange(values, factor)
    len = factor*(length(values)-1)+1
    return range(values[1], values[end], length=len)
end
pull_mpoints = morefrequentrange(plotdata.x, 2)


pw_projections = pw_project_fixed_model.(data.x[plotmap])


# calculations
@time cpw_forward  = pull_through_forward( pull_mpoints)
@time cpw_backward = pull_through_backward(pull_mpoints)


let
    v1 = getproperty.(cpw_forward, :min)
    v2 = getproperty.(cpw_backward, :min)
    plot(pull_mpoints, [v1 v2] .- (v1+v2)/2, seriescolor=[3 4], lw=2, lab=["forwand" "backward"])
    plot!(ylab=L"\Delta \mathcal{L}_\mathrm{ext}", xlab=L"m_{\eta\pi}\,\,(\mathrm{GeV})", title="resolved Barlet branches")
end
savefig(
    joinpath("data", "exp_pro", settings["tag"],
        "cPW_ellh_$(settings["tag"])_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf"))

# 
pw_intensities = map(x->abs2.(x), pw_projections)
cpw_intensities_f = map(x->abs2.(x), getproperty.(cpw_forward, :pars))
cpw_intensities_b = map(x->abs2.(x), getproperty.(cpw_backward, :pars))
let
    plot(layout=grid(3,3), size=(900,900))
    for (i,(L,M)) in enumerate(LMs)
        scatter!(sp=i, plotdata.x, getindex.(plotdata.I, i),
            yerr=getindex.(plotdata.δI, i), xerr=(plotdata.x[2]-plotdata.x[1])/2,
            c=:black, title="LM=$L$M", ms=3,
            lab=i!=1 ? "" : "data",)
        #
        # plot!(sp=i, plotdata.x, getindex.(pw_intensities,i), lab=i!=1 ? "" : "PW projection", l=(2))
        # plot!(sp=i, plotdata.x, getindex.(cpw_intensities, i), lab=i!=1 ? "" : "cPW projection", l=(2))
        plot!(sp=i, pull_mpoints, getindex.(cpw_intensities_f, i), lab=i!=1 ? "" : "cPW projection", l=(2), c=3)
        plot!(sp=i, pull_mpoints, getindex.(cpw_intensities_b, i), lab=i!=1 ? "" : "cPW projection", l=(2), c=4)
        vspan!(sp=i, fitdata.x[[1,end]], lab="", α=0.1, seriescolor=7)
    end
    plot!(xlab=L"m_{\eta\pi}\,\,(\mathrm{GeV})")
end
savefig(
    joinpath("data", "exp_pro", settings["tag"],
        "cPW_invensities_$(settings["tag"])_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf"))

phase_with_resp2(x; δ=zeros(Float64, length(x))) = arg.(x .* conj(x[2]) .* cis.(-δ)) .+ δ
function phase_with_resp2_with_common_shift(complexarray)
    δ = 0.6
    for _ in 1:2
        naivearg = map(x->phase_with_resp2(x; δ=δ), complexarray)
        δ = mean(naivearg)
    end
    return map(x->phase_with_resp2(x; δ=δ), complexarray)
end

pw_phases = phase_with_resp2_with_common_shift(pw_projections)
cpw_phases_f = phase_with_resp2_with_common_shift(getproperty.(cpw_forward, :pars))
cpw_phases_b = phase_with_resp2_with_common_shift(getproperty.(cpw_backward, :pars))
let
    plot(layout=grid(3,3), size=(900,900))
    for (i,(L,M)) in enumerate(LMs)
        i==2 && continue
        scatter!(sp=i, plotdata.x, getindex.(plotdata.ϕ, i),
            yerr=getindex.(plotdata.δϕ, i), xerr=(plotdata.x[2]-plotdata.x[1])/2,
            c=:black, title="LM=$L$M", ms=3,
            lab=i!=2 ? "" : "data",)
        #
        # plot!(sp=i, plotdata.x, getindex.(pw_phases,i), lab=i!=2 ? "" : "PW projection", l=(2))
        plot!(sp=i, pull_mpoints, getindex.(cpw_phases_f,i), lab=i!=2 ? "" : "cPW projection", l=(2), c=3)
        plot!(sp=i, pull_mpoints, getindex.(cpw_phases_b,i), lab=i!=2 ? "" : "cPW projection", l=(2), c=4)
        vspan!(sp=i, fitdata.x[[1,end]], lab="", α=0.1, seriescolor=7)
    end
    plot!(xlab=L"m_{\eta\pi}\,\,(\mathrm{GeV})")
end
savefig(
    joinpath("data", "exp_pro", settings["tag"],
        "cPW_phases_$(settings["tag"])_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf"))

summaryline(t::NamedTuple) = [t.min, real.(t.pars)..., imag.(t.pars)...]

packed_fitresutls_f = hcat(pull_mpoints, hcat(summaryline.(cpw_forward)...)')
packed_fitresutls_b = hcat(pull_mpoints, hcat(summaryline.(cpw_backward)...)')

writedlm(joinpath("data", "exp_pro", "pull_pws_$(settings["tag"]).txt"), [packed_fitresutls_f; packed_fitresutls_b])
