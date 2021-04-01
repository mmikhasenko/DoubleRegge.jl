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
using Plots
import Plots.PlotMeasures.mm
theme(:wong2; size=(500,350), bottom_margin=5mm)
# 
using DoubleRegge

#                            _|            
#    _|_|_|    _|_|      _|_|_|    _|_|    
#  _|        _|    _|  _|    _|  _|_|_|_|  
#  _|        _|    _|  _|    _|  _|        
#    _|_|_|    _|_|      _|_|_|    _|_|_|  

# settings_file = joinpath("data", "exp_pro","fit-results_bottom-Po_Np=3.toml")
settings_file = joinpath("data", "exp_pro", "a2Po-f2Po-a2f2-f2f2_opposite-sign", "fit-results_a2Po-f2Po-a2f2-f2f2_opposite-sign_Np=4_alpha=0.8.toml")
! isfile(settings_file) && error("no file")
# 
parsed = TOML.parsefile(settings_file)
@unpack settings, fit_results = parsed
# 
setsystem!(Symbol(settings["system"]))

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

constrained_pw_projection_fixed_model(m, init_pars) = (@show m;
    constrained_pw_projection((cosθ,ϕ)->intensity(m,cosθ,ϕ), init_pars, compass_ηπ_LMs))
#
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
plotmap = map(x->inlims(x.x, (2.4,3.0)), data)
plotdata = data[plotmap]
# 
#
pw_projections = pw_project_fixed_model.(data.x[plotmap])
writedlm(joinpath("data", "exp_pro", "pws_$(settings["tag"]).txt"),
    hcat(unfold.(pw_projections)...))

# @time cPWs_starting_from_compass_pw =
#     [constrained_pw_projection_fixed_model(x,a) for (x,a) in zip(plotdata.x, plotdata.amps)]
# writedlm(joinpath("data", "exp_pro", "constrained_pws_$(settings["tag"])_starting_from_compass_pw.txt"),
#     hcat(unfold.(cPWs_starting_from_compass_pw)...))
# # 
@time cPWs_starting_from_pw = 
    [constrained_pw_projection_fixed_model(x,a) for (x,a) in zip(plotdata.x, pw_projections)]
writedlm(joinpath("data", "exp_pro", "constrained_pws_$(settings["tag"])_starting_from_pw.txt"),
    hcat(unfold.(getproperty.(cPWs_starting_from_pw, :pars))...))
#
# cPWs_starting_from_pw
v = readdlm(joinpath("data", "exp_pro", "constrained_pws_$(settings["tag"])_starting_from_pw.txt"))
#     hcat(unfold.(getproperty.(cPWs_starting_from_pw, :pars))...))
cPWs_starting_from_pw = [fold(v[i,:]) for i in 1:size(v,1)]

# 
pw_intensities = map(x->abs2.(x), pw_projections)
# cpw_intensities = map(x->abs2.(x.pars), cPWs_starting_from_pw)
cpw_intensities = map(x->abs2.(x), [fold(v[i,:]) for i in 1:size(v,1)])


let
    plot(layout=grid(3,3), size=(900,900))
    for (i,(L,M)) in enumerate(LMs)
        scatter!(sp=i, plotdata.x, getindex.(plotdata.I, i),
            yerr=getindex.(plotdata.δI, i), xerr=(plotdata.x[2]-plotdata.x[1])/2,
            c=:black, title="LM=$L$M", ms=3,
            lab=i!=1 ? "" : "data",)
        #
        plot!(sp=i, plotdata.x, getindex.(pw_intensities,i), lab=i!=1 ? "" : "PW projection", l=(2))
        plot!(sp=i, plotdata.x, getindex.(cpw_intensities,i), lab=i!=1 ? "" : "cPW projection", l=(2))
        vspan!(sp=i, fitdata.x[[1,end]], lab="", α=0.1, seriescolor=7)
    end
    plot!(xlab="m(ηπ) (GeV)")
end

phase_with_resp2(x) = arg.(x .* conj(x[2]))
pw_phases = map(phase_with_resp2, pw_projections)
cpw_phases = map(phase_with_resp2, getproperty.(cPWs_starting_from_pw, :pars))
# cpw_phases′ = map(phase_with_resp2, cPWs_starting_from_compass_pw)

function shift(L,M) 
    (L,M) == (1,1) && return 2π
    (L,M) == (3,1) && return 2π
    return 0.0
end

let
    plot(layout=grid(3,3), size=(900,900))
    for (i,(L,M)) in enumerate(LMs)
        scatter!(sp=i, plotdata.x, getindex.(plotdata.ϕ, i),
            yerr=getindex.(plotdata.δϕ, i), xerr=(plotdata.x[2]-plotdata.x[1])/2,
            c=:black, title="LM=$L$M", ms=3,
            lab=i!=1 ? "" : "data",)
        #
        plot!(sp=i, plotdata.x, getindex.(pw_phases,i) .+ shift(L,M), lab=i!=1 ? "" : "PW projection", l=(2))
        plot!(sp=i, plotdata.x, getindex.(cpw_phases,i), lab=i!=1 ? "" : "cPW projection", l=(2))
        # vspan!(sp=i, fitdata.x[[1,end]], lab="", α=0.1, seriescolor=7)
    end
    plot!(xlab="m(ηπ) (GeV)")
end
