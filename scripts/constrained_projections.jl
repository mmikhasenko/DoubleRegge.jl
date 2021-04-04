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

# # # # # # # # # # # # # # # # # # # # 
# 
tag = "etapi_a2Po-f2Po-a2f2-f2f2_opposite-sign"
# 
# # # # # # # # # # # # # # # # # # # # 

settings_file = fitsfolder(tag, "fit-results.toml")
! isfile(settings_file) && error("no file")
parsed = TOML.parsefile(settings_file)
@unpack settings, fit_results = parsed

# 
setsystem!(Symbol(settings["system"]))
description = settings["system"] == "compass_ηπ" ? description_ηπ :
                (settings["system"] == "compass_η′π" ? description_η′π :
                    error("unknown system $(settings["system"])"))


# build model
const model = build_model(
    sixexchages[settings["exchanges"]],
    settings["t2"],
    settings["scale_α"])
const fixed_pars = fit_results["fit_minimizer"]
fixed_model(m,cosθ,ϕ) = model(m,cosθ,ϕ; pars=fixed_pars)
intensity(m, cosθ, ϕ) = abs2(fixed_model(m, cosθ, ϕ))*q(m)

# get data
data = read_data(settings["pathtodata"], description)
# fit range
fitdata = filter(data) do x
    inlims(x.x, settings["fitrange"])
end
# plot 
plotdata = filter(data) do x
    inlims(x.x, (2.4,3.0))
end

# 
const used_LMs = Array(data[1].Iϕ.LMs)
function pw_project_fixed_model(m)
    amplitude(cosθ,ϕ) = fixed_model(m,cosθ,ϕ)*sqrt(q(m))
    pws = [pw_project(amplitude,L,M) for (L,M) in used_LMs]
    return pws
end
constrained_pw_projection_fixed_model(m, init_pars) = (@show m;
    constrained_pw_projection((cosθ,ϕ)->intensity(m,cosθ,ϕ), init_pars, used_LMs))

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
pull_through_backward(mpoints) = reverse(pull_through_forward(reverse(mpoints)))
#
function morefrequentrange(values, factor)
    len = factor*(length(values)-1)+1
    return range(values[1], values[end], length=len)
end

#  _|              _|                                    
#      _|_|_|    _|_|_|_|    _|_|    _|_|_|      _|_|_|  
#  _|  _|    _|    _|      _|_|_|_|  _|    _|  _|_|      
#  _|  _|    _|    _|      _|        _|    _|      _|_|  
#  _|  _|    _|      _|_|    _|_|_|  _|    _|  _|_|_|    

pull_mpoints = morefrequentrange(plotdata.x, 5)
# 
pw_projections = pw_project_fixed_model.(pull_mpoints)
writedlm(fitsfolder(tag,"PWs.txt"), hcat(unfold.(pw_projections)...))

@time cPWs_starting_from_pw = 
    [constrained_pw_projection_fixed_model(x,a) for (x,a) in zip(pull_mpoints, pw_projections)]
writedlm(fitsfolder(tag,"cPWs.txt"),
    hcat(unfold.(getproperty.(cPWs_starting_from_pw, :pars))...))
#

# ambibuities
@time cpw_forward  = pull_through_forward( pull_mpoints)
writedlm(fitsfolder(tag,"cPWs_forward.txt"),
    hcat(unfold.(getproperty.(cpw_forward, :pars))...))
# 
@time cpw_backward = pull_through_backward(pull_mpoints)
writedlm(fitsfolder(tag,"cPWs_backward.txt"),
    hcat(unfold.(getproperty.(cpw_backward, :pars))...))

#            _|              _|      _|      _|                      
#  _|_|_|    _|    _|_|    _|_|_|_|_|_|_|_|      _|_|_|      _|_|_|  
#  _|    _|  _|  _|    _|    _|      _|      _|  _|    _|  _|    _|  
#  _|    _|  _|  _|    _|    _|      _|      _|  _|    _|  _|    _|  
#  _|_|_|    _|    _|_|        _|_|    _|_|  _|  _|    _|    _|_|_|  
#  _|                                                            _|  
#  _|                                                        _|_|    

function read_PW_matrices(filename, LMs)
    matrix = readdlm(filename)
    # 
    amplitudevectors = [fold(matrix[:,i]) for i in 1:size(matrix,2)]
    return TwoBodyPartialWaves.(Ref(LMs), amplitudevectors)
end

PWs = changerepresentation.(read_PW_matrices(fitsfolder(tag,"PWs.txt"), used_LMs); iref=2)
cPWs = changerepresentation.(read_PW_matrices(fitsfolder(tag,"cPWs.txt"), used_LMs); iref=2)
cPWs_f = changerepresentation.(read_PW_matrices(fitsfolder(tag,"cPWs_forward.txt"), used_LMs); iref=2)
cPWs_b = changerepresentation.(read_PW_matrices(fitsfolder(tag,"cPWs_backward.txt"), used_LMs); iref=2)
# 
data_intensities = [(((plotdata.Iϕ)..:PWs)..i)..:I for i in 1:length(used_LMs)]
pw_intensities = [((PWs..:PWs)..i)..:I for i in 1:length(PWs[1].LMs)]
cpw_intensities = [((cPWs..:PWs)..i)..:I for i in 1:length(cPWs[1].LMs)]
cpw_f_intensities = [((cPWs_f..:PWs)..i)..:I for i in 1:length(cPWs[1].LMs)]
cpw_b_intensities = [((cPWs_b..:PWs)..i)..:I for i in 1:length(cPWs[1].LMs)]

let
    N = 3
    M = div(length(used_LMs)-1,3)+1
    plot(layout=grid(M,N), size=(300*N,300*M))
    # 
    for (i,(L,M)) in enumerate(used_LMs)
        scatter!(sp=i, plotdata.x, data_intensities[i],
            xerr=(plotdata.x[2]-plotdata.x[1])/2,
            c=:black, title="LM=$L$M", ms=3,
            lab=i!=1 ? "" : "data",)
        #
        plot!(sp=i, pull_mpoints, pw_intensities[i], lab=i!=1 ? "" : "PW projection", l=(2))
        plot!(sp=i, pull_mpoints, cpw_intensities[i], lab=i!=1 ? "" : "cPW projection", l=(2))
        # 
        plot!(sp=i, pull_mpoints, cpw_f_intensities[i], lab="", l=(1,:gray))
        plot!(sp=i, pull_mpoints, cpw_b_intensities[i], lab="", l=(1,:gray))
        # vspan!(sp=i, fitdata.x[[1,end]], lab="", α=0.1, seriescolor=7)
    end
    plot!(xlab="m(ηπ) (GeV)")
end


data_phases = [[p.Iϕ.PWs[i].ϕ for p in plotdata] for i in 1:length(used_LMs)]
data_phases_adj = meanshiftbyperiod.(data_phases, 0)
# 
pw_phases = alignperiodicsequence.([((PWs..:PWs)..i)..:ϕ for i in 1:length(used_LMs)])
pw_phases_adj = meanshiftbyperiod.(pw_phases, mean.(data_phases_adj)..:val)
# 
cpw_phases = alignperiodicsequence.([((cPWs..:PWs)..i)..:ϕ for i in 1:length(used_LMs)])
cpw_phases_adj = meanshiftbyperiod.(cpw_phases, mean.(data_phases_adj)..:val)
# 
cpw_f_phases = alignperiodicsequence.([((cPWs_f..:PWs)..i)..:ϕ for i in 1:length(used_LMs)])
cpw_f_phases_adj = meanshiftbyperiod.(cpw_f_phases, mean.(data_phases_adj)..:val)
# 
cpw_b_phases = alignperiodicsequence.([((cPWs..:PWs)..i)..:ϕ for i in 1:length(used_LMs)])
cpw_b_phases_adj = meanshiftbyperiod.(cpw_b_phases, mean.(data_phases_adj)..:val)

let
    N = 3
    M = div(length(used_LMs)-1,3)+1
    plot(layout=grid(M,N), size=(300*N,300*M))
    # 
    for (i,(L,M)) in enumerate(used_LMs)
        scatter!(sp=i, plotdata.x, data_phases_adj[i],
            xerr=(plotdata.x[2]-plotdata.x[1])/2,
            c=:black, title="LM=$L$M", ms=3,
            lab=i!=1 ? "" : "data",)
        #
        plot!(sp=i, pull_mpoints, pw_phases_adj[i], lab=i!=1 ? "" : "PW projection", l=(2))
        plot!(sp=i, pull_mpoints, cpw_phases_adj[i], lab=i!=1 ? "" : "cPW projection", l=(2))
        #
        plot!(sp=i, pull_mpoints, cpw_f_phases_adj[i], lab="", l=(1,:gray))
        plot!(sp=i, pull_mpoints, cpw_b_phases_adj[i], lab="", l=(1,:gray))
        # vspan!(sp=i, fitdata.x[[1,end]], lab="", α=0.1, seriescolor=7)
    end
    plot!(xlab="m(ηπ) (GeV)")
end
# 
writedlm(fitsfolder(tag,"data_ajusted.txt"),
    hcat(plotdata.x,
        hcat(data_intensities...)..:val,
        hcat(data_intensities...)..:err,
        hcat(data_phases_adj...)..:val,
        hcat(data_phases_adj...)..:err))
# 
writedlm(fitsfolder(tag,"pw_ajusted.txt"),
    hcat(plotdata.x,
        hcat(pw_intensities...),
        hcat(cpw_intensities...),
        hcat(pw_phases_adj...),
        hcat(cpw_phases_adj...)))
