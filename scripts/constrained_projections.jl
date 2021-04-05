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
using LaTeXStrings
# 
using DoubleRegge


#                                _|                      
#  _|  _|_|    _|_|      _|_|_|      _|_|_|      _|_|    
#  _|_|      _|_|_|_|  _|        _|  _|    _|  _|_|_|_|  
#  _|        _|        _|        _|  _|    _|  _|        
#  _|          _|_|_|    _|_|_|  _|  _|_|_|      _|_|_|  
#                                    _|                  
#                                    _|                  

function intensity(PWs::Vector{TwoBodyPartialWaves{N, NamedTuple{(:I,:ϕ),V}}} where N where V, i::Integer)
    (((PWs)..:PWs)..i)..:I
end
function phase(PWs::Vector{TwoBodyPartialWaves{N, NamedTuple{(:I,:ϕ),V}}} where N where V, i::Integer)
    phases = alignperiodicsequence(((PWs..:PWs)..i)..:ϕ)
    phases_adj = meanshiftbyperiod(phases)
end

@recipe function f(x, PWs::Vector{TwoBodyPartialWaves{N, V}} where N where V <: Number,
    what::Symbol, i::Integer; iref=2)
    PWsIϕ = changerepresentation.(PWs; iref=iref)
    (x,PWsIϕ,what,i)
end

@recipe function f(x, PWs::Vector{TwoBodyPartialWaves{N, NamedTuple{(:I,:ϕ),V}}} where N where V,
    what::Symbol, i::Integer)
    y = []
    label --> ""
    L,M = PWs[1].LMs[i]
    title --> "LM=$L$M"
    if what == :I
        intensities = intensity(PWs, i)
        y = intensities
    end
    if what == :ϕ
        phases_adj = phase(PWs, i)
        if mean(phases_adj) < π/2
            phases_adj .+= 2π
        end
        y = phases_adj .* (180/π)
    end
    (x,y)
end

#                            _|            
#    _|_|_|    _|_|      _|_|_|    _|_|    
#  _|        _|    _|  _|    _|  _|_|_|_|  
#  _|        _|    _|  _|    _|  _|        
#    _|_|_|    _|_|      _|_|_|    _|_|_|  


# # # # # # # # # # # # # # # # # # # # 
# 
tag = "etappi_a2Po-f2Po-a2f2-f2f2_opposite-sign"
# 
# # # # # # # # # # # # # # # # # # # # 

settings_file = fitsfolder(tag, "fit-results.toml")
! isfile(settings_file) && error("no file")
parsed = TOML.parsefile(settings_file)
@unpack settings, fit_results = parsed

# 
setsystem!(Symbol(settings["system"]))
description = Dict()
xlab = ""

if settings["system"] == "compass_ηπ"
    description = description_ηπ
    xlab = L"m_{\eta\pi}\,\,(\mathrm{GeV})"
elseif settings["system"] == "compass_η′π"
    description = description_η′π
    xlab = L"m_{\eta'\pi}\,\,(\mathrm{GeV})"
else
    error("unknown system $(settings["system"])")
end


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


PWs = read_PW_matrices(fitsfolder(tag,"PWs.txt"), used_LMs)
cPWs = read_PW_matrices(fitsfolder(tag,"cPWs.txt"), used_LMs)
cPWs_f = read_PW_matrices(fitsfolder(tag,"cPWs_forward.txt"), used_LMs)
cPWs_b = read_PW_matrices(fitsfolder(tag,"cPWs_backward.txt"), used_LMs)
#

function ambiguousPWs(PWs::TwoBodyPartialWaves{N, V} where N where V <: Number)
    L1_indices = used_LMs..2 .== 1
    ba = bartlettambiguities(PWs.PWs[L1_indices])
    return [TwoBodyPartialWaves(Vector(PWs.LMs),
        [L1_indices[i] ? a[sum(L1_indices[1:i])] : PWs.PWs[i] for i in 1:length(used_LMs)])
        for a in ba]
    #
    
end
reorder(vecofvec) = [getindex.(vecofvec,i) for i in 1:length(vecofvec[1])]

cPWs_all = reorder(ambiguousPWs.(cPWs_f))

let
    N = 3
    M = div(length(used_LMs)-1,3)+1
    plot(layout=grid(M,N), size=(300*N,300*M), frame=:box, grid=false, ylim=(0,:auto),
        xlab=xlab, ylab=L"\mathrm{Intensity}\,/\,40\,\mathrm{MeV}")
    # 
    for i in 1:M*N
        if i > length(used_LMs) 
            plot!(sp=i, axis=false, ticks=false, xlab="", ylab="")
            continue
        end
        # 
        plot!.(sp=i, Ref(pull_mpoints), cPWs_all, :I, i, lab="", l=(1, :gray, 0.3))
    #     #
        plot!(sp=i, pull_mpoints, PWs, :I, i, lab=i!=1 ? "" : "PW projection", l=(2,:red,:dash))
    #     # plot!(sp=i, pull_mpoints, cPWs_f, :I, i, lab=i!=1 ? "" : "cPW projection", l=(2,:red))
    #     # 
        scatter!(sp=i, plotdata.x, plotdata.Iϕ, :I, i,
            xerr=(plotdata.x[2]-plotdata.x[1])/2,
            c=:black, ms=3,
            lab=i!=1 ? "" : "data")
    end
    plot!()
end
savefig(plotsdir(tag, "intensities_with_bartlett.pdf"))

let
    N = 3
    M = div(length(used_LMs)-1,3)+1
    plot(layout=grid(M,N), size=(300*N,300*M), frame=:box, grid=false,
        xlab=xlab, ylab=L"\mathrm{Phase}\,\,(\mathrm{deg})")
    # 
    for i in 1:M*N
        if i == 2 || i > length(used_LMs) 
            plot!(sp=i, axis=false, ticks=false, xlab="", ylab="")
            continue
        end
        # 
        plot!.(sp=i, Ref(pull_mpoints), cPWs_all, :ϕ, i, lab="", l=(1, :gray, 0.3))
        #
        plot!(sp=i, pull_mpoints, PWs, :ϕ, i, lab=i!=1 ? "" : "PW projection", l=(2,:red,:dash))
        # plot!(sp=i, pull_mpoints, cPWs, :ϕ, i, lab=i!=1 ? "" : "cPW projection", l=(2,:red))
        #
        scatter!(sp=i, plotdata.x, plotdata.Iϕ, :ϕ, i,
            xerr=(plotdata.x[2]-plotdata.x[1])/2,
            c=:black, ms=3,
            lab=i!=1 ? "" : "data",)
        #
    end
    plot!()
end
savefig(plotsdir(tag, "phases_with_bartlett.pdf"))


wavesbinsmatrix(func, PWs::Vector{TwoBodyPartialWaves{N, NamedTuple{(:I,:ϕ),V}}} where N where V) = 
    hcat([func(PWs, i) for i in 1:length(used_LMs)]...)
#

writedlm(fitsfolder(tag,"data_ajusted.txt"),
    hcat(plotdata.x,
        wavesbinsmatrix(intensity, plotdata.Iϕ)..:val,
        wavesbinsmatrix(intensity, plotdata.Iϕ)..:err,
        wavesbinsmatrix(phase, plotdata.Iϕ)..:val,
        wavesbinsmatrix(phase, plotdata.Iϕ)..:err))
#
writedlm(fitsfolder(tag,"pw_ajusted.txt"),
    hcat(pull_mpoints,
        wavesbinsmatrix(intensity, changerepresentation.(PWs)),
        wavesbinsmatrix(intensity, changerepresentation.(cPWs)),
        wavesbinsmatrix(intensity, changerepresentation.(cPWs_f)),
        wavesbinsmatrix(intensity, changerepresentation.(cPWs_b))))
