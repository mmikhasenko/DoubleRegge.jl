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
using LinearAlgebra
# 
using Measurements


#                                _|                      
#  _|  _|_|    _|_|      _|_|_|      _|_|_|      _|_|    
#  _|_|      _|_|_|_|  _|        _|  _|    _|  _|_|_|_|  
#  _|        _|        _|        _|  _|    _|  _|        
#  _|          _|_|_|    _|_|_|  _|  _|_|_|      _|_|_|  
#                                    _|                  
#                                    _|                  

function intensity(mass_PWs::Vector{TwoBodyPartialWaveIϕs{N,V}} where {N,V}, i::Integer)
    (((mass_PWs)..:PWs)..i)..:I
end
function phase(mass_PWs::Vector{TwoBodyPartialWaveIϕs{N,V}} where {N,V}, i::Integer)
    phases = alignperiodicsequence(((mass_PWs..:PWs)..i)..:ϕ)
    phases_adj = meanshiftbyperiod(phases)
end

@recipe function f(x, mass_PWs::Vector{TwoBodyPartialWaveAs{N,T}} where {N,T},
    what::Symbol, i::Integer; iref=2)
    mass_PWsIϕ = changerepresentation.(mass_PWs; iref=iref)
    (x, mass_PWsIϕ, what, i)
end

@recipe function f(x, mass_PWs::Vector{TwoBodyPartialWaveIϕs{N,V}} where {N,V},
    what::Symbol, i::Integer)
    y = []
    label --> ""
    L,M = mass_PWs[1].LMs[i]
    title --> "LM=$L$M"
    if what == :I
        intensities = intensity(mass_PWs, i)
        y = intensities
    end
    if what == :ϕ
        phases_adj = phase(mass_PWs, i)
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
# 
pull_mpoints = morefrequentrange(plotdata.x, 5)

#  _|              _|                                    
#      _|_|_|    _|_|_|_|    _|_|    _|_|_|      _|_|_|  
#  _|  _|    _|    _|      _|_|_|_|  _|    _|  _|_|      
#  _|  _|    _|    _|      _|        _|    _|      _|_|  
#  _|  _|    _|      _|_|    _|_|_|  _|    _|  _|_|_|    

# 
pw_projections = pw_project_fixed_model.(pull_mpoints)
writedlm(fitsfolder(tag,"PWs.txt"), hcat(unfold.(pw_projections)...))

@time mass_cPWs_starting_from_pw = 
    [constrained_pw_projection_fixed_model(x,a) for (x,a) in zip(pull_mpoints, pw_projections)]
writedlm(fitsfolder(tag,"cPWs.txt"),
    hcat(unfold.(getproperty.(mass_cPWs_starting_from_pw, :pars))...))
#

# ambibuities
@time cpw_forward  = pull_through_forward( pull_mpoints)
writedlm(fitsfolder(tag,"cPWs_forward.txt"),
    hcat(unfold.(getproperty.(cpw_forward, :pars))...))
# 
@time cpw_backward = pull_through_backward(pull_mpoints)
writedlm(fitsfolder(tag,"cPWs_backward.txt"),
    hcat(unfold.(getproperty.(cpw_backward, :pars))...))


#                            _|        _|                      _|    _|      _|                      
#    _|_|_|  _|_|_|  _|_|    _|_|_|          _|_|_|  _|    _|      _|_|_|_|        _|_|      _|_|_|  
#  _|    _|  _|    _|    _|  _|    _|  _|  _|    _|  _|    _|  _|    _|      _|  _|_|_|_|  _|_|      
#  _|    _|  _|    _|    _|  _|    _|  _|  _|    _|  _|    _|  _|    _|      _|  _|            _|_|  
#    _|_|_|  _|    _|    _|  _|_|_|    _|    _|_|_|    _|_|_|  _|      _|_|  _|    _|_|_|  _|_|_|    
#                                                _|                                                  
#                                            _|_|                                                    

function read_PW_matrices(filename, LMs)
    matrix = readdlm(filename)
    amplitudevectors = [fold(matrix[:,i]) for i in 1:size(matrix,2)]
    return TwoBodyPartialWaves.(Ref(LMs), amplitudevectors)
end

mass_PWs = read_PW_matrices(fitsfolder(tag,"PWs.txt"), used_LMs)
mass_cPWs = read_PW_matrices(fitsfolder(tag,"cPWs.txt"), used_LMs)
mass_cPWs_f = read_PW_matrices(fitsfolder(tag,"cPWs_forward.txt"), used_LMs)
mass_cPWs_b = read_PW_matrices(fitsfolder(tag,"cPWs_backward.txt"), used_LMs)
#
# 
function ambiguousPWs(PWs::TwoBodyPartialWaveAs)
    L1_indices = used_LMs..2 .== 1
    ba = bartlettambiguities(PWs.PWs[L1_indices])
    return [TwoBodyPartialWaves(Vector(PWs.LMs),
        [L1_indices[i] ? a[sum(L1_indices[1:i])] : PWs.PWs[i] for i in 1:length(used_LMs)])
        for a in ba]
end

function ambiguity_tracking(init_index, sets)
    N = length(sets)
    indices = [init_index]
    intensities = [reorder([s..:PWs for k in 1:length(used_LMs)]) for s in sets]
    for i in 2:N
        Ii = intensities[i-1][indices[end]]
        set_i = intensities[i]
        v,ind = findmin(map(x->norm(x.-Ii), set_i))
        push!(indices, ind)
    end
    getindex.(sets, indices)
end
# 
pre_PWs = ambiguousPWs.(mass_cPWs_f)
mass_cPWs_all = [ambiguity_tracking(k, pre_PWs) for k in 1:length(pre_PWs[1])]


#                                _|                      _|            
#    _|_|_|  _|_|_|      _|_|_|  _|  _|    _|    _|_|_|        _|_|_|  
#  _|    _|  _|    _|  _|    _|  _|  _|    _|  _|_|      _|  _|_|      
#  _|    _|  _|    _|  _|    _|  _|  _|    _|      _|_|  _|      _|_|  
#    _|_|_|  _|    _|    _|_|_|  _|    _|_|_|  _|_|_|    _|  _|_|_|    
#                                          _|                          
#                                      _|_|                            

function ellh(intensity_cosθϕ, PWs, LMs)
    function integrand(cosθ,ϕ,pars)
        Id = intensity_cosθϕ(cosθ, ϕ)
        Im = abs2(sum(p*Psi(L,M,cosθ,ϕ) for (p,(L,M)) in zip(pars,LMs)))
        Im ≈ 0.0 && (Im=nextfloat(0.0))
        return -Id*log(Im)
    end
    f(pars) = sum(abs2, pars) + integrate_dcosθdϕ((cosθ,ϕ)->integrand(cosθ,ϕ,fold(pars)))[1]
    return f(unfold(PWs))
end

function full_llh(xv, mass_PWs)
    return sum(ellh(
        (cosθ,ϕ)->intensity(x,cosθ,ϕ),
        y.PWs, y.LMs) for (x, y) in zip(xv, mass_PWs))
end

@time llh_ambiguities = full_llh.(Ref(pull_mpoints), mass_cPWs_all)
histogram(llh_ambiguities .- llh_ambiguities[1])
# 

function chi2(Imodel, Idata, i)
    χ2 = 0
    #
    I_m = intensity(Imodel, i)
    I_d = intensity(Idata, i)
    χ2 += sum(((I_d..:val) .- I_m).^2 ./ (I_d..:err).^2)
    (i==2 || i==6) && return χ2
    # phase
    ϕ_d = phase(Idata, i)
    phases = alignperiodicsequence(((Imodel..:PWs)..i)..:ϕ)
    ϕ_m = meanshiftbyperiod(phases, mean(ϕ_d).val)
    χ2 += sum(((ϕ_d..:val) .- ϕ_m).^2 ./ (ϕ_d..:err).^2)
    return χ2
end
full_chi2(Imodel, Idata) = sum(chi2(Imodel, Idata, i) for i in 1:length(used_LMs))

@time chi2_all = [full_chi2(changerepresentation.(expansion[1:5:end]; iref=2), plotdata.Iϕ)
    for expansion in mass_cPWs_all]
#
histogram(chi2_all, bins=100)
_, best_ambiguity_i = findmin(chi2_all)


#            _|              _|      _|      _|                      
#  _|_|_|    _|    _|_|    _|_|_|_|_|_|_|_|      _|_|_|      _|_|_|  
#  _|    _|  _|  _|    _|    _|      _|      _|  _|    _|  _|    _|  
#  _|    _|  _|  _|    _|    _|      _|      _|  _|    _|  _|    _|  
#  _|_|_|    _|    _|_|        _|_|    _|_|  _|  _|    _|    _|_|_|  
#  _|                                                            _|  
#  _|                                                        _|_|    


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
        # plot!.(sp=i, Ref(pull_mpoints), mass_cPWs_all, :I, i, lab="", l=(1, :gray, 0.3))
        plot!(sp=i, pull_mpoints, mass_cPWs_all[best_ambiguity_i], :I, i, lab=i!=1 ? "" : "cPW projection", l=(2, :blue))
    #     #
        plot!(sp=i, pull_mpoints, mass_PWs, :I, i, lab=i!=1 ? "" : "PW projection", l=(2,:red,:dash))
    #     # plot!(sp=i, pull_mpoints, mass_cPWs_f, :I, i, lab=i!=1 ? "" : "cPW projection", l=(2,:red))
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
        # plot!.(sp=i, Ref(pull_mpoints), mass_cPWs_all, :ϕ, i, lab="", l=(1, :gray, 0.3))
        plot!(sp=i, pull_mpoints, mass_cPWs_all[best_ambiguity_i], :ϕ, i, lab=i!=1 ? "" : "cPW projection", l=(2, :blue))
        #
        plot!(sp=i, pull_mpoints, mass_PWs, :ϕ, i, lab=i!=1 ? "" : "PW projection", l=(2,:red,:dash))
        # plot!(sp=i, pull_mpoints, mass_cPWs, :ϕ, i, lab=i!=1 ? "" : "cPW projection", l=(2,:red))
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

wavesbinsmatrix(func, PWs::Vector{TwoBodyPartialWaveIϕs{N,V}} where {N,V}) = 
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
        wavesbinsmatrix(intensity, changerepresentation.(mass_PWs; iref=2)),
        wavesbinsmatrix(intensity, changerepresentation.(mass_cPWs_all[best_ambiguity_i]; iref=2))))


sumintensityoverbins(mass_PWs::Vector{TwoBodyPartialWaveIϕs{N,V}} where {N,V}) = 
    sum(sum, intensity(mass_PWs, i) for i in 1:length(used_LMs))

Itot_PWs = sumintensityoverbins(changerepresentation.(mass_PWs[1:5:end]; iref=2))
Itot_cPWs = sumintensityoverbins(changerepresentation.(mass_cPWs[1:5:end]; iref=2))
Itot_data = sumintensityoverbins(plotdata.Iϕ)

to_toml(v::Measurement) = [v.val, v.err]
# 
open(fitsfolder(tag, "constrained.toml"),"w") do io
    TOML.print(to_toml, io, Dict(
            "settings"=>settings,
            "fit_results"=>fit_results,
            "constrained_fits" => Dict(
                "Itot_PWs"  => Itot_PWs,
                "Itot_cPWs" => Itot_cPWs,
                "Itot_data" => to_toml(Itot_data),
                "pw_fraction" => Itot_PWs/Itot_cPWs)
            ))
end