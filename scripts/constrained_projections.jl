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
tag = "etappi_a2Po-f2Po-a2f2-f2f2_opposite-sign"
# 
# # # # # # # # # # # # # # # # # # # # 

settings_file = fitsfolder(tag, "fit-results.toml")
! isfile(settings_file) && error("no file")
parsed = TOML.parsefile(settings_file)
@unpack settings, fit_results = parsed

# 
setsystem!(Symbol(settings["system"]))


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


#  _|              _|                                    
#      _|_|_|    _|_|_|_|    _|_|    _|_|_|      _|_|_|  
#  _|  _|    _|    _|      _|_|_|_|  _|    _|  _|_|      
#  _|  _|    _|    _|      _|        _|    _|      _|_|  
#  _|  _|    _|      _|_|    _|_|_|  _|    _|  _|_|_|    

pw_projections = pw_project_fixed_model.(plotdata.x)
writedlm(fitsfolder(tag,"PWs.txt"), hcat(unfold.(pw_projections)...))

@time cPWs_starting_from_pw = 
    [constrained_pw_projection_fixed_model(x,a) for (x,a) in zip(plotdata.x, pw_projections)]
writedlm(fitsfolder(tag,"cPWs.txt"),
    hcat(unfold.(getproperty.(cPWs_starting_from_pw, :pars))...))
#

                                                                   
#            _|              _|      _|      _|                      
#  _|_|_|    _|    _|_|    _|_|_|_|_|_|_|_|      _|_|_|      _|_|_|  
#  _|    _|  _|  _|    _|    _|      _|      _|  _|    _|  _|    _|  
#  _|    _|  _|  _|    _|    _|      _|      _|  _|    _|  _|    _|  
#  _|_|_|    _|    _|_|        _|_|    _|_|  _|  _|    _|    _|_|_|  
#  _|                                                            _|  
#  _|                                                        _|_|    

v = readdlm(fitsfolder(tag,"cPWs.txt"))
#     hcat(unfold.(getproperty.(cPWs_starting_from_pw, :pars))...))
cPWs_starting_from_pw = [fold(v[i,:]) for i in 1:size(v,1)]

# 
pw_intensities = map(x->abs2.(x), pw_projections)
cpw_intensities = map(x->abs2.(x.pars), cPWs_starting_from_pw)
cpw_intensities = map(x->abs2.(x), [fold(v[i,:]) for i in 1:size(v,1)])



# getindex.((plotdata.Iϕ..:PWs), :I)
let
    plot(layout=grid(3,3), size=(900,900))
    for (i,(L,M)) in enumerate(used_LMs)
        scatter!(sp=i, plotdata.x, [p.Iϕ.PWs[i].I for p in plotdata],
            xerr=(plotdata.x[2]-plotdata.x[1])/2,
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
    for (i,(L,M)) in enumerate(used_LMs)
        scatter!(sp=i, plotdata.x, [p.Iϕ.PWs[i].ϕ for p in plotdata],
            xerr=(plotdata.x[2]-plotdata.x[1])/2,
            c=:black, title="LM=$L$M", ms=3,
            lab=i!=1 ? "" : "data",)
        #
        plot!(sp=i, plotdata.x, getindex.(pw_phases,i) .+ shift(L,M), lab=i!=1 ? "" : "PW projection", l=(2))
        plot!(sp=i, plotdata.x, getindex.(cpw_phases,i), lab=i!=1 ? "" : "cPW projection", l=(2))
        # vspan!(sp=i, fitdata.x[[1,end]], lab="", α=0.1, seriescolor=7)
    end
    plot!(xlab="m(ηπ) (GeV)")
end
