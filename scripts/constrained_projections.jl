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
setsystem!(:compass_ηπ)


#  _|_|_|  _|_|      _|_|    _|      _|    _|_|    
#  _|    _|    _|  _|    _|  _|      _|  _|_|_|_|  
#  _|    _|    _|  _|    _|    _|  _|    _|        
#  _|    _|    _|    _|_|        _|        _|_|_|  

fold(longp) = (Np = div(length(longp),2); longp[1:Np] .+ 1im .* longp[(Np+1):end])
unfold(p) = vcat(real.(p), imag.(p))
# test
(v = rand(Complex{Float64},10); prod(fold(unfold(v)) .== v))
# 
function constrained_pw_projection(intensity_cosθϕ, init_pars, LMs)
    
    function integrand(cosθ,ϕ,pars)
        Id = intensity_cosθϕ(cosθ, ϕ)
        Im = abs2(sum(p*Psi(L,M,cosθ,ϕ) for (p,(L,M)) in zip(pars,LMs)))
        Im ≈ 0.0 && (Im=nextfloat(0.0))
        return -Id*log(Im)
    end
    f(pars) = sum(abs2, pars) + integrate_dcosθdϕ((cosθ,ϕ)->integrand(cosθ,ϕ,fold(pars)))[1]
    #
    f = Optim.optimize(f, unfold(init_pars), BFGS(), # f′!, 
                Optim.Options(show_trace = true, iterations=100))
    return fold(f.minimizer)
end

#                            _|            
#    _|_|_|    _|_|      _|_|_|    _|_|    
#  _|        _|    _|  _|    _|  _|_|_|_|  
#  _|        _|    _|  _|    _|  _|        
#    _|_|_|    _|_|      _|_|_|    _|_|_|  

settings_file = joinpath("data", "exp_pro","fit-results_bottom-Po_Np=3.toml")
# settings_file = joinpath("data", "exp_pro","fit-results_a2Po-f2f2-PoPo_Np=3.toml")
! isfile(settings_file) && error("no file")

# 
parsed = TOML.parsefile(settings_file)
settings = parsed["settings"]
fit_results = parsed["fit_results"]
# 

# fit
const exchanges = sixexchages[settings["exchanges"]]
const t2 = settings["t2"]
const α = settings["scale_α"]
const model = build_model(exchanges, settings["t2"], settings["scale_α"])
const fixed_pars = fit_results["fit_minimizer"]
fixed_model(m,cosθ,ϕ) = model(m,cosθ,ϕ; pars=fixed_pars)
intensity(m, cosθ, ϕ) = abs2(fixed_model(m, cosθ, ϕ))*q(m)

function pw_project_fixed_model(m)
    amplitude(cosθ,ϕ) = fixed_model(m,cosθ,ϕ)*sqrt(q(m))
    pws = [pw_project(amplitude,L,M) for (L,M) in LMs]
    return pws
end

constrained_pw_projection_fixed_model(m, init_pars) =
    constrained_pw_projection((cosθ,ϕ)->intensity(m,cosθ,ϕ), init_pars, compass_ηπ_LMs)

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


@time cPWs_starting_from_compass_pw =
    [constrained_pw_projection_fixed_model(x,a) for (x,a) in zip(plotdata.x, plotdata.amps)]
writedlm(joinpath("data", "exp_pro", "constrained_pws_$(settings["tag"])_starting_from_compass_pw.txt"),
    hcat(unfold.(cPWs_starting_from_compass_pw)...))
# # 
@time cPWs_starting_from_pw = 
    [constrained_pw_projection_fixed_model(x,a) for (x,a) in zip(plotdata.x, pw_projections)]
writedlm(joinpath("data", "exp_pro", "constrained_pws_$(settings["tag"])_starting_from_pw.txt"),
    hcat(unfold.(cPWs_starting_from_pw)...))
# 
pw_intensities = map(x->abs2.(x), pw_projections)
cpw_intensities = map(x->abs2.(x), cPWs_starting_from_pw)
# # 
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





# # projections
# cpw_projections = map(x->vcat(x...), cPWs)
# # 

# # test
# const testLM = [(L=1, M=1), (L=3,M=1), (L=3,M=2)]
# const testA = copy(rand(Complex{Float64}, length(testLM)))
# test_intensity(cosθ,ϕ) = abs2(recamp(cosθ,ϕ,testA,testLM))
# # prod(constrained_pw_projection(test_intensity, testA, testLM) .≈ testA)
# # hcat(constrained_pw_projection(test_intensity, testA, testLM), testA)



# authomatic derivative
# function f′(pars)
#     dint(cosθ,ϕ) = ForwardDiff.gradient(
#             p->integrand(cosθ,ϕ, fold(p)),
#         unfold(pars))
#     integral = integrate_dcosθdϕ(dint; dims=2*length(pars))
#     fold(integral)
# end
# f′!(stor,pars) = copyto!(stor, f′(pars))


#
# 
# @btime cuhre((x,f)->f[1]=abs2(fixed_model(2.1, 2x[1]-1,π*(2x[2]-1))),2,1)[1] # 18.1 ms
# @btime integrate_dcosθdϕ((cosθ,ϕ)->abs2(fixed_model(2.1, cosθ,ϕ))) # 18.1 ms

# # @profview fixed_model(2.1,0.3,0.3) # 6.3 ms


# let intensity_cosθϕ = (cosθ,ϕ)->abs2(fixed_model(2.1, cosθ,ϕ))
#     function integrand(cosθ,ϕ,pars)
#         Id = intensity_cosθϕ(cosθ, ϕ)
#         # Id = abs2(fixed_model(2.1, cosθ,ϕ))
#         Im = abs2(sum(p*Psi(L,M,cosθ,ϕ) for (p,(L,M)) in zip(pars,LMs)))
#         Im ≈ 0.0 && (Im=nextfloat(0.0))
#         # return sum(abs2, pars)/(4π) - Id * log(Im)
#         return Im - Id * log(Im)
#     end
#     # 
#     f(pars) = integrate_dcosθdϕ((cosθ,ϕ)->integrand(cosθ,ϕ,pars))[1]
#     @btime $f(rand(14))
#     # f(rand(14))
# end


# let intensity_cosθϕ = (cosθ,ϕ)->abs2(fixed_model(2.1, cosθ,ϕ)), init_pars=pw_projections[1]
#     function integrand(cosθ,ϕ,pars)
#         Id = intensity_cosθϕ(cosθ, ϕ)
#         Im = abs2(sum(p*Psi(L,M,cosθ,ϕ) for (p,(L,M)) in zip(pars,LMs)))
#         Im ≈ 0.0 && (Im=nextfloat(0.0))
#         return -Id*log(Im)
#     end
#     #
#     f(pars) = sum(abs2, pars) + integrate_dcosθdϕ((cosθ,ϕ)->integrand(cosθ,ϕ,fold(pars)))[1]
#     # 
#     f = Optim.optimize(f, unfold(init_pars), BFGS(),
#                 Optim.Options(show_trace = true, iterations=50))
#     return f.minimizer
# end

↑
# @time
# constrained_pw_projection_fixed_model(plotdata.x[2], pw_projections[2])
# # 

# hcat(sum.(pw_intensities),sum.(cpw_intensities))


# @show ↑

# const testA = copy(data.amps[44])
# const testLM = copy(LMs)
# test_intensity(cosθ,ϕ) = abs2(recamp(cosθ,ϕ,testA,testLM))
# let intensity_cosθϕ = test_intensity, init_pars=data.amps[44], LMs=testLM
    
#     function integrand(cosθ,ϕ,pars)
#         Id = intensity_cosθϕ(cosθ, ϕ)
#         Im = abs2(sum(p*Psi(L,M,cosθ,ϕ) for (p,(L,M)) in zip(pars,LMs)))
#         return Im - Id * log(Im)
#     end
#     # 
#     f(pars) = integrate_dcosθdϕ((cosθ,ϕ)->integrand(cosθ,ϕ,pars))[1]
#     # 
#     function f′(pars)
#         dint(cosθ,ϕ) = ForwardDiff.gradient(
#                 p->integrand(cosθ,ϕ, fold(p)),
#             unfold(pars))
#         integral = integrate_dcosθdϕ(dint; dims=2*length(pars))
#         fold(integral)
#     end
#     f′!(stor,pars) = copyto!(stor, f′(pars))
#     #
#     f = Optim.optimize(f, f′!, init_pars, BFGS(), #, f′!
#                 Optim.Options(show_trace = true))
#     return f.minimizer
# end

# 600*15
#  = 6000s + 3000s
#  = 100m + 50 m
#  = 3h
# 
# x[s]/240 = [x]h
# 
# 
# let pars0 = rand(7)
#     hcat(dellh(pars0), numerical_grad(ellh,pars0))
# end
#
# function numerical_grad(f,p)
#     e0 = f(p)
#     Np = length(p)
#     ϵ = 1e-7
#     nd = [((f(p .+ ϵ .* [k==i for i in 1:Np]) - e0) / ϵ) for k in 1:Np]
#     return nd
# end
