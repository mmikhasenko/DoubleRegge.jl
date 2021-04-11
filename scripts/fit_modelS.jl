# using DrWatson
# @quickactivate "DoubleRegge"
using Pkg
Pkg.activate(joinpath(@__DIR__,".."))

using Plots
theme(:wong; size=(500,350))

using QuadGK
using Cuba
# 
using Optim
using TypedTables
using ForwardDiff
#
using TOML
using StaticArrays
using Measurements
using DelimitedFiles
# 
using DoubleRegge

# # # # # # # # # # # # # # # # # # # # 
# 
tag = "etappi_a2Po-a2f2-f2f2-PoPo_opposite-sign"
# 
# # # # # # # # # # # # # # # # # # # # 

# 
settings_file = fitsfolder(tag, "settings.toml")
! isfile(settings_file) && error("no file")
settings = TOML.parsefile(settings_file)

# 
setsystem!(Symbol(settings["system"]))
description = settings["system"] == "compass_ηπ" ? description_ηπ :
                (settings["system"] == "compass_η′π" ? description_η′π :
                    error("unknown system $(settings["system"])"))

# get data
data = read_data(settings["pathtodata"], description)
fitdata = filter(data) do x
    inlims(x.x, settings["fitrange"])
end

# build model
const model = build_model(
    sixexchages[settings["exchanges"]],
    settings["t2"],
    settings["scale_α"])
#
# ellh fit functions
function integrand(cosθ,ϕ,pars)
    Id = abs2.(recamp.(cosθ, ϕ, fitdata.amps))
    Am = model.(fitdata.x, cosθ, ϕ; pars=pars)
    Im = abs2.(Am) .* q.(fitdata.x)
    # @show pars
    return sum(Im .- Id .* log.(Im))
end
ellh(pars) = integrate_dcosθdϕ((cosθ,ϕ)->integrand(cosθ,ϕ,pars))[1]
# ellh(settings["initial_pars"])

integrand′(cosθ,ϕ,pars) = ForwardDiff.gradient(p->integrand(cosθ,ϕ,p), pars)
ellh′(pars) = integrate_dcosθdϕ((cosθ,ϕ)->integrand′(cosθ,ϕ,pars), dims=length(pars))
ellh′!(stor,pars) = (stor .= ellh′(pars)) 
# ellh′!, 
# fit
ft = Optim.optimize(ellh, settings["initial_pars"], BFGS(),
               Optim.Options(show_trace = true, g_tol=1e-4, iterations=15))
#

# save results
fit_results = Dict(
    "fit_converged" => ft.x_converged,
    "fit_minimizer" => ft.minimizer,
    "fit_minimum" => ft.minimum)
    
output_name = fitsfolder(tag, "fit-results.toml");

open(output_name,"w") do io
    TOML.print(io, Dict(
            "settings"=>settings,
            "fit_results"=>fit_results))
end

