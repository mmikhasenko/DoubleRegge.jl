using DrWatson
@macroexpand(@quickactivate "DoubleRegge")
# 
using Plots
# 
using Optim
using TypedTables
using ForwardDiff
using TOML
# 
using DoubleRegge

length(ARGS) < 1 && error("provide config file as the first argument!")

settings_file = ARGS[1]
! isfile(settings_file) && error("no file")

parsed = TOML.parsefile(settings_file)
settings = parsed["settings"]

setsystem!(Symbol(settings["system"]))

# data
const LMs = compass_ηπ_LMs
data = Table(x_IδI_ϕδϕ_compass_ηπ(settings["pathtodata"]))
amplitudes = [sqrt.(is) .* cis.(ϕs) for (is,ϕs) in zip(data.I, data.ϕ)]

# range
fitrangemap = map(x->inlims(x.x, settings["fitrange"]), data)
fitdata = Table(data[fitrangemap], amps = amplitudes[fitrangemap])

# fit
const exchanges = sixexchages[settings["exchanges"]]
const model = build_model(exchanges, settings["t2"], settings["scale_α"]; s2shift=settings["s2_shift"])
intensity(m, cosθ, ϕ; pars) = abs2(model(m, cosθ, ϕ; pars=pars))*q(m)

function integrand(cosθ,ϕ,pars)
    Id = abs2.(recamp.(cosθ, ϕ, fitdata.amps, Ref(LMs)))
    Am = model.(fitdata.x, cosθ, ϕ; pars=pars)
    Im = abs2.(Am) .* q.(fitdata.x)
    return sum(Im .- Id .* log.(Im))
end
ellh(pars) = integrate_dcosθdϕ((cosθ,ϕ)->integrand(cosθ,ϕ,pars))[1]

# 
integrand′(cosθ,ϕ,pars) = ForwardDiff.gradient(p->integrand(cosθ,ϕ,p), pars)
ellh′(pars) = integrate_dcosθdϕ((cosθ,ϕ)->integrand′(cosθ,ϕ,pars), dims=length(pars))
ellh′!(stor,pars) = (stor .= ellh′(pars)) 

ft = Optim.optimize(ellh, ellh′!, settings["initial_pars"], BFGS(),
               Optim.Options(show_trace = true, g_tol=1e-4, iterations=15))
#
fit_results = Dict(
    "fit_converged" => ft.x_converged,
    "fit_minimizer" => ft.minimizer,
    "fit_minimum" => ft.minimum)

output_name = joinpath("fit-results.toml");

open(output_name,"w") do io
    TOML.print(io, Dict(
            "settings"=>settings,
            "fit_results"=>fit_results))
end