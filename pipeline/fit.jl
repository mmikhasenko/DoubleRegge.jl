using DrWatson
@quickactivate "DoubleRegge"
# 
using Plots
# 
using Optim
using TypedTables
using ForwardDiff
using TOML
# 
using DoubleRegge

push!(ARGS, "pipeline/fit_settings.toml")

length(ARGS) < 1 && error("provide config file as the first argument!")
settings_file = ARGS[1]
!isfile(settings_file) && error("no file")

parsed = TOML.parsefile(settings_file)
settings = parsed["settings"]

const reaction_system = getproperty(DoubleRegge, Symbol(settings["system"]))
const LMs = reaction_system.LMs

# data
data = read_data(joinpath(@__DIR__, settings["pathtodata"]), reaction_system)

amplitudes = [sqrt.(is) .* cis.(ϕs) for (is, ϕs) in zip(data.I, data.ϕ)]

# range
fitrangemap = map(x -> inlims(x.x, settings["fitrange"]), data)
fitdata = Table(data[fitrangemap], amps = amplitudes[fitrangemap])

# fit
const exchanges = six_exchanges[settings["exchanges"]]
const model = DoubleReggeModel(
    exchanges,
    settings["t2"],
    settings["scale_α"],
    reaction_system,
    settings["initial_pars"];
    s2shift = get(settings, "s2_shift", 0.0),
)
intensity(m, cosθ, ϕ; pars) = abs2(amplitude(with_parameters(model, pars), m, cosθ, ϕ)) * q(m, reaction_system)

function integrand(cosθ, ϕ, pars)
    trial_model = with_parameters(model, pars)
    Id = abs2.(recamp.(cosθ, ϕ, fitdata.amps, Ref(LMs)))
    Am = amplitude.(Ref(trial_model), fitdata.x, cosθ, ϕ)
    Im = abs2.(Am) .* q.(fitdata.x, Ref(reaction_system))
    return sum(Im .- Id .* log.(Im))
end
ellh(pars) = integrate_dcosθdϕ((cosθ, ϕ) -> integrand(cosθ, ϕ, pars))[1]

# 
integrand′(cosθ, ϕ, pars) = ForwardDiff.gradient(p -> integrand(cosθ, ϕ, p), pars)
ellh′(pars) = integrate_dcosθdϕ((cosθ, ϕ) -> integrand′(cosθ, ϕ, pars), dims = length(pars))
ellh′!(stor, pars) = (stor .= ellh′(pars))

ft = Optim.optimize(ellh, ellh′!, settings["initial_pars"], BFGS(),
    Optim.Options(show_trace = true, g_tol = 1e-4, iterations = 15))
#
fit_results = Dict(
    "fit_converged" => ft.x_converged,
    "fit_minimizer" => ft.minimizer,
    "fit_minimum" => ft.minimum)

output_name = joinpath("fit-results.toml");

open(output_name, "w") do io
    TOML.print(io, Dict(
        "settings" => settings,
        "fit_results" => fit_results))
end
