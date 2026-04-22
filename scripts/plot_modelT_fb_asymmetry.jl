using DrWatson
@quickactivate "DoubleRegge"

using TOML
using Plots
import Plots.PlotMeasures.mm
theme(:wong2; size = (700, 350), bottom_margin = 5mm)

using DoubleRegge

settings_file = get(ENV, "DR_SETTINGS", "data/exp_pro/my_model/model-settings.toml")
isfile(settings_file) || error(
    "Settings file not found: $settings_file\nSet DR_SETTINGS env var to override.",
)
parsed = TOML.parsefile(settings_file)

const config = load_modelT_config(parsed)
const settings = config.settings
const tag = settings["tag"]
const reaction_system = config.reaction_system

mkpath(joinpath("data", "exp_pro", tag))
const model = if haskey(ENV, "DR_T2")
    TDoubleReggeModel(
        config.model.exchanges,
        parse(Float64, ENV["DR_T2"]),
        reaction_system,
        config.model.pars,
    )
else
    config.model
end
const fixed_t2 = model.t2

channel_label(system) = system.name == :compass_η′π ? "η′π" : "ηπ"
mass_threshold(system) = let ch = system.channel
    ch.m1 + ch.m2
end

const default_m_min = Float64(get(settings, "m_min", mass_threshold(reaction_system) + 0.01))
const default_m_max = Float64(get(settings, "m_max", 3.0))
const default_m_points = Int(get(settings, "m_points", 81))

const m_min = parse(Float64, get(ENV, "DR_M_MIN", string(default_m_min)))
const m_max = parse(Float64, get(ENV, "DR_M_MAX", string(default_m_max)))
const m_points = parse(Int, get(ENV, "DR_M_POINTS", string(default_m_points)))

m_min < m_max || error("Need DR_M_MIN < DR_M_MAX, got $m_min >= $m_max")
m_points >= 2 || error("Need at least 2 mass points, got $m_points")

fixed_model(m, cosθ, ϕ) = amplitude(model, m, cosθ, ϕ)
intensity_density(m, cosθ, ϕ) = abs2(fixed_model(m, cosθ, ϕ)) * q(m, reaction_system)

model_integral_forward(m) =
    integrate_dcosθdϕ((cosθ, ϕ) -> intensity_density(m, cosθ, ϕ), (0, 1))[1]
model_integral_backward(m) =
    integrate_dcosθdϕ((cosθ, ϕ) -> intensity_density(m, cosθ, ϕ), (-1, 0))[1]

asymmetry(f, b) = (f - b) / (f + b)

const masses = collect(range(m_min, m_max; length = m_points))
@info "Built model" tag = tag system = reaction_system.name t2 = fixed_t2 npars = length(model.pars)
@info "Evaluating forward/backward asymmetry" m_min = m_min m_max = m_max m_points = m_points

const forward = model_integral_forward.(masses)
const backward = model_integral_backward.(masses)
const fb_asymmetry = [asymmetry(f, b) for (f, b) in zip(forward, backward)]

outfile = joinpath("data", "exp_pro", tag, "fb-asymmetry_$(tag).pdf")
mass_label = channel_label(reaction_system)

let
    plot(
        size = (900, 350),
        layout = grid(1, 2),
        title = ["forward/backward intensity" "forward-backward asymmetry"],
        xlab = ["m($mass_label) (GeV)" "m($mass_label) (GeV)"],
        ylab = ["intensity" "(F-B) / (F+B)"],
    )
    plot!(sp = 1, masses, forward, lw = 2, lab = "forward", seriescolor = 3)
    plot!(sp = 1, masses, backward, lw = 2, lab = "backward", seriescolor = 4)
    plot!(sp = 2, masses, fb_asymmetry, lw = 2, lab = "", ylim = (-1, 1))
    savefig(outfile)
end

@info "Wrote" outfile
