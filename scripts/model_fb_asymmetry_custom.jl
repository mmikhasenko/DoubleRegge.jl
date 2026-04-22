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

const settings = parsed["settings"]
const tag = settings["tag"]
const scalar_α = Float64(settings["scalar_α"])
const s2shift = Float64(get(settings, "s2shift", 0.0))
const reaction_system = getproperty(DoubleRegge, Symbol(settings["system"]))
const fixed_t2 = parse(Float64, get(ENV, "DR_T2", string(get(settings, "t2", -0.1))))

mkpath(joinpath("data", "exp_pro", tag))

trajs = Dict(
    k => trajectory(Float64(v["slope"]), Float64(v["intercept"]))
    for (k, v) in parsed["trajectories"]
)

verts = Dict(
    k => Vertex(trajs[v["trajectory"]], Float64(v["b"]), Float64(v["tau"]))
    for (k, v) in parsed["vertices"]
)

const exchanges = ReggeExchange[
    ReggeExchange(
        verts[ex["top"]],
        verts[ex["bot"]],
        Bool(ex["eta_forward"]),
        ex["label"],
    )
    for ex in parsed["exchanges"]
]

const pars = Float64.(parsed["fit_results"]["fit_minimizer"])
const model = DoubleReggeModel(
    exchanges,
    fixed_t2,
    scalar_α,
    reaction_system,
    pars;
    s2shift = s2shift,
)

channel_label(system) = system.name == :compass_η′π ? "η′π" : "ηπ"
mass_threshold(system) = let ch = system.channel
    ch.m1 + ch.m2
end

const m_min = parse(Float64, get(ENV, "DR_M_MIN", string(mass_threshold(reaction_system) + 0.01)))
const m_max = parse(Float64, get(ENV, "DR_M_MAX", "3.0"))
const m_points = parse(Int, get(ENV, "DR_M_POINTS", "81"))

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
@info "Built model" tag = tag system = reaction_system.name t2 = fixed_t2 npars = length(pars)
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
