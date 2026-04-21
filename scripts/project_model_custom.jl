using DrWatson
@quickactivate "DoubleRegge"
using TOML
using UnROOT
using Plots
import Plots.PlotMeasures.mm
theme(:wong2; size = (500, 350), bottom_margin = 5mm)
using LaTeXStrings
using DoubleRegge
using Statistics

# ─── TOML ────────────────────────────────────────────────────────────────────

settings_file = get(ENV, "DR_SETTINGS", "data/exp_pro/my_model/model-settings.toml")
isfile(settings_file) || error("Settings file not found: $settings_file\nSet DR_SETTINGS env var to override.")
parsed = TOML.parsefile(settings_file)

const settings   = parsed["settings"]
const tag        = settings["tag"]
const scalar_α   = Float64(settings["scalar_α"])
const s2shift    = Float64(get(settings, "s2shift", 0.0))
const pathtoevents = settings["pathtoevents"]
const treename   = get(settings, "treename", "ProgramVariables")

mkpath(joinpath("data", "exp_pro", tag))

# ─── Build model from TOML ───────────────────────────────────────────────────

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
const model = DoubleReggeModel(exchanges, 0.0, scalar_α, compass_ηπ, pars; s2shift = s2shift)

# ─── Load events from ROOT ───────────────────────────────────────────────────

function load_events(path::String, treename::String)
    f = ROOTFile(path)
    ev = treename * "/Event/"
    rd(b)    = UnROOT.array(f, ev * b)
    rdvar(b) = UnROOT.array(f, ev * "var/var." * b)
    cosθ_v = rd("CosThetaGJ")
    ϕ_v    = rd("PhiGJ")
    s12_v  = rdvar("s12")
    s13_v  = rdvar("s13")
    s23_v  = rdvar("s23")
    t1_v   = rdvar("t1_1")
    tπ_v   = rdvar("t1_2")
    t2_v   = rdvar("t2")
    s_v    = rdvar("s")
    K_v    = rdvar("KFactor")
    return [
        EventKinematics(
            s_v[i], s12_v[i], s13_v[i], s23_v[i],
            t1_v[i], tπ_v[i], t2_v[i],
            cosθ_v[i], ϕ_v[i], K_v[i],
        )
        for i in eachindex(s12_v)
    ]
end

@info "Loading events from $pathtoevents …"
const all_events = load_events(pathtoevents, treename)
@info "Loaded $(length(all_events)) events"

# ─── Mass binning ────────────────────────────────────────────────────────────

const ms_all = sqrt.(getfield.(all_events, :s12))
const bin_step = Float64(get(settings, "bin_step", 0.1))
const bin_edges = collect(range(minimum(ms_all), maximum(ms_all); step = bin_step))
const bin_centers = (bin_edges[1:end-1] .+ bin_edges[2:end]) ./ 2

function bin_events(events, bin_edges)
    return [
        filter(ev -> bin_edges[i] ≤ sqrt(ev.s12) < bin_edges[i+1], events)
        for i in 1:length(bin_edges)-1
    ]
end

@info "Binning events into $(length(bin_centers)) mass bins …"
const binned = bin_events(all_events, bin_edges)

# ─── MC partial wave projection ──────────────────────────────────────────────

function pw_project_mc(model, events, L, M)
    isempty(events) && return 0.0 + 0.0im
    return mean(ev -> amplitude(model, ev) * conj(Psi(L, M, ev.cosθ, ev.ϕ)), events)
end

const LMs = compass_ηπ.LMs  # [(1,1),(2,1),(2,2),(3,1),(4,1),(5,1),(6,1)]

@info "Computing partial wave projections …"
const pw_amplitudes = [
    [pw_project_mc(model, binned[b], L, M) for (L, M) in LMs]
    for b in eachindex(bin_centers)
]
const pw_intensities = [abs2.(cs) for cs in pw_amplitudes]
const total_intensity = sum.(pw_intensities)

# ─── Plots ───────────────────────────────────────────────────────────────────

function saveto(filename)
    joinpath("data", "exp_pro", tag, filename)
end

# Partial wave projections
let
    n = length(LMs)
    ncols = 3
    nrows = ceil(Int, n / ncols)
    plot(layout = grid(nrows, ncols), size = (900, 300 * nrows))
    for (i, (L, M)) in enumerate(LMs)
        plot!(
            sp = i,
            bin_centers,
            getindex.(pw_intensities, i),
            lab = "",
            lw = 2,
            title = "LM=$L$M",
            xlab = i > n - ncols ? "m(ηπ) (GeV)" : "",
            ylab = "intensity",
        )
    end
    savefig(saveto("pw-projections_$tag.pdf"))
end

# Total intensity
let
    plot(
        bin_centers,
        total_intensity,
        lw = 2,
        lab = "",
        xlab = "m(ηπ) (GeV)",
        ylab = "intensity",
        title = tag,
    )
    savefig(saveto("total-intensity_$tag.pdf"))
end

@info "Done. PDFs written to data/exp_pro/$tag/"
