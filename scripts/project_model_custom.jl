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

const config = load_modelT_config(parsed)
const settings   = config.settings
const tag        = settings["tag"]
const reaction_system = config.reaction_system
const pathtoevents = settings["pathtoevents"]
const treename   = get(settings, "treename", "ProgramVariables")

mkpath(joinpath("data", "exp_pro", tag))

const model = config.model

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
        TEventKinematics(
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
    binned = [TEventKinematics[] for _ in 1:(length(bin_edges)-1)]
    for ev in events
        m = sqrt(ev.s12)
        idx = searchsortedlast(bin_edges, m)
        if 1 <= idx < length(bin_edges)
            push!(binned[idx], ev)
        end
    end
    return binned
end

@info "Binning events into $(length(bin_centers)) mass bins …"
const binned = bin_events(all_events, bin_edges)

# ─── MC partial wave projection ──────────────────────────────────────────────

function pw_project_mc(model, events, LMs)
    isempty(events) && return fill(0.0 + 0.0im, length(LMs))
    amps = amplitude.(Ref(model), events)
    return [
        mean(i -> amps[i] * conj(Psi(L, M, events[i].cosθ, events[i].ϕ)), eachindex(events))
        for (L, M) in LMs
    ]
end

const LMs = reaction_system.LMs

@info "Computing partial wave projections …"
const pw_amplitudes = [
    pw_project_mc(model, binned[b], LMs)
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
