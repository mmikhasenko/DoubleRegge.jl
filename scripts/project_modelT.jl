using DrWatson
@quickactivate "DoubleRegge"

using TOML
using Printf
using DelimitedFiles
using Plots
import Plots.PlotMeasures.mm
theme(:wong2; size=(700, 350), bottom_margin=5mm)

using DoubleRegge

# ─── Config ──────────────────────────────────────────────────────────────────
#
# Load modelT exactly the same way as `scripts/plot_modelT_fb_asymmetry.jl`,
# so that this projection script and the asymmetry script use a consistent
# model definition. Environment variables can override the default mass grid
# and t2 if desired.

settings_file = get(ENV, "DR_SETTINGS", "data/exp_pro/my_model/model-settings.toml")
isfile(settings_file) || error(
    "Settings file not found: $settings_file\nSet DR_SETTINGS env var to override.",
)
parsed = TOML.parsefile(settings_file)

const config = load_modelT_config(parsed)
const settings = config.settings
const tag = settings["tag"]
const reaction_system = config.reaction_system
const LMs = reaction_system.LMs

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

# ─── Mass grid ───────────────────────────────────────────────────────────────

channel_label(system) = system.name == :compass_η′π ? "η′π" : "ηπ"
mass_threshold(system) =
    let ch = system.channel
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

const masses = collect(range(m_min, m_max; length=m_points))

# ─── Amplitude on (cosθ, ϕ) at fixed m ──────────────────────────────────────
#
# We project the amplitude weighted by √q(m), so that |a_{LM}(m)|² directly
# gives the intensity contribution of the (L,M) partial wave (matching the
# convention used in `scripts/project_model.jl`).

fixed_model(m, cosθ, ϕ) = amplitude(model, m, cosθ, ϕ)
fixed_model_sqrtq(m, cosθ, ϕ) = fixed_model(m, cosθ, ϕ) * sqrt(q(m, reaction_system))
intensity_density(m, cosθ, ϕ) = abs2(fixed_model_sqrtq(m, cosθ, ϕ))

model_integral(m) = integrate_dcosθdϕ((cosθ, ϕ) -> intensity_density(m, cosθ, ϕ))[1]

function pw_project_at(m::Float64, L::Int, M::Int)
    return pw_project((cosθ, ϕ) -> fixed_model_sqrtq(m, cosθ, ϕ), L, M)
end

# ─── Run projections ─────────────────────────────────────────────────────────

@info "Built modelT" tag = tag system = reaction_system.name t2 = fixed_t2 npars = length(model.pars)
@info "Projecting onto partial waves" m_min = m_min m_max = m_max m_points = m_points LMs = LMs

const pw_amplitudes = Matrix{ComplexF64}(undef, length(masses), length(LMs))
for (i, m) in enumerate(masses)
    for (j, (L, M)) in enumerate(LMs)
        pw_amplitudes[i, j] = pw_project_at(m, L, M)
    end
    @info "projected" i m
end

const pw_intensities = abs2.(pw_amplitudes)        # |a_{LM}(m)|²
const pw_total = sum(pw_intensities; dims=2) |> vec  # Σ_{LM} |a_{LM}|²
const full_total = model_integral.(masses)           # ∫ |A|² dcosθ dϕ

# fraction carried by the *expanded* (L,M) basis; remainder lives in L>Lmax.
const fraction_in_basis = pw_total ./ full_total

# ─── Save numeric results ────────────────────────────────────────────────────

saveto(filename) = joinpath("data", "exp_pro", tag, filename)
lm_label(L, M) = "$(L)$(M)"

# CSV with: m, Re/Im for each LM, |a|² for each LM, sum(|a|²), full integral
let
    header = String["m"]
    for (L, M) in LMs
        push!(header, "Re_a_$(lm_label(L, M))")
        push!(header, "Im_a_$(lm_label(L, M))")
    end
    for (L, M) in LMs
        push!(header, "I_$(lm_label(L, M))")
    end
    push!(header, "I_sum_basis")
    push!(header, "I_total")

    rows = Matrix{Float64}(undef, length(masses), length(header))
    rows[:, 1] = masses
    col = 2
    for j in 1:length(LMs)
        rows[:, col] = real.(pw_amplitudes[:, j])
        col += 1
        rows[:, col] = imag.(pw_amplitudes[:, j])
        col += 1
    end
    for j in 1:length(LMs)
        rows[:, col] = pw_intensities[:, j]
        col += 1
    end
    rows[:, col] = pw_total
    col += 1
    rows[:, col] = full_total

    csv_path = saveto("pw-projections_$(tag).csv")
    open(csv_path, "w") do io
        println(io, join(header, ","))
        writedlm(io, rows, ',')
    end
    @info "Wrote" csv_path
end

# ─── Plots ───────────────────────────────────────────────────────────────────

const ch_label = channel_label(reaction_system)

# Per-wave intensity panels
let
    n = length(LMs)
    ncols = min(3, n)
    nrows = ceil(Int, n / ncols)
    plot(layout=grid(nrows, ncols), size=(300 * ncols, 260 * nrows))
    for (i, (L, M)) in enumerate(LMs)
        plot!(
            sp=i,
            masses,
            pw_intensities[:, i],
            lw=2,
            lab="",
            title="LM=$(L)$(M)",
            xlab=i > n - ncols ? "m($ch_label) (GeV)" : "",
            ylab="|a_{LM}|²",
        )
    end
    savefig(saveto("pw-projections_$(tag).pdf"))
end

# Total intensity: full integral vs sum over expanded basis
let
    plot(
        size=(700, 350),
        xlab="m($ch_label) (GeV)",
        ylab="intensity",
        title="modelT total intensity — $tag",
    )
    plot!(masses, full_total, lw=2, lab="full ∫|A|² dcosθ dϕ", seriescolor=1)
    plot!(masses, pw_total, lw=2, ls=:dash, lab="Σ_{LM} |a_{LM}|² (expanded basis)", seriescolor=2)
    savefig(saveto("pw-total_$(tag).pdf"))
end

# Odd/even/higher-wave fractions (mirrors scripts/project_model.jl)
let
    filtodd = isodd.(getindex.(LMs, 1))
    fOdd = vec(sum(pw_intensities[:, filtodd]; dims=2)) ./ full_total
    fEven = vec(sum(pw_intensities[:, .!filtodd]; dims=2)) ./ full_total
    fHigher = 1 .- fraction_in_basis

    plot(
        size=(700, 350),
        xlab="m($ch_label) (GeV)",
        ylab="fraction",
        title="modelT wave-group fractions — $tag",
        ylim=(0, 1),
        leg=:outerright,
    )
    plot!(masses, fEven, lw=2, lab="Even L (expanded)")
    plot!(masses, fOdd, lw=2, lab="Odd L (expanded)")
    plot!(masses, fHigher, lw=2, lab="L outside expanded basis")
    savefig(saveto("pw-fractions_$(tag).pdf"))
end

@info "Done. Outputs in data/exp_pro/$(tag)/"
