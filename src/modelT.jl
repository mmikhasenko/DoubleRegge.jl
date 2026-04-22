struct TVertex{T}
    α::T
    b::Float64
    τ::Float64
end

form_factor(v::TVertex, t) = exp(v.b * t)

struct TReggeExchange{T1,T2,S<:AbstractString}
    top::TVertex{T1}
    bot::TVertex{T2}
    η_forward::Bool
    label::S
end

# ─── Kinematics ──────────────────────────────────────────────────────────────
#
# Two orthogonal representations of a production-decay event are used:
#
#   KinematicsGJ – Gottfried-Jackson angles:  (s, s1, cosθ, ϕ, t2)
#                  the "model-side" input, parametrising the phase space.
#
#   KinematicsM  – Mandelstam invariants + K: (s, s1, s13, s23, t1, tπ, t2, K)
#                  the minimal set the amplitude actually needs.
#
# KinematicsM can be built from KinematicsGJ (needs a reaction_system), or
# directly from externally provided invariants (e.g. MC events on disk).

struct KinematicsGJ
    s::Float64
    s1::Float64
    cosθ::Float64
    ϕ::Float64
    t2::Float64
end

struct KinematicsM
    s::Float64
    s1::Float64    # s12 in 3-body labeling: (η π) invariant mass squared
    s13::Float64   # sηp
    s23::Float64   # sπp
    t1::Float64    # momentum transfer to η
    tπ::Float64    # momentum transfer to π
    t2::Float64    # momentum transfer to recoil proton
    K::Float64     # kinematic K-factor (carries the sin ϕ dependence)
end

function KinematicsM(gj::KinematicsGJ, reaction_system)
    vars = (s=gj.s, s1=gj.s1, cosθ=gj.cosθ, ϕ=gj.ϕ, t2=gj.t2)
    return KinematicsM(
        gj.s,
        gj.s1,
        sηp(vars, reaction_system),
        sπp(vars, reaction_system),
        t1(vars, reaction_system),
        tπ(vars, reaction_system),
        gj.t2,
        Kfactor(vars, reaction_system),
    )
end

# ─── Model ───────────────────────────────────────────────────────────────────

struct TDoubleReggeModel{E<:AbstractVector,P,R}
    exchanges::E
    t2::Float64
    scalar_α::Float64
    reaction_system::R
    pars::P
end

function TDoubleReggeModel(
    exchanges,
    t2::Real,
    scalar_α::Real,
    reaction_system,
    pars,
)
    length(exchanges) == length(pars) || throw(
        DimensionMismatch("expected $(length(exchanges)) parameters, got $(length(pars))"),
    )
    return TDoubleReggeModel(
        exchanges,
        Float64(t2),
        Float64(scalar_α),
        reaction_system,
        pars,
    )
end

with_parameters(model::TDoubleReggeModel, pars) = TDoubleReggeModel(
    model.exchanges,
    model.t2,
    model.scalar_α,
    model.reaction_system,
    pars,
)

# ─── Amplitude ───────────────────────────────────────────────────────────────

"""
    modelDR(exchange::TReggeExchange, kin::KinematicsM; α′)

Core single-diagram amplitude. Works purely on Mandelstam + K: no reaction
system is needed here. The choice of `(t, s2)` inside a diagram is driven by
`exchange.η_forward`.
"""
function modelDR(exchange::TReggeExchange, kin::KinematicsM; α′::Float64=0.9)
    top, bot = exchange.top, exchange.bot
    t = exchange.η_forward ? kin.t1 : kin.tπ
    s2 = exchange.η_forward ? kin.s23 : kin.s13

    α1 = top.α(t)
    α2 = bot.α(kin.t2)

    prefactor =
        -kin.K * sf_gamma(1 - α1) * sf_gamma(1 - α2) *
        (α′ * kin.s1)^α1 * (α′ * s2)^α2 / (α′ * kin.s)
    ff = form_factor(top, t) * form_factor(bot, kin.t2)

    η = kin.s / (α′ * kin.s1 * s2)
    τ1, τ2 = top.τ, bot.τ
    ξ1 = ξ(τ1, α1)
    ξ21 = ξ(τ2, τ1, α2, α1)
    ξ2 = ξ(τ2, α2)
    ξ12 = ξ(τ1, τ2, α1, α2)
    vertex = η^α1 * ξ1 * ξ21 * V(α1, α2, η) +
             η^α2 * ξ2 * ξ12 * V(α2, α1, η)

    return ff * prefactor * vertex
end

"""
    modelDR(exchange, gj::KinematicsGJ, reaction_system; α′)

Convenience wrapper that converts GJ → Mandelstam on the fly.
"""
modelDR(
    exchange::TReggeExchange,
    gj::KinematicsGJ,
    reaction_system;
    α′::Float64=0.9,
) = modelDR(exchange, KinematicsM(gj, reaction_system); α′=α′)

# ─── Full amplitude (sum over diagrams) ──────────────────────────────────────

function amplitude(model::TDoubleReggeModel, kin::KinematicsM)
    generator = (
        p * modelDR(exchange, kin; α′=model.scalar_α)
        for (p, exchange) in zip(model.pars, model.exchanges)
    )
    return mysum(typeof(1im * model.pars[1]), generator)
end

amplitude(model::TDoubleReggeModel, gj::KinematicsGJ) =
    amplitude(model, KinematicsM(gj, model.reaction_system))

function amplitude(model::TDoubleReggeModel, m, cosθ, ϕ)
    gj = KinematicsGJ(model.reaction_system.s0, m^2, cosθ, ϕ, model.t2)
    return amplitude(model, gj)
end

# ─── TOML loader ─────────────────────────────────────────────────────────────

function _get_tmodel_couplings(parsed, diagrams)
    if haskey(parsed, "couplings")
        couplings = parsed["couplings"]
        return Float64[
            Float64(couplings[get(diagram, "coupling", diagram["label"])]) for diagram in diagrams
        ]
    end
    if haskey(parsed, "fit_results") && haskey(parsed["fit_results"], "fit_minimizer")
        return Float64.(parsed["fit_results"]["fit_minimizer"])
    end
    throw(ArgumentError("expected either [couplings] or [fit_results].fit_minimizer in modelT config"))
end

function load_modelT_config(parsed; settings_key::AbstractString="settings")
    settings = parsed[settings_key]
    reaction_system = getproperty(DoubleRegge, Symbol(settings["system"]))

    trajectories = Dict(
        key => trajectory(Float64(value["slope"]), Float64(value["intercept"]))
        for (key, value) in parsed["trajectories"]
    )

    vertices = Dict(
        key => TVertex(
            trajectories[value["trajectory"]],
            Float64(get(value, "b", 0.0)),
            Float64(get(value, "tau", 1.0)),
        ) for (key, value) in parsed["vertices"]
    )

    diagrams = haskey(parsed, "diagrams") ? parsed["diagrams"] :
               haskey(parsed, "exchanges") ? parsed["exchanges"] :
               throw(ArgumentError("expected [[diagrams]] or [[exchanges]] in modelT config"))

    exchanges = TReggeExchange[
        TReggeExchange(
            vertices[diagram["top"]],
            vertices[diagram["bot"]],
            Bool(diagram["eta_forward"]),
            diagram["label"],
        ) for diagram in diagrams
    ]

    pars = _get_tmodel_couplings(parsed, diagrams)
    length(exchanges) == length(pars) || throw(
        DimensionMismatch("expected $(length(exchanges)) couplings, got $(length(pars))"),
    )

    model = TDoubleReggeModel(
        exchanges,
        Float64(get(settings, "t2", 0.0)),
        Float64(settings["scalar_α"]),
        reaction_system,
        pars,
    )

    return (
        model=model,
        settings=settings,
        reaction_system=reaction_system,
        trajectories=trajectories,
        vertices=vertices,
        diagrams=diagrams,
        couplings=pars,
    )
end

load_modelT_from_toml(path::AbstractString; settings_key::AbstractString="settings") =
    load_modelT_config(TOML.parsefile(path); settings_key=settings_key)
