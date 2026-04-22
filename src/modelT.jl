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

TReggeExchange(α_top, α_bot, η_forward::Bool, label::AbstractString) = TReggeExchange(
    TVertex(α_top, 0.0, 1.0),
    TVertex(α_bot, 0.0, 1.0),
    η_forward,
    label,
)

struct TKinematics
    s::Float64
    s1::Float64
    s2::Float64
    t::Float64
    t2::Float64
end

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

_model_vars(model::TDoubleReggeModel, m, cosθ, ϕ) =
    (s = model.reaction_system.s0, s1 = m^2, cosθ = cosθ, ϕ = ϕ, t2 = model.t2)

function _modelTR_core(
    top::TVertex,
    bot::TVertex,
    kin::TKinematics,
    K;
    α′::Float64,
)
    α1 = top.α(kin.t)
    α2 = bot.α(kin.t2)
    prefactor =
        -K * sf_gamma(1 - α1) * sf_gamma(1 - α2) * (α′ * kin.s1)^α1 *
        (α′ * kin.s2)^α2 / (α′ * kin.s)
    ff = form_factor(top, kin.t) * form_factor(bot, kin.t2)
    η = kin.s / (α′ * kin.s1 * kin.s2)
    τ1, τ2 = top.τ, bot.τ
    vertex = 0.0im
    ξ1 = ξ(τ1, α1)
    ξ21 = ξ(τ2, τ1, α2, α1)
    vertex += η^α1 * ξ1 * ξ21 * V(α1, α2, η)
    ξ2 = ξ(τ2, α2)
    ξ12 = ξ(τ1, τ2, α1, α2)
    vertex += η^α2 * ξ2 * ξ12 * V(α2, α1, η)
    return ff * prefactor * vertex
end

function modelDR(
    exchange::TReggeExchange,
    vars,
    reaction_system;
    α′::Float64 = 0.9,
)
    @unpack s, s1, t2 = vars
    t = exchange.η_forward ? t1(vars, reaction_system) : tπ(vars, reaction_system)
    s2 = exchange.η_forward ? sπp(vars, reaction_system) : sηp(vars, reaction_system)
    K = Kfactor(vars, reaction_system)
    kin = TKinematics(s, s1, s2, t, t2)
    return _modelTR_core(exchange.top, exchange.bot, kin, K; α′ = α′)
end

function amplitude(model::TDoubleReggeModel, m, cosθ, ϕ)
    vars = _model_vars(model, m, cosθ, ϕ)
    generator = (
        p * modelDR(exchange, vars, model.reaction_system; α′ = model.scalar_α)
        for (p, exchange) in zip(model.pars, model.exchanges)
    )
    return mysum(typeof(1im * model.pars[1]), generator)
end

struct TEventKinematics
    s::Float64
    s12::Float64
    s13::Float64
    s23::Float64
    t1::Float64
    tπ::Float64
    t2::Float64
    cosθ::Float64
    ϕ::Float64
    K::Float64
end

function modelDR(exchange::TReggeExchange, ev::TEventKinematics; α′::Float64 = 0.9)
    t = exchange.η_forward ? ev.t1 : ev.tπ
    s2 = exchange.η_forward ? ev.s23 : ev.s13
    kin = TKinematics(ev.s, ev.s12, s2, t, ev.t2)
    return _modelTR_core(exchange.top, exchange.bot, kin, ev.K; α′ = α′)
end

function amplitude(model::TDoubleReggeModel, ev::TEventKinematics)
    generator = (
        p * modelDR(exchange, ev; α′ = model.scalar_α)
        for (p, exchange) in zip(model.pars, model.exchanges)
    )
    return mysum(typeof(1im * model.pars[1]), generator)
end

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

function load_modelT_config(parsed; settings_key::AbstractString = "settings")
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
        model = model,
        settings = settings,
        reaction_system = reaction_system,
        trajectories = trajectories,
        vertices = vertices,
        diagrams = diagrams,
        couplings = pars,
    )
end

load_modelT_from_toml(path::AbstractString; settings_key::AbstractString = "settings") =
    load_modelT_config(TOML.parsefile(path); settings_key = settings_key)
