const αT_a2 = trajectory(0.917, 0.44)
const αT_f2 = trajectory(0.917, 0.44)
const αT_ℙ = trajectory(0.25, 1.08)
const αT_a2′ = trajectory(0.917, -0.56)
const αT_π1 = trajectory(0.68, -0.74)

struct TVertex{T}
    α::T
    b::Float64
    τ::Float64
end

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

Base.getindex(exchange::TReggeExchange, i::Int) =
    i == 1 ? exchange.top.α :
    i == 2 ? exchange.bot.α :
    i == 3 ? exchange.η_forward :
    i == 4 ? exchange.label :
    throw(BoundsError(exchange, i))

Base.firstindex(::TReggeExchange) = 1
Base.lastindex(::TReggeExchange) = 4

struct TDoubleReggeModel{E<:AbstractVector,P,R}
    exchanges::E
    t2::Float64
    scalar_α::Float64
    reaction_system::R
    pars::P
    s2shift::Float64
end

function TDoubleReggeModel(
    exchanges,
    t2::Real,
    scalar_α::Real,
    reaction_system,
    pars;
    s2shift::Real = 0.0,
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
        Float64(s2shift),
    )
end

with_parameters(model::TDoubleReggeModel, pars) = TDoubleReggeModel(
    model.exchanges,
    model.t2,
    model.scalar_α,
    model.reaction_system,
    pars;
    s2shift = model.s2shift,
)

_model_vars(model::TDoubleReggeModel, m, cosθ, ϕ) =
    (s = model.reaction_system.s0, s1 = m^2, cosθ = cosθ, ϕ = ϕ, t2 = model.t2)

function _modelTR_core(
    α1,
    α2,
    s,
    s1,
    s2,
    t,
    t2,
    K;
    α′::Float64 = 0.9,
    s2shift::Float64 = 0.0,
    τ1::Float64 = 1.0,
    τ2::Float64 = 1.0,
    b_top::Float64 = 0.0,
    b_bot::Float64 = 0.0,
)
    prefactor =
        -K * sf_gamma(1 - α1) * sf_gamma(1 - α2) * (α′ * s1)^α1 *
        (α′ * (s2 + s2shift))^α2 / (α′ * s)
    form_factor = exp(b_top * t) * exp(b_bot * t2)
    η = s / (α′ * s1 * s2)
    vertex = 0.0im
    ξ1 = ξ(τ1, α1)
    ξ21 = ξ(τ2, τ1, α2, α1)
    vertex += η^α1 * ξ1 * ξ21 * V(α1, α2, η)
    ξ2 = ξ(τ2, α2)
    ξ12 = ξ(τ1, τ2, α1, α2)
    vertex += η^α2 * ξ2 * ξ12 * V(α2, α1, η)
    return form_factor * prefactor * vertex
end

function modelDR(
    α1oft::trajectory,
    α2oft::trajectory,
    vars,
    system,
    ::Type{TDoubleReggeModel};
    η_forward::Bool = true,
    α′::Float64 = 0.9,
    s2shift::Float64 = 0.0,
    τ1::Float64 = 1.0,
    τ2::Float64 = 1.0,
    b_top::Float64 = 0.0,
    b_bot::Float64 = 0.0,
)
    @unpack s, s1, cosθ, ϕ, t2 = vars
    t = η_forward ? t1(vars, system) : tπ(vars, system)
    s2 = η_forward ? sπp(vars, system) : sηp(vars, system)
    α1 = α1oft(t)
    α2 = α2oft(t2)
    K = Kfactor(vars, system)
    return _modelTR_core(
        α1,
        α2,
        s,
        s1,
        s2,
        t,
        t2,
        K;
        α′ = α′,
        s2shift = s2shift,
        τ1 = τ1,
        τ2 = τ2,
        b_top = b_top,
        b_bot = b_bot,
    )
end

modelDR(
    exchange::TReggeExchange,
    vars,
    reaction_system;
    α′::Float64 = 0.9,
    s2shift::Float64 = 0.0,
) = modelDR(
    exchange.top.α,
    exchange.bot.α,
    vars,
    reaction_system,
    TDoubleReggeModel;
    η_forward = exchange.η_forward,
    α′ = α′,
    s2shift = s2shift,
    τ1 = exchange.top.τ,
    τ2 = exchange.bot.τ,
    b_top = exchange.top.b,
    b_bot = exchange.bot.b,
)

function amplitude(model::TDoubleReggeModel, m, cosθ, ϕ)
    vars = _model_vars(model, m, cosθ, ϕ)
    generator = (
        p * modelDR(
            exchange,
            vars,
            model.reaction_system;
            α′ = model.scalar_α,
            s2shift = model.s2shift,
        ) for (p, exchange) in zip(model.pars, model.exchanges)
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

function modelDR(
    exchange::TReggeExchange,
    ev::TEventKinematics;
    α′::Float64 = 0.9,
    s2shift::Float64 = 0.0,
)
    t = exchange.η_forward ? ev.t1 : ev.tπ
    s2 = exchange.η_forward ? ev.s23 : ev.s13
    α1 = exchange.top.α(t)
    α2 = exchange.bot.α(ev.t2)
    return _modelTR_core(
        α1,
        α2,
        ev.s,
        ev.s12,
        s2,
        t,
        ev.t2,
        ev.K;
        α′ = α′,
        s2shift = s2shift,
        τ1 = exchange.top.τ,
        τ2 = exchange.bot.τ,
        b_top = exchange.top.b,
        b_bot = exchange.bot.b,
    )
end

function amplitude(model::TDoubleReggeModel, ev::TEventKinematics)
    generator = (
        p * modelDR(exchange, ev; α′ = model.scalar_α, s2shift = model.s2shift)
        for (p, exchange) in zip(model.pars, model.exchanges)
    )
    return mysum(typeof(1im * model.pars[1]), generator)
end

const vT_ππ1η = TVertex(αT_π1, -0.0882, -1.0)
const vT_πa2η = TVertex(αT_a2, -0.192, 1.0)
const vT_πa2′η = TVertex(αT_a2′, 38.0, 1.0)
const vT_πf2π = TVertex(αT_f2, 1.432, 1.0)
const vT_πℙπ = TVertex(αT_ℙ, 0.681, 1.0)
const vT_pf2p = TVertex(αT_f2, 1.054, 1.0)
const vT_pℙp = TVertex(αT_ℙ, 1.468, 1.0)

const ten_exchanges_T = TReggeExchange[
    TReggeExchange(vT_ππ1η, vT_pf2p, true, "π1/f2"),
    TReggeExchange(vT_ππ1η, vT_pℙp, true, "π1/ℙ"),
    TReggeExchange(vT_πa2η, vT_pf2p, true, "a2/f2"),
    TReggeExchange(vT_πa2η, vT_pℙp, true, "a2/ℙ"),
    TReggeExchange(vT_πa2′η, vT_pf2p, true, "a2′/f2"),
    TReggeExchange(vT_πa2′η, vT_pℙp, true, "a2′/ℙ"),
    TReggeExchange(vT_πf2π, vT_pf2p, false, "f2/f2"),
    TReggeExchange(vT_πf2π, vT_pℙp, false, "f2/ℙ"),
    TReggeExchange(vT_πℙπ, vT_pf2p, false, "ℙ/f2"),
    TReggeExchange(vT_πℙπ, vT_pℙp, false, "ℙ/ℙ"),
]
