
V(x, y, η) = gamma(x-y) / gamma(1-y) * hypergeom(1-x, 1-x+y, -1/η) #pFq([1-x],[1-x+y],-1/η) #

ξ(τ, α) = (τ + cis(-π*α))/2
ξ(τ1, τ2, α1, α2) = (τ1*τ2 + cis(-π*(α1-α2)))/2

# trajectories
struct trajectory{T}
    slope::T
    intercept::T
end
(α::trajectory)(t) = α.intercept + α.slope*t
#
const α_a2 = trajectory(0.917,  0.44)
const α_f2 = trajectory(0.917,  0.44)
const α_ℙ  = trajectory(0.25,   1.08)
const α_a2′ = trajectory(0.917, -0.56)
const α_π1  = trajectory(0.68,  -0.74)

# vertex: pairs a Regge trajectory with a form-factor slope b and signature τ
struct Vertex{T}
    α::T        # Regge trajectory
    b::Float64  # form factor slope: exp(b * t) at this vertex
    τ::Float64  # signature factor: +1 (normal) or -1 (π1 top vertex)
end

struct ReggeExchange{T1,T2,S<:AbstractString}
    top::Vertex{T1}
    bot::Vertex{T2}
    η_forward::Bool
    label::S
end

struct DoubleReggeModel{E<:AbstractVector,P,R}
    exchanges::E
    t2::Float64
    scalar_α::Float64
    reaction_system::R
    pars::P
    s2shift::Float64
end

function DoubleReggeModel(
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
    return DoubleReggeModel(
        exchanges,
        Float64(t2),
        Float64(scalar_α),
        reaction_system,
        pars,
        Float64(s2shift),
    )
end

with_parameters(model::DoubleReggeModel, pars) = DoubleReggeModel(
    model.exchanges,
    model.t2,
    model.scalar_α,
    model.reaction_system,
    pars;
    s2shift = model.s2shift,
)

_model_vars(model::DoubleReggeModel, m, cosθ, ϕ) =
    (s = model.reaction_system.s0, s1 = m^2, cosθ = cosθ, ϕ = ϕ, t2 = model.t2)

# --- existing path: fixed kinematics computed from (m, cosθ, ϕ) ---

function modelDR(
    α1oft::trajectory,
    α2oft::trajectory,
    vars,
    system;
    η_forward::Bool = true,
    α′::Float64 = 0.9,
    s2shift::Float64 = 0.0,
    τ1::Float64 = 1.0,
    τ2::Float64 = 1.0,
    b_top::Float64 = 0.0,
    b_bot::Float64 = 0.0,
)
    #
    @unpack s, s1, cosθ, ϕ, t2 = vars
    #
    t = η_forward ? t1(vars, system) : tπ(vars, system)
    s2 = η_forward ? sπp(vars, system) : sηp(vars, system)
    #
    α1 = α1oft(t)
    α2 = α2oft(t2)
    #
    K = Kfactor(vars, system)
    prefactor =
        -K * sf_gamma(1-α1) * sf_gamma(1-α2) * (α′*s1)^α1 * (α′*(s2+s2shift))^α2 / (α′*s)
    #
    form_factor = exp(b_top * t) * exp(b_bot * t2)
    #
    # vertex functions
    η = s/(α′*s1*s2)
    vertex = 0.0im
    ξ1 = ξ(τ1, α1)
    ξ21 = ξ(τ2, τ1, α2, α1)
    vertex += η^α1*ξ1*ξ21*V(α1, α2, η)
    ξ2 = ξ(τ2, α2)
    ξ12 = ξ(τ1, τ2, α1, α2)
    vertex += η^α2*ξ2*ξ12*V(α2, α1, η)
    return form_factor * prefactor * vertex
end

modelDR(
    exchange::ReggeExchange,
    vars,
    reaction_system;
    α′::Float64 = 0.9,
    s2shift::Float64 = 0.0,
) = modelDR(
    exchange.top.α,
    exchange.bot.α,
    vars,
    reaction_system;
    η_forward = exchange.η_forward,
    α′ = α′,
    s2shift = s2shift,
    τ1 = exchange.top.τ,
    τ2 = exchange.bot.τ,
    b_top = exchange.top.b,
    b_bot = exchange.bot.b,
)

function amplitude(model::DoubleReggeModel, m, cosθ, ϕ)
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

# --- new path: pre-calculated event kinematics ---

struct EventKinematics
    s::Float64    # total s (var.s)
    s12::Float64  # m²(ηπ)  (var.s12)
    s13::Float64  # m²(ηp)  (var.s13)
    s23::Float64  # m²(πp)  (var.s23)
    t1::Float64   # momentum transfer to η (var.t1_1)
    tπ::Float64   # momentum transfer to π (var.t1_2)
    t2::Float64   # momentum transfer to proton (var.t2)
    cosθ::Float64 # CosThetaGJ
    ϕ::Float64    # PhiGJ
    K::Float64    # KFactor
end

function modelDR(
    exchange::ReggeExchange,
    ev::EventKinematics;
    α′::Float64 = 0.9,
    s2shift::Float64 = 0.0,
)
    t  = exchange.η_forward ? ev.t1  : ev.tπ
    s2 = exchange.η_forward ? ev.s23 : ev.s13
    α1 = exchange.top.α(t)
    α2 = exchange.bot.α(ev.t2)
    prefactor =
        -ev.K * sf_gamma(1-α1) * sf_gamma(1-α2) *
        (α′*ev.s12)^α1 * (α′*(s2+s2shift))^α2 / (α′*ev.s)
    form_factor = exp(exchange.top.b * t) * exp(exchange.bot.b * ev.t2)
    τ1 = exchange.top.τ
    τ2 = exchange.bot.τ
    η_var = ev.s / (α′ * ev.s12 * s2)
    vertex = 0.0im
    ξ1  = ξ(τ1, α1);  ξ21 = ξ(τ2, τ1, α2, α1)
    vertex += η_var^α1 * ξ1 * ξ21 * V(α1, α2, η_var)
    ξ2  = ξ(τ2, α2);  ξ12 = ξ(τ1, τ2, α1, α2)
    vertex += η_var^α2 * ξ2 * ξ12 * V(α2, α1, η_var)
    return form_factor * prefactor * vertex
end

function amplitude(model::DoubleReggeModel, ev::EventKinematics)
    generator = (
        p * modelDR(exchange, ev; α′ = model.scalar_α, s2shift = model.s2shift)
        for (p, exchange) in zip(model.pars, model.exchanges)
    )
    return mysum(typeof(1im * model.pars[1]), generator)
end

# --- named vertices for ten_exchanges ---

const v_ππ1η  = Vertex(α_π1,  -0.0882,  -1.0)  # top, η_forward=true,  τ=-1
const v_πa2η  = Vertex(α_a2,  -0.192,    1.0)  # top, η_forward=true
const v_πa2′η = Vertex(α_a2′,  38.0,     1.0)  # top, η_forward=true
const v_πf2π  = Vertex(α_f2,   1.432,    1.0)  # top, η_forward=false
const v_πℙπ   = Vertex(α_ℙ,    0.681,    1.0)  # top, η_forward=false
const v_pf2p  = Vertex(α_f2,   1.054,    1.0)  # bottom
const v_pℙp   = Vertex(α_ℙ,    1.468,    1.0)  # bottom

const ten_exchanges = ReggeExchange[
    ReggeExchange(v_ππ1η,  v_pf2p, true,  "π1/f2"),
    ReggeExchange(v_ππ1η,  v_pℙp,  true,  "π1/ℙ"),
    ReggeExchange(v_πa2η,  v_pf2p, true,  "a2/f2"),
    ReggeExchange(v_πa2η,  v_pℙp,  true,  "a2/ℙ"),
    ReggeExchange(v_πa2′η, v_pf2p, true,  "a2′/f2"),
    ReggeExchange(v_πa2′η, v_pℙp,  true,  "a2′/ℙ"),
    ReggeExchange(v_πf2π,  v_pf2p, false, "f2/f2"),
    ReggeExchange(v_πf2π,  v_pℙp,  false, "f2/ℙ"),
    ReggeExchange(v_πℙπ,   v_pf2p, false, "ℙ/f2"),
    ReggeExchange(v_πℙπ,   v_pℙp,  false, "ℙ/ℙ"),
]

# --- six_exchanges rebuilt with Vertex wrappers (b=0, τ=+1): identical physics ---

const six_exchanges = ReggeExchange[
    ReggeExchange(Vertex(α_a2, 0.0, 1.0), Vertex(α_ℙ,  0.0, 1.0), true,  "a2/ℙ"),
    ReggeExchange(Vertex(α_a2, 0.0, 1.0), Vertex(α_f2, 0.0, 1.0), true,  "a2/f2"),
    ReggeExchange(Vertex(α_f2, 0.0, 1.0), Vertex(α_ℙ,  0.0, 1.0), false, "f2/ℙ"),
    ReggeExchange(Vertex(α_f2, 0.0, 1.0), Vertex(α_f2, 0.0, 1.0), false, "f2/f2"),
    ReggeExchange(Vertex(α_ℙ,  0.0, 1.0), Vertex(α_ℙ,  0.0, 1.0), false, "ℙ/ℙ"),
    ReggeExchange(Vertex(α_ℙ,  0.0, 1.0), Vertex(α_f2, 0.0, 1.0), false, "ℙ/f2"),
]
