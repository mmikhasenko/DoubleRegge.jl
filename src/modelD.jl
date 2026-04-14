
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
const α_a2 = trajectory(0.9, 1-0.9*0.77^2)
const α_f2 = trajectory(0.89, 0.47)
const α_ℙ = trajectory(0.25, 1.08)
#
struct ReggeExchange{T1,T2,S<:AbstractString}
    α_top::T1
    α_bot::T2
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

function modelDR(
    α1oft::trajectory,
    α2oft::trajectory,
    vars,
    system;
    η_forward::Bool = true,
    α′::Float64 = 0.9,
    s2shift::Float64 = 0.0,
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
    # vertex funcions
    τ1 = τ2 = 1.0
    η = s/(α′*s1*s2)
    vertex = 0.0im
    # if η_forward
    ξ1 = ξ(τ1, α1)
    ξ21 = ξ(τ2, τ1, α2, α1)
    vertex += η^α1*ξ1*ξ21*V(α1, α2, η)
    # else
    ξ2 = ξ(τ2, α2)
    ξ12 = ξ(τ1, τ2, α1, α2)
    vertex += η^α2*ξ2*ξ12*V(α2, α1, η)
    # end
    return prefactor*vertex
end

modelDR(
    exchange::ReggeExchange,
    vars,
    reaction_system;
    α′::Float64 = 0.9,
    s2shift::Float64 = 0.0,
) = modelDR(
    exchange.α_top,
    exchange.α_bot,
    vars,
    reaction_system;
    η_forward = exchange.η_forward,
    α′ = α′,
    s2shift = s2shift,
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

const six_exchanges = ReggeExchange[
    ReggeExchange(α_a2, α_ℙ, true, "a2/ℙ"),
    ReggeExchange(α_a2, α_f2, true, "a2/f2"),
    ReggeExchange(α_f2, α_ℙ, false, "f2/ℙ"),
    ReggeExchange(α_f2, α_f2, false, "f2/f2"),
    ReggeExchange(α_ℙ, α_ℙ, false, "ℙ/ℙ"),
    ReggeExchange(α_ℙ, α_f2, false, "ℙ/f2"),
];
