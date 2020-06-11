
V(x,y,η) = sf_gamma(x-y) / sf_gamma(1-y) * hypergeom(1-x,1-x+y,-1/η)

ξ(τ,α) = (τ + cis(-π*α))/2
ξ(τ1,τ2,α1,α2) = (τ1*τ2 + cis(-π*(α1-α2)))/2

const α′ = 0.9;

# trajectories
struct trajectory
    slope::Real
    intercept::Real
end
(α::trajectory)(t) = α.intercept + α.slope*t
#
const α_a2 = trajectory(0.9, 1-0.9*0.77^2)
const α_f2 = trajectory(0.89, 0.47)
const α_ℙ  = trajectory(0.25, 1.08)
#
function modelDR(α1oft::trajectory, α2oft::trajectory, vars; η_forward::Bool=true)
    @unpack s, s1, cosθ, ϕ, t2 = vars
    #
    t = η_forward ? t1(vars) : tπ(vars)
    s2 = η_forward ? sπp(vars) : sηp(vars)
    #
    α1 = α1oft(t)
    α2 = α2oft(t2)
    #
    K = Kfactor(vars)
    prefactor = -K*sf_gamma(1-α1)*sf_gamma(1-α2) *
        (α′*s1)^α1 * (α′*s2)^α2 / (α′*s)
    #
    # vertex funcions
    τ1 = τ2 = 1.0
    η = s/(α′*s1*s2)
    vertex = 0.0im
    # if η_forward
    ξ1= ξ(τ1,α1)
    ξ21 = ξ(τ2,τ1,α2,α1)
    vertex += η^α1*ξ1*ξ21*V(α1,α2,η)
    # else
    ξ2 = ξ(τ2,α2)
    ξ12 = ξ(τ1,τ2,α1,α2)
    vertex += η^α2*ξ2*ξ12*V(α2,α1,η)
    # end
    #
    return prefactor*vertex
end