λ(x,y,z) = x^2+y^2+z^2-2x*y-2y*z-2z*x
function t1(vars, system=G.system)
    @unpack s, s1, cosθ, ϕ, t2 = vars
    @unpack mb², mt², mr², m1², m2² = system^2
    return mb²+m1²-(s1+mb²-t2)*(s1+m1²-m2²)/(2*s1)+sqrt(λ(s1,t2,mb²)*λ(s1,m1²,m2²))/(2*s1) * cosθ
end
#
function tπ(vars, system=G.system)
    @unpack s, s1, cosθ, ϕ, t2 = vars
    @unpack mb², mt², mr², m1², m2² = system^2
    return mb²+m2²-(s1+mb²-t2)*(s1+m2²-m1²)/(2*s1)-sqrt(λ(s1,t2,mb²)*λ(s1,m1²,m2²))/(2*s1) * cosθ
end
#
function cosθ2(vars, system=G.system) # to be checked mp^2 -> mr²
    @unpack s, s1, cosθ, ϕ, t2 = vars
    @unpack mb², mt², mr², m1², m2² = system^2
    return -(2s1*(s1-s-t2+mr²)+(s1-t2+mb²)*(s-s1-mr²)) / sqrt(λ(s,s1,mr²)*λ(t2,s1,mb²))
end
#
function cosθ1(vars, system=G.system)
    @unpack s1, t1, t2 = vars
    @unpack mb², mt², mr², m1², m2² = system^2
    return (2s1*(t1-mb²-m1²)+(s1-t2+mb²)*(s1+m1²-m2²)) / sqrt(λ(s1,t2,mb²)*λ(s1,m1²,m2²))
end

#
function sπp(vars, system=G.system)
    @unpack s, s1, cosθ, ϕ, t2 = vars
    @unpack mb², mt², mr², m1², m2² = system^2
    sinθ = sqrt(1-cosθ^2)
    _cosθ2 = cosθ2(vars)
    abs(_cosθ2) > 1 && error("|cos|>1: cos=$(_cosθ2), vars: ", vars)
    sinθ2 = sqrt(1-_cosθ2^2)
    #
    return mr²+m2²+(s1+m2²-m1²)*(s-s1-mr²) / (2s1) -
        sqrt(λ(s1,m2²,m1²)*λ(s1,s,mr²)) / (2s1) *
            (sinθ2*sinθ*cos(ϕ)+_cosθ2*cosθ)
end

function sηp(vars, system=G.system)
    @unpack s, s1, cosθ, ϕ, t2 = vars
    @unpack mb², mt², mr², m1², m2² = system^2
    sinθ = sqrt(1-cosθ^2)
    _cosθ2 = cosθ2(vars)
    sinθ2 = sqrt(1-_cosθ2^2)
    #
    return mr²+m2²+(s1+m1²-m2²) * (s-s1-mr²) / (2s1) +
        sqrt(λ(s1,m2²,m1²)*λ(s1,s,mr²)) / (2s1) *
            (sinθ2*sinθ*cos(ϕ)+_cosθ2*cosθ)
end

function Kfactor(vars, system=G.system)
    @unpack s, s1, cosθ, ϕ, t2 = vars
    @unpack mb², mt², mr², m1², m2² = system^2
    #
    sinθ = sqrt(1-cosθ^2)
    #
    _cosθ2 = cosθ2(vars)
    sinθ2 = sqrt(1-_cosθ2^2)
    #
    p2 = sqrt(λ(s1,s,mr²)/(4*s1));
    q  = sqrt(λ(s1,t2,mb²)/(4*s1));
    k  = sqrt(λ(s1,m1²,m2²)/(4*s1));
    #
    K = 4*sqrt(s1) * p2 * q * k * sinθ2 * sinθ * sin(ϕ)
    return K
end

q(m, system=G.system) = sqrt(λ(m^2, system.m1^2, system.m2^2))/m
