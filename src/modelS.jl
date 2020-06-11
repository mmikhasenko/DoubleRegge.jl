const β = 9.0;

λ(x,y,z) = x^2+y^2+z^2-2x*y-2y*z-2z*x
function t1(vars)
    @unpack s, s1, cosθ, ϕ, t2 = vars
    return mπ2+mη2-(s1-t2+mπ2)*(s1+mη2-mπ2)/(2*s1)+sqrt(λ(s1,t2,mπ2)*λ(s1,mη2,mπ2))/(2*s1) * cosθ
end
function tπ(vars)
    @unpack s, s1, cosθ, ϕ, t2 = vars
    return 2mπ2-(s1-t2+mπ2)*(s1-mη2+mπ2)/(2s1)-sqrt(λ(s1,t2,mπ2)*λ(s1,mη2,mπ2))/(2*s1) * cosθ
end
function cosθ2(vars)
    @unpack s, s1, cosθ, ϕ, t2 = vars
    return - (2*s1*(s1-s-t2+mp2)+(s1-t2+mπ2)*(s-s1-mp2)) / sqrt( λ(s,s1,mp2) * λ(t2,s1,mπ2) )
end
#
α(t) = 0.5 + 0.9*t
#
function modelS(x, vars)
    @unpack s, s1, cosθ, ϕ, t2 = vars
    sinθ = sqrt(1-cosθ^2)
    #
    _cosθ2 = cosθ2(vars)
    sinθ2 = sqrt(1-_cosθ2^2)
    #
    p2 = sqrt(λ(s,s1,mp2)/(4*s));
    q  = sqrt(λ(s,s1,mp2)/(4*s));
    k  = sqrt(λ(s,s1,mp2)/(4*s));
    #
    K = 4*sqrt(s1) * p2 * q * k * sinθ2 * sinθ * sin(ϕ)
    _t1 = t1(vars); _tπ = tπ(vars)
    exp(β*t2/2) * K * ( exp(α(_t1)) - x[1]*exp(α(_tπ)) )
end
# plot
# plot(cosθ->A((s = s0, s1 = 4.0^2, cosθ = cosθ, ϕ = π/4, t2 = -0.45)), -1, 1)
