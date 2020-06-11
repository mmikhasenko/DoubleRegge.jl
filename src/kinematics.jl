λ(x,y,z) = x^2+y^2+z^2-2x*y-2y*z-2z*x
function t1(vars)
    @unpack s, s1, cosθ, ϕ, t2 = vars
    return mπ2+mη2-(s1+mπ2-t2)*(s1+mη2-mπ2)/(2*s1)+sqrt(λ(s1,t2,mπ2)*λ(s1,mη2,mπ2))/(2*s1) * cosθ
end
#
function tπ(vars)
    @unpack s, s1, cosθ, ϕ, t2 = vars
    return mπ2+mπ2-(s1+mπ2-t2)*(s1+mπ2-mη2)/(2*s1)-sqrt(λ(s1,t2,mπ2)*λ(s1,mη2,mπ2))/(2*s1) * cosθ
end
#
function cosθ2(vars)
    @unpack s, s1, cosθ, ϕ, t2 = vars
    return - (2*s1*(s1-s-t2+mp2)+(s1-t2+mπ2)*(s-s1-mp2)) / sqrt( λ(s,s1,mp2) * λ(t2,s1,mπ2) )
end
#
function sπp(vars)
    @unpack s, s1, cosθ, ϕ, t2 = vars
    sinθ = sqrt(1-cosθ^2)
    _cosθ2 = cosθ2(vars)
    sinθ2 = sqrt(1-_cosθ2^2)
    #
    return mp2+mπ2+(s1+mπ2-mη2)*(s-s1-mp2) / (2s1) -
        sqrt(λ(s1,mπ2,mη2)*λ(s1,s,mp2)) / (2s1) *
            (sinθ2*sinθ*cos(ϕ)+_cosθ2*cosθ)
end

function sηp(vars)
    @unpack s, s1, cosθ, ϕ, t2 = vars
    sinθ = sqrt(1-cosθ^2)
    _cosθ2 = cosθ2(vars)
    sinθ2 = sqrt(1-_cosθ2^2)
    #
    return mp2+mπ2+(s1+mη2-mπ2) * (s-s1-mp2) / (2s1) +
        sqrt(λ(s1,mπ2,mη2)*λ(s1,s,mp2)) / (2s1) *
            (sinθ2*sinθ*cos(ϕ)+_cosθ2*cosθ)
end

function Kfactor(vars)
    @unpack s, s1, cosθ, ϕ, t2 = vars
    #
    sinθ = sqrt(1-cosθ^2)
    #
    _cosθ2 = cosθ2(vars)
    sinθ2 = sqrt(1-_cosθ2^2)
    #
    p2 = sqrt(λ(s1,s,mp2)/(4*s1));
    q  = sqrt(λ(s1,t2,mπ2)/(4*s1));
    k  = sqrt(λ(s1,mη2,mπ2)/(4*s1));
    #
    K = 4*sqrt(s1) * p2 * q * k * sinθ2 * sinθ * sin(ϕ)
    return K
end
