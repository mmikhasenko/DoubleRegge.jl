function modelS(x, vars)
    @unpack s, s1, cosθ, ϕ, t2 = vars
    _t1 = t1(vars); _tπ = tπ(vars)
    K = Kfactor(vars)
    #
    α(t) = 0.5 + 0.9*t
    β = 9.0;
    exp(β*t2/2) * K * ( exp(α(_t1)) - x[1]*exp(α(_tπ)) )
end
