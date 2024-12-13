using DrWatson
@quickactivate "DoubleRegge"

using Test
using DoubleRegge

setsystem!(:compass_ηπ)

let
    @test typeof(α_a2) == trajectory
    @test typeof(α_f2) == trajectory
    @test typeof(α_ℙ) == trajectory
end


let
    vars = (s = G.s0, s1 = 4.0^2, cosθ = 0.4, ϕ = π/4, t2 = -0.45)
    v = modelDR(α_a2, α_ℙ, vars; η_forward=true)
    @test v != 0.0
end

let
    vars = (s = DoubleRegge.G.s0, s1 = 4.0^2, cosθ = 0.4, ϕ = π/4, t2 = -0.45)
    @test DoubleRegge.sπp(vars) != DoubleRegge.sηp(vars)
end

#
get_control_variables(vars) =
    (sπp = DoubleRegge.sπp(vars),
     sηp = DoubleRegge.sηp(vars),
     tη = DoubleRegge.t1(vars),
     tπ = DoubleRegge.tπ(vars),
     K = DoubleRegge.Kfactor(vars),
     A = modelDR(α_a2, α_ℙ, vars; η_forward=true))

let
    vars = (s = DoubleRegge.G.s0, s1 = 2.32^2, cosθ = cos(π/3), ϕ = π/4, t2 = -0.2)
    contv = get_control_variables(vars)
    vincv = (52.88, 300.37, -1.32, -3.92,  211.22, 8.859 + 13.81im)
    # [println("$k: $c $v") for (k,c,v) in zip(keys(contv),contv,vincv)]
    @test sum(abs2(c-v) for (c,v) in zip(contv,vincv)) < 0.1
end

let
    vars = (s = DoubleRegge.G.s0, s1 = 2.32^2, cosθ = cos(2π/3), ϕ = -π/4, t2 = -0.2)
    contv = get_control_variables(vars)
    vincv = (205.89, 147.37, -3.936, -1.3, -211.22,  0.4789 - 2.848im)
    # [println("$k: $c $v") for (k,c,v) in zip(keys(contv),contv,vincv)]
    @test sum(abs2(c-v) for (c,v) in zip(contv,vincv)) < 0.1
end  #
