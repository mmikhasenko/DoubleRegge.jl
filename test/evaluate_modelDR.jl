using DrWatson
@quickactivate "DoubleRegge"

using Test
using DoubleRegge

const ηπ_system = compass_ηπ

let
    @test α_a2 isa trajectory
    @test α_f2 isa trajectory
    @test α_ℙ isa trajectory
end


let
    vars = (s = ηπ_system.s0, s1 = 4.0^2, cosθ = 0.4, ϕ = π/4, t2 = -0.45)
    v = modelDR(α_a2, α_ℙ, vars, ηπ_system; η_forward = true)
    @test v != 0.0
end

let
    exchange = ReggeExchange(α_a2, α_ℙ, true, "a2/ℙ")
    vars = (s = ηπ_system.s0, s1 = 4.0^2, cosθ = 0.4, ϕ = π/4, t2 = -0.45)
    direct = modelDR(α_a2, α_ℙ, vars, ηπ_system; η_forward = true, α′ = 0.85, s2shift = 0.1)
    typed = modelDR(exchange, vars, ηπ_system; α′ = 0.85, s2shift = 0.1)
    @test typed ≈ direct
    @test exchange[1] === α_a2
    @test exchange[2] === α_ℙ
    @test exchange[3] == true
    @test exchange[4] == "a2/ℙ"
end

let
    vars = (s = ηπ_system.s0, s1 = 4.0^2, cosθ = 0.4, ϕ = π/4, t2 = -0.45)
    @test DoubleRegge.sπp(vars, ηπ_system) != DoubleRegge.sηp(vars, ηπ_system)
end

#
get_control_variables(vars, reaction_system) = (
    sπp = DoubleRegge.sπp(vars, reaction_system),
    sηp = DoubleRegge.sηp(vars, reaction_system),
    tη = DoubleRegge.t1(vars, reaction_system),
    tπ = DoubleRegge.tπ(vars, reaction_system),
    K = DoubleRegge.Kfactor(vars, reaction_system),
    A = modelDR(α_a2, α_ℙ, vars, reaction_system; η_forward = true),
)

let
    vars = (s = ηπ_system.s0, s1 = 2.32^2, cosθ = cos(π/3), ϕ = π/4, t2 = -0.2)
    contv = get_control_variables(vars, ηπ_system)
    vincv = (52.88, 300.37, -1.32, -3.92, 211.22, 8.859 + 13.81im)
    # [println("$k: $c $v") for (k,c,v) in zip(keys(contv),contv,vincv)]
    @test sum(abs2(c-v) for (c, v) in zip(contv, vincv)) < 0.1
end

let
    vars = (s = ηπ_system.s0, s1 = 2.32^2, cosθ = cos(2π/3), ϕ = -π/4, t2 = -0.2)
    contv = get_control_variables(vars, ηπ_system)
    vincv = (205.89, 147.37, -3.936, -1.3, -211.22, 0.4789 - 2.848im)
    # [println("$k: $c $v") for (k,c,v) in zip(keys(contv),contv,vincv)]
    @test sum(abs2(c-v) for (c, v) in zip(contv, vincv)) < 0.1
end  #

@testset "DoubleReggeModel parameter binding" begin
    pars = [1.0, -0.5]
    model = DoubleReggeModel(six_exchanges[1:2], -0.2, 0.8, ηπ_system, pars; s2shift = 0.05)
    updated = with_parameters(model, [2.0, 3.0])
    @test updated.pars == [2.0, 3.0]
    @test updated.exchanges === model.exchanges
    @test updated.t2 == model.t2
    @test updated.scalar_α == model.scalar_α
    @test updated.reaction_system === model.reaction_system
    @test updated.s2shift == model.s2shift
end
