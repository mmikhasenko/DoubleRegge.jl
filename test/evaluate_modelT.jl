using DrWatson
@quickactivate "DoubleRegge"

using Test
using DoubleRegge

const ηπ_system_T = compass_ηπ

@testset "modelT parameter binding" begin
    pars = [1.0, -0.5]
    model = TDoubleReggeModel(ten_exchanges_T[1:2], -0.2, 0.8, ηπ_system_T, pars; s2shift = 0.05)
    updated = with_parameters(model, [2.0, 3.0])
    @test updated.pars == [2.0, 3.0]
    @test updated.exchanges === model.exchanges
    @test updated.t2 == model.t2
    @test updated.scalar_α == model.scalar_α
    @test updated.reaction_system === model.reaction_system
    @test updated.s2shift == model.s2shift
end

@testset "modelT exchange compatibility" begin
    exchange = TReggeExchange(αT_a2, αT_ℙ, true, "a2/ℙ")
    vars = (s = ηπ_system_T.s0, s1 = 4.0^2, cosθ = 0.4, ϕ = π / 4, t2 = -0.45)
    direct = modelDR(
        αT_a2,
        αT_ℙ,
        vars,
        ηπ_system_T,
        TDoubleReggeModel;
        η_forward = true,
        α′ = 0.85,
        s2shift = 0.1,
    )
    typed = modelDR(exchange, vars, ηπ_system_T; α′ = 0.85, s2shift = 0.1)
    @test typed ≈ direct
    @test exchange[1] === αT_a2
    @test exchange[2] === αT_ℙ
    @test exchange[3] == true
    @test exchange[4] == "a2/ℙ"
end

@testset "modelT EventKinematics path" begin
    vars = (s = ηπ_system_T.s0, s1 = 2.32^2, cosθ = cos(π / 3), ϕ = π / 4, t2 = -0.2)
    ev = TEventKinematics(
        vars.s,
        vars.s1,
        DoubleRegge.sηp(vars, ηπ_system_T),
        DoubleRegge.sπp(vars, ηπ_system_T),
        DoubleRegge.t1(vars, ηπ_system_T),
        DoubleRegge.tπ(vars, ηπ_system_T),
        vars.t2,
        vars.cosθ,
        vars.ϕ,
        DoubleRegge.Kfactor(vars, ηπ_system_T),
    )
    exchange = ten_exchanges_T[3]
    @test modelDR(exchange, vars, ηπ_system_T; α′ = 0.8, s2shift = 0.1) ≈
          modelDR(exchange, ev; α′ = 0.8, s2shift = 0.1)
end

@testset "modelT amplitude runs" begin
    model = TDoubleReggeModel(ten_exchanges_T, -0.2, 0.8, ηπ_system_T, ones(length(ten_exchanges_T)))
    amp = amplitude(model, 2.3, 0.4, π / 5)
    @test isfinite(real(amp))
    @test isfinite(imag(amp))
end
