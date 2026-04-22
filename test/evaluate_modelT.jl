using DrWatson
@quickactivate "DoubleRegge"

using Test
using DoubleRegge
using TOML

const ηπ_system_T = compass_ηπ

@testset "modelT parameter binding" begin
    parsed = TOML.parsefile(joinpath(@__DIR__, "..", "data", "exp_pro", "my_model", "model-settings.toml"))
    config = load_modelT_config(parsed)
    pars = [1.0, -0.5]
    model = TDoubleReggeModel(config.model.exchanges[1:2], -0.2, 0.8, ηπ_system_T, pars)
    updated = with_parameters(model, [2.0, 3.0])
    @test updated.pars == [2.0, 3.0]
    @test updated.exchanges === model.exchanges
    @test updated.t2 == model.t2
    @test updated.scalar_α == model.scalar_α
    @test updated.reaction_system === model.reaction_system
end

@testset "modelT form factor dispatch" begin
    α_top = trajectory(0.917, 0.44)
    α_bot = trajectory(0.25, 1.08)
    top = TVertex(α_top, -0.5, 1.0)
    bot = TVertex(α_bot, 0.3, -1.0)
    @test form_factor(top, -0.4) ≈ exp(-0.5 * -0.4)
    @test form_factor(bot, -0.1) ≈ exp(0.3 * -0.1)
end

@testset "modelT exchange runs" begin
    α_top = trajectory(0.917, 0.44)
    α_bot = trajectory(0.25, 1.08)
    exchange = TReggeExchange(α_top, α_bot, true, "a2/ℙ")
    @test exchange.top.α === α_top
    @test exchange.bot.α === α_bot
    @test exchange.η_forward == true
    @test exchange.label == "a2/ℙ"
    vars = (s = ηπ_system_T.s0, s1 = 4.0^2, cosθ = 0.4, ϕ = π / 4, t2 = -0.45)
    amp = modelDR(exchange, vars, ηπ_system_T; α′ = 0.85)
    @test isfinite(real(amp))
    @test isfinite(imag(amp))
end

@testset "modelT TKinematics + _modelTR_core" begin
    α_top = trajectory(0.917, 0.44)
    α_bot = trajectory(0.25, 1.08)
    vars = (s = ηπ_system_T.s0, s1 = 4.0^2, cosθ = 0.4, ϕ = π / 4, t2 = -0.45)
    t = DoubleRegge.t1(vars, ηπ_system_T)
    s2 = DoubleRegge.sπp(vars, ηπ_system_T)
    K = DoubleRegge.Kfactor(vars, ηπ_system_T)
    kin = TKinematics(vars.s, vars.s1, s2, t, vars.t2)

    exchange = TReggeExchange(α_top, α_bot, true, "a2/ℙ")
    direct = DoubleRegge._modelTR_core(exchange.top, exchange.bot, kin, K; α′ = 0.85)
    via_vars = modelDR(exchange, vars, ηπ_system_T; α′ = 0.85)
    @test direct ≈ via_vars

    b_top, b_bot = -0.5, 0.3
    ex_ff = TReggeExchange(
        TVertex(α_top, b_top, 1.0),
        TVertex(α_bot, b_bot, 1.0),
        true,
        "a2/ℙ-ff",
    )
    amp_ff = modelDR(ex_ff, vars, ηπ_system_T; α′ = 0.85)
    @test amp_ff ≈ via_vars * form_factor(ex_ff.top, t) * form_factor(ex_ff.bot, vars.t2)
end

@testset "modelT EventKinematics path" begin
    parsed = TOML.parsefile(joinpath(@__DIR__, "..", "data", "exp_pro", "my_model", "model-settings.toml"))
    config = load_modelT_config(parsed)
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
    exchange = config.model.exchanges[3]
    @test modelDR(exchange, vars, ηπ_system_T; α′ = 0.8) ≈
          modelDR(exchange, ev; α′ = 0.8)
end

@testset "modelT amplitude runs" begin
    parsed = TOML.parsefile(joinpath(@__DIR__, "..", "data", "exp_pro", "my_model", "model-settings.toml"))
    config = load_modelT_config(parsed)
    model = TDoubleReggeModel(
        config.model.exchanges,
        -0.2,
        0.8,
        ηπ_system_T,
        ones(length(config.model.exchanges)),
    )
    amp = amplitude(model, 2.3, 0.4, π / 5)
    @test isfinite(real(amp))
    @test isfinite(imag(amp))
end

@testset "modelT config loader" begin
    parsed = TOML.parsefile(joinpath(@__DIR__, "..", "data", "exp_pro", "my_model", "model-settings.toml"))
    config = load_modelT_config(parsed)
    @test config.reaction_system === compass_ηπ
    @test length(config.model.exchanges) == 10
    @test length(config.model.pars) == 10
    @test config.model.pars[1] == parsed["couplings"]["c_pi1_f2"]
    @test config.diagrams[1]["coupling"] == "c_pi1_f2"
end
