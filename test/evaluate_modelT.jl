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
    exchange = TReggeExchange(TVertex(α_top, 0.0, 1.0), TVertex(α_bot, 0.0, 1.0), true, "a2/ℙ")
    @test exchange.top.α === α_top
    @test exchange.bot.α === α_bot
    @test exchange.η_forward == true
    @test exchange.label == "a2/ℙ"
    gj = KinematicsGJ(ηπ_system_T.s0, 4.0^2, 0.4, π / 4, -0.45)
    amp = modelDR(exchange, gj, ηπ_system_T; α′ = 0.85)
    @test isfinite(real(amp))
    @test isfinite(imag(amp))
end

@testset "modelT KinematicsM ≡ KinematicsGJ path" begin
    α_top = trajectory(0.917, 0.44)
    α_bot = trajectory(0.25, 1.08)
    gj = KinematicsGJ(ηπ_system_T.s0, 4.0^2, 0.4, π / 4, -0.45)
    kin = KinematicsM(gj, ηπ_system_T)

    exchange = TReggeExchange(TVertex(α_top, 0.0, 1.0), TVertex(α_bot, 0.0, 1.0), true, "a2/ℙ")
    direct   = modelDR(exchange, kin;  α′ = 0.85)
    via_gj   = modelDR(exchange, gj, ηπ_system_T; α′ = 0.85)
    @test direct ≈ via_gj

    # form factors multiply in at both vertices
    b_top, b_bot = -0.5, 0.3
    ex_ff = TReggeExchange(
        TVertex(α_top, b_top, 1.0),
        TVertex(α_bot, b_bot, 1.0),
        true,
        "a2/ℙ-ff",
    )
    amp_ff = modelDR(ex_ff, kin; α′ = 0.85)
    @test amp_ff ≈ direct * form_factor(ex_ff.top, kin.t1) * form_factor(ex_ff.bot, kin.t2)
end

@testset "modelT KinematicsM direct path" begin
    parsed = TOML.parsefile(joinpath(@__DIR__, "..", "data", "exp_pro", "my_model", "model-settings.toml"))
    config = load_modelT_config(parsed)
    gj = KinematicsGJ(ηπ_system_T.s0, 2.32^2, cos(π / 3), π / 4, -0.2)
    kin = KinematicsM(gj, ηπ_system_T)
    exchange = config.model.exchanges[3]
    @test modelDR(exchange, gj, ηπ_system_T; α′ = 0.8) ≈
          modelDR(exchange, kin; α′ = 0.8)
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

    # all three amplitude entry points agree
    gj  = KinematicsGJ(ηπ_system_T.s0, 2.3^2, 0.4, π / 5, model.t2)
    kin = KinematicsM(gj, ηπ_system_T)
    @test amplitude(model, gj)  ≈ amp
    @test amplitude(model, kin) ≈ amp
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
