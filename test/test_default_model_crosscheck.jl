using JSON
using TOML

load_model_crosscheck_fixture() =
    JSON.parsefile(joinpath(@__DIR__, "fixtures", "default_model_crosscheck.json"))

function build_model_from_fixture_entry(entry)
    parsed = TOML.parsefile(joinpath(@__DIR__, "..", entry["fit_results_file"]))
    settings = parsed["settings"]
    fit_results = parsed["fit_results"]
    reaction_system = getproperty(DoubleRegge, Symbol(entry["system"]))
    model = DoubleReggeModel(
        six_exchanges[settings["exchanges"]],
        settings["t2"],
        settings["scale_α"],
        reaction_system,
        fit_results["fit_minimizer"];
        s2shift = get(settings, "s2_shift", 0.0),
    )
    return model, reaction_system
end

@testset "Default model JSON cross-check" begin
    fixture = load_model_crosscheck_fixture()
    @test fixture["schema_version"] == 1.0
    @test length(fixture["models"]) == 2

    for model_entry in fixture["models"]
        model, reaction_system = build_model_from_fixture_entry(model_entry)
        for point in model_entry["points"]
            amp = amplitude(model, point["m"], point["cos_theta"], point["phi"])
            expected = point["amplitude_re"] + point["amplitude_im"] * im
            @test isapprox(amp, expected; rtol = 1e-12, atol = 1e-12)

            vars = (
                s = reaction_system.s0,
                s1 = point["m"]^2,
                cosθ = point["cos_theta"],
                ϕ = point["phi"],
                t2 = model.t2,
            )
            @test isfinite(
                abs2(modelDR(α_a2, α_ℙ, vars, reaction_system; η_forward = true)),
            )
        end
    end
end
