using DrWatson
@quickactivate "DoubleRegge"

using JSON
using TOML
using DoubleRegge

const OUTPUT_FILE = joinpath(@__DIR__, "..", "test", "fixtures", "default_model_crosscheck.json")

function model_fixture(tag, reaction_system, points)
    parsed = TOML.parsefile(joinpath(@__DIR__, "..", "data", "exp_pro", tag, "fit-results.toml"))
    settings = parsed["settings"]
    fit_results = parsed["fit_results"]
    model = DoubleReggeModel(
        sixexchages[settings["exchanges"]],
        settings["t2"],
        settings["scale_α"],
        reaction_system,
        fit_results["fit_minimizer"];
        s2shift=get(settings, "s2_shift", 0.0),
    )
    Dict(
        "tag" => tag,
        "system" => String(reaction_system.name),
        "fit_results_file" => joinpath("data", "exp_pro", tag, "fit-results.toml"),
        "t2" => model.t2,
        "points" => [
            Dict(
                "m" => point[:m],
                "cos_theta" => point[:cosθ],
                "phi" => point[:ϕ],
                "amplitude_re" => real(amplitude(model, point[:m], point[:cosθ], point[:ϕ])),
                "amplitude_im" => imag(amplitude(model, point[:m], point[:cosθ], point[:ϕ])),
            ) for point in points
        ],
    )
end

points = [
    (m=2.45, cosθ=0.2, ϕ=0.3),
    (m=2.75, cosθ=-0.55, ϕ=-1.1),
]

fixture = Dict(
    "schema_version" => 1.0,
    "models" => [
        model_fixture("etapi_a2Po-f2Po-a2f2-f2f2_opposite-sign", compass_ηπ, points),
        model_fixture("etappi_a2Po-a2f2-f2f2-PoPo_opposite-sign", compass_η′π, points),
    ],
)

open(OUTPUT_FILE, "w") do io
    JSON.print(io, fixture, 2)
    write(io, '\n')
end

println("wrote $(OUTPUT_FILE)")
