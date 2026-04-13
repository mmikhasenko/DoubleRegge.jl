using DrWatson
@quickactivate "DoubleRegge"

using TOML
using DoubleRegge

const OUTPUT_FILE = joinpath(@__DIR__, "..", "test", "fixtures", "default_model_crosscheck.json")

function json_escape(text::AbstractString)
    replace(text,
        "\\" => "\\\\",
        "\"" => "\\\"",
        "\n" => "\\n",
        "\r" => "\\r",
        "\t" => "\\t")
end

function json_value(x::AbstractString)
    "\"" * json_escape(x) * "\""
end

json_value(x::Symbol) = json_value(String(x))
json_value(x::Bool) = x ? "true" : "false"
json_value(x::Real) = repr(Float64(x))
json_value(x::Nothing) = "null"
json_value(values::AbstractVector) = "[" * join(json_value.(values), ", ") * "]"

function json_value(dict::AbstractDict)
    entries = ["$(json_value(String(key))): $(json_value(value))" for (key, value) in pairs(dict)]
    "{\n" * join("  " .* entries, ",\n") * "\n}"
end

function model_fixture(tag, reaction_system, points)
    parsed = TOML.parsefile(joinpath(@__DIR__, "..", "data", "exp_pro", tag, "fit-results.toml"))
    settings = parsed["settings"]
    fit_results = parsed["fit_results"]
    model = build_model(
        sixexchages[settings["exchanges"]],
        settings["t2"],
        settings["scale_α"],
        reaction_system;
        s2shift = get(settings, "s2_shift", 0.0),
    )
    Dict(
        "tag" => tag,
        "system" => String(reaction_system.name),
        "fit_results_file" => joinpath("data", "exp_pro", tag, "fit-results.toml"),
        "t2" => settings["t2"],
        "points" => [
            Dict(
                "m" => point[:m],
                "cos_theta" => point[:cosθ],
                "phi" => point[:ϕ],
                "amplitude_re" => real(model(point[:m], point[:cosθ], point[:ϕ]; pars = fit_results["fit_minimizer"])),
                "amplitude_im" => imag(model(point[:m], point[:cosθ], point[:ϕ]; pars = fit_results["fit_minimizer"])),
            ) for point in points
        ],
    )
end

points = [
    (m = 2.45, cosθ = 0.2, ϕ = 0.3),
    (m = 2.75, cosθ = -0.55, ϕ = -1.1),
]

fixture = Dict(
    "schema_version" => 1,
    "models" => [
        model_fixture("etapi_a2Po-f2Po-a2f2-f2f2_opposite-sign", compass_ηπ, points),
        model_fixture("etappi_a2Po-a2f2-f2f2-PoPo_opposite-sign", compass_η′π, points),
    ],
)

open(OUTPUT_FILE, "w") do io
    write(io, json_value(fixture))
    write(io, "\n")
end

println("wrote $(OUTPUT_FILE)")
