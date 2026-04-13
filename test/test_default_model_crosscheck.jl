using TOML

function parse_json_fixture(path::AbstractString)
    text = read(path, String)
    index = firstindex(text)

    skip_ws(i) = begin
        while i <= lastindex(text) && isspace(text[i])
            i = nextind(text, i)
        end
        i
    end

    function parse_string(i)
        text[i] == '"' || error("expected string at index $i")
        i = nextind(text, i)
        buffer = IOBuffer()
        while true
            i > lastindex(text) && error("unterminated string")
            c = text[i]
            if c == '"'
                return String(take!(buffer)), nextind(text, i)
            elseif c == '\\'
                i = nextind(text, i)
                escaped = text[i]
                print(buffer,
                    escaped == '"'  ? '"'  :
                    escaped == '\\' ? '\\' :
                    escaped == '/'  ? '/'  :
                    escaped == 'b'  ? '\b' :
                    escaped == 'f'  ? '\f' :
                    escaped == 'n'  ? '\n' :
                    escaped == 'r'  ? '\r' :
                    escaped == 't'  ? '\t' :
                    error("unsupported escape sequence \\$escaped"))
                i = nextind(text, i)
            else
                print(buffer, c)
                i = nextind(text, i)
            end
        end
    end

    function parse_number(i)
        start = i
        allowed = Set(['+', '-', '.', 'e', 'E'])
        while i <= lastindex(text) && (isdigit(text[i]) || text[i] in allowed)
            i = nextind(text, i)
        end
        return parse(Float64, text[start:prevind(text, i)]), i
    end

    function parse_array(i)
        text[i] == '[' || error("expected array at index $i")
        i = skip_ws(nextind(text, i))
        values = Any[]
        if text[i] == ']'
            return values, nextind(text, i)
        end
        while true
            value, i = parse_value(i)
            push!(values, value)
            i = skip_ws(i)
            if text[i] == ']'
                return values, nextind(text, i)
            end
            text[i] == ',' || error("expected comma in array at index $i")
            i = skip_ws(nextind(text, i))
        end
    end

    function parse_object(i)
        text[i] == '{' || error("expected object at index $i")
        i = skip_ws(nextind(text, i))
        pairs = Dict{String, Any}()
        if text[i] == '}'
            return pairs, nextind(text, i)
        end
        while true
            key, i = parse_string(i)
            i = skip_ws(i)
            text[i] == ':' || error("expected colon after key at index $i")
            value, i = parse_value(skip_ws(nextind(text, i)))
            pairs[key] = value
            i = skip_ws(i)
            if text[i] == '}'
                return pairs, nextind(text, i)
            end
            text[i] == ',' || error("expected comma in object at index $i")
            i = skip_ws(nextind(text, i))
        end
    end

    function parse_value(i)
        i = skip_ws(i)
        c = text[i]
        c == '{' && return parse_object(i)
        c == '[' && return parse_array(i)
        c == '"' && return parse_string(i)
        c == 't' && return true, i + 4
        c == 'f' && return false, i + 5
        c == 'n' && return nothing, i + 4
        return parse_number(i)
    end

    value, index = parse_value(index)
    skip_ws(index) > lastindex(text) || error("unexpected trailing data in $path")
    return value
end

load_model_crosscheck_fixture() =
    parse_json_fixture(joinpath(@__DIR__, "fixtures", "default_model_crosscheck.json"))

function build_model_from_fixture_entry(entry)
    parsed = TOML.parsefile(joinpath(@__DIR__, "..", entry["fit_results_file"]))
    settings = parsed["settings"]
    fit_results = parsed["fit_results"]
    reaction_system = getproperty(DoubleRegge, Symbol(entry["system"]))
    model = build_model(
        sixexchages[settings["exchanges"]],
        settings["t2"],
        settings["scale_α"],
        reaction_system;
        s2shift = get(settings, "s2_shift", 0.0),
    )
    return model, reaction_system, fit_results["fit_minimizer"]
end

@testset "Default model JSON cross-check" begin
    fixture = load_model_crosscheck_fixture()
    @test fixture["schema_version"] == 1.0
    @test length(fixture["models"]) == 2

    for model_entry in fixture["models"]
        model, reaction_system, pars = build_model_from_fixture_entry(model_entry)
        for point in model_entry["points"]
            amplitude = model(point["m"], point["cos_theta"], point["phi"]; pars = pars)
            expected = point["amplitude_re"] + point["amplitude_im"] * im
            @test isapprox(amplitude, expected; rtol = 1e-12, atol = 1e-12)

            vars = (
                s = reaction_system.s0,
                s1 = point["m"]^2,
                cosθ = point["cos_theta"],
                ϕ = point["phi"],
                t2 = model_entry["t2"],
            )
            @test isfinite(abs2(modelDR(α_a2, α_ℙ, vars, reaction_system; η_forward = true)))
        end
    end
end
