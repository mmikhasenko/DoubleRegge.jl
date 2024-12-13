# #
# slice(data, property, bin; ind = 2) = [v[bin, ind] for v in getproperty.(data, property)]
# slice(data, bin) = getindex.(data, bin)
# #

strip_errors(_Iϕ::TwoBodyPartialWaveIϕs{Tuple{Measurement{Float64}, Measurement{Float64}}}) =
    TwoBodyPartialWaves(
        _Iϕ.LMs,
        NamedTuple{(:I, :ϕ)}.(zip(
            getproperty.(getindex.(_Iϕ.PWs, :I), :val),
            getproperty.(getindex.(_Iϕ.PWs, :ϕ), :val),
        )) |> SVector,
    )

function bootstrap_band(_Iϕ; nPoints::Int = 101, nSamples::Int = 1000)
    # x ± δx
    I_with_err = getindex.(_Iϕ.PWs, :I)
    ϕ_with_err = getindex.(_Iϕ.PWs, :ϕ)
    # # val and err
    i, δi = getproperty.(I_with_err, :val), getproperty.(I_with_err, :err)
    ϕ, δϕ = getproperty.(ϕ_with_err, :val), getproperty.(ϕ_with_err, :err)
    # computation
    amps = changerepresentation(strip_errors(_Iϕ))
    cosθv = range(-1, 1, nPoints)
    calv = dNdcosθ.(cosθv, Ref(amps))
    # bootstrap
    bstp = map(1:nSamples) do _
        i′ = i .+ randn(length(i)) .* δi
        i′ = map(v -> v < 0 ? 0.0 : v, i′)
        ϕ′ = ϕ + randn(length(i)) .* δϕ
        nt_Iϕ′ = NamedTuple{(:I, :ϕ)}.(zip(i′, ϕ′))
        Iϕ′ = TwoBodyPartialWaves(_Iϕ.LMs, nt_Iϕ′ |> SVector)
        A′ = changerepresentation(Iϕ′)
        dNdcosθ.(cosθv, Ref(A′))
    end
    qls = [quantile(getindex.(bstp, i), [0.16, 0.84]) for i in 1:length(bstp[1])]
    return getindex.(qls, 1), getindex.(qls, 2)
end


# data = read_data(settings["pathtodata"], description_ηπ);
# @test typeof(strip_errors(data.Iϕ[1])) <: TwoBodyPartialWaveIϕs{Tuple{Float64, Float64}}
