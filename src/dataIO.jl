struct TwoBodyPartialWaves{N, T}
    LMs::SVector{N, Tuple{Int, Int}}
    PWs::SVector{N, T}
end

const TwoBodyPartialWaveAs{N, T} = TwoBodyPartialWaves{N, T} where N where T <: Number
const TwoBodyPartialWaveIϕs{N, V} = TwoBodyPartialWaves{N, NamedTuple{(:I, :ϕ), V}} where N where V

TwoBodyPartialWaves(LMs::Vector, PWs::Vector) =
    (N = length(LMs); TwoBodyPartialWaves(SVector{N}(LMs), SVector{N}(PWs)))
#
changerepresentation(expansion::TwoBodyPartialWaveIϕs) =
    TwoBodyPartialWaves(expansion.LMs, map(x -> amplitude(x.I, x.ϕ), expansion.PWs))
changerepresentation(expansion::TwoBodyPartialWaveAs; iref) =
    TwoBodyPartialWaves(expansion.LMs, Iϕ.(expansion.PWs; ref = expansion.PWs[iref]))
# 

strip_errors(_Iϕ::TwoBodyPartialWaveIϕs{Tuple{Measurement{Float64}, Measurement{Float64}}}) =
    TwoBodyPartialWaves(
        _Iϕ.LMs,
        NamedTuple{(:I, :ϕ)}.(zip(
            getproperty.(getindex.(_Iϕ.PWs, :I), :val),
            getproperty.(getindex.(_Iϕ.PWs, :ϕ), :val),
        )) |> SVector,
    )

function read_data(path2data, description)
    LMs = getindex.(description, 1)
    #
    nBins, _ = size(readdlm(joinpath(path2data, description[1][2][1])))
    #
    Is, ϕs, xs = [], [], []
    for (LM, (filename_I, filename_ϕ)) in description
        # 
        table_I = readdlm(joinpath(path2data, filename_I))
        x, I, δI = table_I[:, 1], table_I[:, 2], table_I[:, 3]
        #
        table_ϕ = filename_ϕ != "nothing" ? readdlm(joinpath(path2data, filename_ϕ)) : zeros(nBins, 3)
        table_ϕ .*= π / 180
        ϕ, δϕ = table_ϕ[:, 2], table_ϕ[:, 3]
        #
        push!(Is, I .± δI)
        push!(ϕs, ϕ .± δϕ)
        push!(xs, x)
    end
    x = xs[1]
    Im, ϕm = hcat(Is...), hcat(ϕs...)
    #
    Iϕ = map(1:nBins) do i
        TwoBodyPartialWaves(LMs, NamedTuple{(:I, :ϕ)}.(zip(Im[i, :], ϕm[i, :])))
    end
    #
    A = changerepresentation.(strip_errors.(Iϕ))
    # 
    Table(x = x, Iϕ = Iϕ, amps = A)
end


rand(v::Measurement{T} where T) = v.val + randn() * v.err
function rand(expansion::TwoBodyPartialWaves{N, NamedTuple{(:I, :ϕ), Tuple{Measurement{V}, Measurement{V}}}} where N where V)
    I′ = rand.(expansion.PWs .. :I)
    ϕ′ = rand.(expansion.PWs .. :ϕ)
    TwoBodyPartialWaves(expansion.LMs, NamedTuple{(:I, :ϕ)}.(zip(I′, ϕ′)))
end

# function randA(data)
#     randI = [abs.(I .+ randn(length(δI)) .* δI) for (I, δI) in zip(data.I, data.δI)]
#     randϕ = [abs.(ϕ .+ randn(length(δϕ)) .* δϕ) for (ϕ, δϕ) in zip(data.ϕ, data.δϕ)]
#     randA = [sqrt.(is) .* cis.(ϕs) for (is,ϕs) in zip(randI,randϕ)]
#     return randA
# end



# 
const description_ηπ = [
    (1, 1) => ("EtaPi-1mp.txt", "EtaPi-Ph1mp.txt"),
    (2, 1) => ("EtaPi-2pp.txt", "nothing"),
    (2, 2) => ("EtaPi-2ppM2.txt", "EtaPi-Ph2ppM2.txt"),
    (3, 1) => ("EtaPi-3mp.txt", "EtaPi-Ph3mp.txt"),
    (4, 1) => ("EtaPi-4pp.txt", "EtaPi-Ph4pp.txt"),
    (5, 1) => ("EtaPi-5mp.txt", "EtaPi-Ph5mp.txt"),
    (6, 1) => ("EtaPi-6pp.txt", "EtaPi-Ph6pp.txt"),
]
const description_η′π = [
    (1, 1) => ("EtapPi-1mp.txt", "EtapPi-Ph1mp.txt"),
    (2, 1) => ("EtapPi-2pp.txt", "nothing"),
    (3, 1) => ("EtapPi-3mp.txt", "EtapPi-Ph3mp.txt"),
    (4, 1) => ("EtapPi-4pp.txt", "EtapPi-Ph4pp.txt"),
    (5, 1) => ("EtapPi-5mp.txt", "EtapPi-Ph5mp.txt"),
    (6, 1) => ("EtapPi-6pp.txt", "EtapPi-Ph6pp.txt"),
]
