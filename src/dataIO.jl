struct TwoBodyPartialWaves{N,T}
    LMs::SVector{N,Tuple{Int,Int}}
    PWs::SVector{N,T}
end

TwoBodyPartialWaves(LMs::Vector, PWs::Vector) =
    (N=length(LMs); TwoBodyPartialWaves(SVector{N}(LMs), SVector{N}(PWs)))


function read_data(pathtodata, description)
    LMs = getindex.(description,1)
    #
    Nbins,_ = size(readdlm(joinpath(pathtodata, description[1][2][1])))
    #
    Is, ϕs, xs = [], [], []
    for (LM, (filename_I,filename_ϕ)) in description
        # 
        table_I = readdlm(joinpath(pathtodata, filename_I))
        x,I,δI = table_I[:,1], table_I[:,2], table_I[:,3]
        #
        table_ϕ = filename_ϕ != "nothing" ? readdlm(joinpath(pathtodata, filename_ϕ)) : zeros(Nbins,3)
        table_ϕ .*= π / 180 
        ϕ,δϕ = table_ϕ[:,2], table_ϕ[:,3]
        #
        push!(Is, I .± δI)
        push!(ϕs, ϕ .± δϕ)
        push!(xs, x)
    end
    x = xs[1]
    Im, ϕm = hcat(Is...), hcat(ϕs...)
    #
    Iϕ = [TwoBodyPartialWaves(LMs, NamedTuple{(:I,:ϕ)}.(zip(Im[i,:], ϕm[i,:])))
        for i in 1:Nbins]
    A = [TwoBodyPartialWaves(expansion.LMs,
                             map(x->amplitude(x.I.val, x.ϕ.val), expansion.PWs))
            for expansion in Iϕ]
    Table(x=x, Iϕ = Iϕ, amps = A)
end

typeof((z=1±2,))


rand(v::Measurement{T} where T) = v.val + randn()*v.err
function rand(expansion::TwoBodyPartialWaves{N, NamedTuple{(:I,:ϕ),Tuple{Measurement{V},Measurement{V}}}} where N where V)
    I′ = rand.(expansion.PWs..:I)
    ϕ′ = rand.(expansion.PWs..:ϕ)
    TwoBodyPartialWaves(expansion.LMs, NamedTuple{(:I,:ϕ)}.(zip(I′,ϕ′)))
end
# function randA(data)
#     randI = [abs.(I .+ randn(length(δI)) .* δI) for (I, δI) in zip(data.I, data.δI)]
#     randϕ = [abs.(ϕ .+ randn(length(δϕ)) .* δϕ) for (ϕ, δϕ) in zip(data.ϕ, data.δϕ)]
#     randA = [sqrt.(is) .* cis.(ϕs) for (is,ϕs) in zip(randI,randϕ)]
#     return randA
# end



# 
const description_ηπ = [
    (1,1) => ("EtaPi-1mp.txt",   "EtaPi-Ph1mp.txt"),
    (2,1) => ("EtaPi-2pp.txt",   "nothing"),
    (2,2) => ("EtaPi-2ppM2.txt", "EtaPi-Ph2ppM2.txt"),
    (3,1) => ("EtaPi-3mp.txt",   "EtaPi-Ph3mp.txt"),
    (4,1) => ("EtaPi-4pp.txt",   "EtaPi-Ph4pp.txt"),
    (5,1) => ("EtaPi-5mp.txt",   "EtaPi-Ph5mp.txt"),
    (6,1) => ("EtaPi-6pp.txt",   "EtaPi-Ph6pp.txt")
]
const description_η′π = [
    (1,1) => ("EtapPi-1mp.txt",   "EtapPi-Ph1mp.txt"),
    (2,1) => ("EtapPi-2pp.txt",   "nothing"),
    (3,1) => ("EtapPi-3mp.txt",   "EtapPi-Ph3mp.txt"),
    (4,1) => ("EtapPi-4pp.txt",   "EtapPi-Ph4pp.txt"),
    (5,1) => ("EtapPi-5mp.txt",   "EtapPi-Ph5mp.txt"),
    (6,1) => ("EtapPi-6pp.txt",   "EtapPi-Ph6pp.txt")
]
