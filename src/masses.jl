
const standard_masses = (mπ = 0.13957018, mp = 0.93827203, mη = 0.54751, mη′ = 0.95778)

@with_kw struct TwoParticleDiffraction
    mb::Float64 = -1.0
    mt::Float64 = -1.0
    mr::Float64 = -1.0
    m1::Float64 = -1.0
    m2::Float64 = -1.0
end

const channel_compass_ηπ = let
    @unpack mπ, mη, mp = standard_masses
    TwoParticleDiffraction(mb = mπ, mt = mp, mr = mp, m1 = mη, m2 = mπ)
end
const channel_compass_η′π = let
    @unpack mπ, mp, mη′ = standard_masses
    TwoParticleDiffraction(mb = mπ, mt = mp, mr = mp, m1 = mη′, m2 = mπ)
end

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

struct ReactionSystem
    name::Symbol
    channel::TwoParticleDiffraction
    s0::Float64
    description::Vector{Pair{Tuple{Int,Int},Tuple{String,String}}}
    LMs::Vector{Tuple{Int,Int}}
end

ReactionSystem(
    channel::TwoParticleDiffraction,
    s0::Real;
    name::Symbol = :custom,
    description::Vector{Pair{Tuple{Int,Int},Tuple{String,String}}} = Pair{
        Tuple{Int,Int},
        Tuple{String,String},
    }[],
) = ReactionSystem(name, channel, Float64(s0), description, getindex.(description, 1))

const compass_Eb = 190
const compass_mb = 0.93827203
const compass_s0 = compass_mb^2 + 2 * compass_mb * compass_Eb

const compass_ηπ = ReactionSystem(
    channel_compass_ηπ,
    compass_s0;
    name = :compass_ηπ,
    description = description_ηπ,
)
const compass_η′π = ReactionSystem(
    channel_compass_η′π,
    compass_s0;
    name = :compass_η′π,
    description = description_η′π,
)

diffraction(system::ReactionSystem) = system.channel
diffraction(system::TwoParticleDiffraction) = system

#    _|_|_|  _|    _|    _|_|_|    _|_|_|  _|  _|_|  
#  _|_|      _|    _|  _|    _|  _|    _|  _|_|      
#      _|_|  _|    _|  _|    _|  _|    _|  _|        
#  _|_|_|      _|_|_|    _|_|_|    _|_|_|  _|        
#                            _|                      
#                        _|_|                        

struct TwoParticleDiffraction²
    mb²::Float64
    mt²::Float64
    mr²::Float64
    m1²::Float64
    m2²::Float64
end

import Base: ^
function ^(tpd::TwoParticleDiffraction, k)
    (k != 2) && error("only ^2 is defined")
    @unpack mb, mt, mr, m1, m2 = tpd
    return TwoParticleDiffraction²(mb^2, mt^2, mr^2, m1^2, m2^2)
end
