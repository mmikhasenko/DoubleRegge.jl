# const mπ = 0.13957; const mπ2 = mπ^2;
# const mη = 0.54782; const mη2 = mη^2;
# const mp = 0.938;   const mp2 = mp^2;
# const s0 = 19^2;

# const mπ = 0.13957018; const mπ2 = mπ^2;
# const mη = 0.54751; const mη2 = mη^2;
# const mp = 0.93827203;   const mp2 = mp^2;
# const s0 = DoubleRegge.mp2 + 2*DoubleRegge.mp*190;


@with_kw struct TwoParticleDiffraction
    mb::Float64 = -1.0
    mt::Float64 = -1.0
    mr::Float64 = -1.0
    m1::Float64 = -1.0
    m2::Float64 = -1.0
end

const compass_ηπ = let
    mπ = 0.13957018; mπ2 = mπ^2;
    mη = 0.54751; mη2 = mη^2;
    mp = 0.93827203;   mp2 = mp^2;
    TwoParticleDiffraction(mb = mπ, mt=mp, mr=mp, m1=mη, m2=mπ)
end
const compass_η′π = let
    mπ = 0.13957018; mπ2 = mπ^2;
    mη = 0.54751; mη2 = mη^2;
    mp = 0.93827203;   mp2 = mp^2;
    mη′ = 0.95778; mη′2 = mη′^2
    TwoParticleDiffraction(mb = mπ, mt=mp, mr=mp, m1=mη′, m2=mπ)
end

mutable struct setup
    system::TwoParticleDiffraction
    s0::Float64
end

const G = setup(TwoParticleDiffraction(), 0.0)

setsystem!(tpd::TwoParticleDiffraction, s0::T where T<:Real) = (G.system = tpd; G.s0 = s0)
#
const compass_Eb = 190
const compass_mb = 0.93827203
const compass_s0 = compass_mb^2 + 2*compass_mb*compass_Eb
# 
function setsystem!(s::Symbol)
    s == :compass_ηπ  && setsystem!(compass_ηπ, compass_s0)
    s == :compass_η′π && setsystem!(compass_η′π, compass_s0)
    return
end

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
    (k!=2) && error("only ^2 is defined")
    @unpack mb, mt, mr, m1, m2 = tpd
    return TwoParticleDiffraction²(mb^2, mt^2, mr^2, m1^2, m2^2)
end