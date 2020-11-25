using DrWatson
@quickactivate "DoubleRegge"
using TOML
using TypedTables
using QuadGK
using Cuba
using ForwardDiff
using Optim
using DelimitedFiles
using BenchmarkTools
using Cuba


function randkinvars(pars, setup = G)
    @unpack abst2min, abst2max, b, s1min, s1max = pars
    s = setup.s0
    # s1min = (setup.system.m1 + setup.system.m2)^2
    s1 = s1min + (s1max-s1min)*rand()
    cosθ = 2rand()-1
    ϕ = π*(2rand()-1)
    t2 = -abst2min+log(rand()/b)/b
    abs(t2) > abst2max && return randkinvars(pars, setup)
    # 
    A = exp(t2*b)
    ρ2 = sqrt(λ(s1, setup.system.m1^2, setup.system.m2^2))/s1
    invwt = A
    w = ρ2
    # 
    return (; s, s1, cosθ, ϕ, t2, w, invwt)
end

l3 = sixexchages[3][4]
l4 = sixexchages[4][4]
model3(vars) = modelDR(sixexchages[3][1], sixexchages[3][2], vars;
    η_forward=sixexchages[3][3],  α′=0.8)
model4(vars) = modelDR(sixexchages[4][1], sixexchages[4][2], vars;
    η_forward=sixexchages[4][3],  α′=0.8)
# 


#  _|_|_|  _|_|      _|_|_|  
#  _|    _|    _|  _|        
#  _|    _|    _|  _|        
#  _|    _|    _|    _|_|_|  


pars = (abst2min=0.1, abst2max=1.0, b=5, s1min=2.4^2, s1max=3.0^2)
pseudodata = Table([randkinvars(pars) for _ in 1:1_000_000])

# stephist(sqrt.(pseudodata.s1))
# stephist(sqrt.(pseudodata.s1), weights = pseudodata.w)
# stephist(pseudodata.t2, weights = pseudodata.w)

pseudodata_with_mw = Table(pseudodata,
    sπp = DoubleRegge.sπp.(pseudodata),
    sηp = DoubleRegge.sηp.(pseudodata),
    m3 = abs2.(model3.(pseudodata)) .* getproperty.(pseudodata, :w),
    m4 = abs2.(model4.(pseudodata)) .* getproperty.(pseudodata, :w))

p1 = let
    plot(xlab=L"m_{\eta\pi}\,(\mathrm{GeV})", ylab="normalized entries")
    stephist!( sqrt.(pseudodata_with_mw.s1), weights=pseudodata_with_mw.m3, norm=true, lab=l3,c=2)
    stephist!(sqrt.(pseudodata_with_mw.s1), weights=pseudodata_with_mw.m4, norm=true, lab=l4,c=3)
end
p2 = let bins=range(3,15,length=50)
    plot(xlab=L"m_{\eta p}\,(\mathrm{GeV})", ylab="normalized entries")
    stephist!(sqrt.(pseudodata_with_mw.sηp), weights=pseudodata_with_mw.w ./ pseudodata_with_mw.invwt, norm=true, lab="phase space", bins=bins)
    stephist!(sqrt.(pseudodata_with_mw.sηp), weights=pseudodata_with_mw.m3 ./ pseudodata_with_mw.invwt, norm=true, lab=l3, bins=bins)
    stephist!(sqrt.(pseudodata_with_mw.sηp), weights=pseudodata_with_mw.m4 ./ pseudodata_with_mw.invwt, norm=true, lab=l4, bins=bins)
end
p3 = let bins=range(14,18.5,length=50)
    plot(xlab=L"m_{\pi p}\,(\mathrm{GeV})", ylab="normalized entries", leg=:topleft)
    stephist!(sqrt.(pseudodata_with_mw.sπp), weights=pseudodata_with_mw.w ./ pseudodata_with_mw.invwt, norm=true, lab="phase space", bins=bins)
    stephist!(sqrt.(pseudodata_with_mw.sπp), weights=pseudodata_with_mw.m3 ./ pseudodata_with_mw.invwt, norm=true, lab=l3, bins=bins)
    stephist!(sqrt.(pseudodata_with_mw.sπp), weights=pseudodata_with_mw.m4 ./ pseudodata_with_mw.invwt, norm=true, lab=l4, bins=bins)
end

plot(p1,p2,p3, layout=grid(1,3), size=(1200,350))
savefig(joinpath("plots", "bottom_reggeon_mc.pdf"))
