using DrWatson
@quickactivate "DoubleRegge"
using TOML
using TypedTables
using Cuba
using Optim
using DelimitedFiles

using Plots
using LaTeXStrings
theme(:wong2)
import Plots.PlotMeasures.mm

using DoubleRegge

settings = Dict(
    "pathtodata" => joinpath("data","exp_raw","PLB_shifted"),
    "fitrange" => [2.4, 2.96],
    "s2_shift" => 30.0,
    "scale_α" => 0.8
)
setsystem!(:compass_ηπ)

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
    η_forward=sixexchages[3][3],  α′=settings["scale_α"])
model4(vars) = modelDR(sixexchages[4][1], sixexchages[4][2], vars;
    η_forward=sixexchages[4][3],  α′=settings["scale_α"])
# 
model5(vars) = modelDR(sixexchages[3][1], sixexchages[3][2], vars;
    η_forward=sixexchages[3][3],  α′=settings["scale_α"], s2shift=settings["s2_shift"])
#

const LMs = compass_ηπ_LMs
const data = Table(x_IδI_ϕδϕ_compass_ηπ(settings["pathtodata"]))
const amplitudes = [sqrt.(is) .* cis.(ϕs) for (is,ϕs) in zip(data.I, data.ϕ)]
# range
fitrangemap = map(x->inlims(x.x, settings["fitrange"]), data)
const fitdata = Table(data[fitrangemap], amps = amplitudes[fitrangemap])

function guess_bin(mηπ)
    dx = fitdata.x[2]-fitdata.x[1]
    return Int(div((mηπ - (fitdata.x[1]-dx/2)), dx)) + 1
end

function intensity_data(vars)
    @unpack s1, cosθ, ϕ = vars
    mηπ = sqrt(s1)
    # 
    bin = guess_bin(mηπ)
    bin > length(fitdata) && error("$(mηπ) > $(fitdata.x[end])")
    return abs2(recamp(cosθ, ϕ, fitdata.amps[bin], LMs))
end

model6(vars) = vars.cosθ > 0 ? intensity_data(vars) : 0.0 # forward
model7(vars) = vars.cosθ < 0 ? intensity_data(vars) : 0.0 # backward

#  _|_|_|  _|_|      _|_|_|  
#  _|    _|    _|  _|        
#  _|    _|    _|  _|        
#  _|    _|    _|    _|_|_|  


pars = (abst2min=0.1, abst2max=1.0, b=5,
    s1min=settings["fitrange"][1]^2,
    s1max=settings["fitrange"][2]^2)
pseudodata = Table([randkinvars(pars) for _ in 1:1_000_000])

# stephist(sqrt.(pseudodata.s1))
# stephist(sqrt.(pseudodata.s1), weights = pseudodata.w)
# stephist(pseudodata.t2, weights = pseudodata.w)

pseudodata_with_mw = Table(pseudodata,
    sπp = DoubleRegge.sπp.(pseudodata),
    sηp = DoubleRegge.sηp.(pseudodata),
    m3 = abs2.(model3.(pseudodata)) .* getproperty.(pseudodata, :w),
    m4 = abs2.(model4.(pseudodata)) .* getproperty.(pseudodata, :w),
    m5 = abs2.(model5.(pseudodata)) .* getproperty.(pseudodata, :w),
    m7 = abs2.(model7.(pseudodata)) .* getproperty.(pseudodata, :w))

p1 = let
    plot(xlab=L"m_{\eta\pi}\,(\mathrm{GeV})", ylab="normalized entries")
    scatterhist!(sqrt.(pseudodata_with_mw.s1), weights=pseudodata_with_mw.m7, norm=true, lab="data",c=1)
    stephist!(sqrt.(pseudodata_with_mw.s1), weights=pseudodata_with_mw.m3, norm=true, lab=l3,c=2)
    stephist!(sqrt.(pseudodata_with_mw.s1), weights=pseudodata_with_mw.m4, norm=true, lab=l4,c=3)
    stephist!(sqrt.(pseudodata_with_mw.s1), weights=pseudodata_with_mw.m5, norm=true, lab=l3*"+shift",c=4)
end
p2 = let bins=range(3,15,length=50)
    plot(xlab=L"m_{\eta p}\,(\mathrm{GeV})", ylab="normalized entries")
    scatterhist!(sqrt.(pseudodata_with_mw.sηp), weights=pseudodata_with_mw.m7 ./ pseudodata_with_mw.invwt, norm=true, lab="data backward", bins=bins)
    # stephist!(sqrt.(pseudodata_with_mw.sηp), weights=pseudodata_with_mw.w ./ pseudodata_with_mw.invwt, norm=true, lab="phase space", bins=bins)
    stephist!(sqrt.(pseudodata_with_mw.sηp), weights=pseudodata_with_mw.m3 ./ pseudodata_with_mw.invwt, norm=true, lab=l3, bins=bins)
    stephist!(sqrt.(pseudodata_with_mw.sηp), weights=pseudodata_with_mw.m4 ./ pseudodata_with_mw.invwt, norm=true, lab=l4, bins=bins)
    stephist!(sqrt.(pseudodata_with_mw.sηp), weights=pseudodata_with_mw.m5 ./ pseudodata_with_mw.invwt, norm=true, lab=l3*"+shift", bins=bins)
end
p3 = let bins=range(14,18.5,length=50)
    plot(xlab=L"m_{\pi p}\,(\mathrm{GeV})", ylab="normalized entries", leg=:topleft)
    scatterhist!(sqrt.(pseudodata_with_mw.sπp), weights=pseudodata_with_mw.m7 ./ pseudodata_with_mw.invwt, norm=true, lab="data", bins=bins)
    # stephist!(sqrt.(pseudodata_with_mw.sπp), weights=pseudodata_with_mw.w ./ pseudodata_with_mw.invwt, norm=true, lab="phase space", bins=bins)
    stephist!(sqrt.(pseudodata_with_mw.sπp), weights=pseudodata_with_mw.m3 ./ pseudodata_with_mw.invwt, norm=true, lab=l3, bins=bins)
    stephist!(sqrt.(pseudodata_with_mw.sπp), weights=pseudodata_with_mw.m4 ./ pseudodata_with_mw.invwt, norm=true, lab=l4, bins=bins)
    stephist!(sqrt.(pseudodata_with_mw.sπp), weights=pseudodata_with_mw.m5 ./ pseudodata_with_mw.invwt, norm=true, lab=l3*"+shift", bins=bins)
end

plot(p1,p2,p3, layout=grid(1,3), size=(1200,350), bottom_margin=5mm, left_margin=5mm)
savefig(joinpath("plots", "bottom_reggeon_mc.pdf"))
