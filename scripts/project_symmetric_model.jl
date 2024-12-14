using DrWatson
@quickactivate "DoubleRegge"
using TOML
using TypedTables
using QuadGK
using Cuba
# 
using Plots
import Plots.PlotMeasures.mm
theme(:wong2; size = (500, 350), bottom_margin = 5mm)
# 
using DoubleRegge
setsystem!(:compass_ηπ)
# 
using Statistics
using LinearAlgebra

#                            _|            
#    _|_|_|    _|_|      _|_|_|    _|_|    
#  _|        _|    _|  _|    _|  _|_|_|_|  
#  _|        _|    _|  _|    _|  _|        
#    _|_|_|    _|_|      _|_|_|    _|_|_|  

const pfr = [0.5, -0.5]
# 
setsystem!(:compass_ηπ)

# fit
const exchanges = sixexchages[[1, 3]]
const model = build_model(exchanges, -0.2, 0.8)
fixed_model(m, cosθ, ϕ; pars = pfr) = model(m, cosθ, ϕ; pars = pars)
fixed_model_sqrtq(m, cosθ, ϕ; pars = pfr) = fixed_model(m, cosθ, ϕ; pars = pars) * sqrt(q(m))
intensity(m, cosθ, ϕ; pars = pfr) = abs2(fixed_model_sqrtq(m, cosθ, ϕ; pars = pars))
# 
function pw_project_fixed_model(m, used_LMs)
    amplitude(cosθ, ϕ) = fixed_model_sqrtq(m, cosθ, ϕ) * sqrt(q(m))
    @show used_LMs
    return map(used_LMs) do (L, M)
        pw_project(amplitude, L, M)
    end
end

# 
model_integral(m)          = integrate_dcosθdϕ((cosθ, ϕ) -> abs2(fixed_model(m, cosθ, ϕ)))[1] * q(m)
model_integral_forward(m)  = integrate_dcosθdϕ((cosθ, ϕ) -> abs2(fixed_model(m, cosθ, ϕ)), (0, 1))[1] * q(m)
model_integral_backward(m) = integrate_dcosθdϕ((cosθ, ϕ) -> abs2(fixed_model(m, cosθ, ϕ)), (-1, 0))[1] * q(m)

# data
const data_folder = "data/exp_raw/PLB_shifted"
data = read_data(data_folder, description_ηπ)
# const LMs = compass_ηπ_LMs
# data = Table(x_IδI_ϕδϕ_compass_ηπ(settings["pathtodata"]))
# amplitudes = [sqrt.(is) .* cis.(ϕs) for (is,ϕs) in zip(data.I, data.ϕ)]
# data = Table(data, amps=amplitudes)
# # fit range
# fitrangemap = map(x->inlims(x.x, settings["fitrange"]), data)
# fitdata = data[fitrangemap]
# # plot 
plotmap = map(x -> inlims(x.x, (2.2, 3.0)), data)
plotdata = data[plotmap]

# intensity
intensity_in_bins = model_integral.(plotdata.x)

# asymmetry
intensity_forward_in_bins = model_integral_forward.(plotdata.x)
intensity_backward_in_bins = model_integral_backward.(plotdata.x)
data_forward_in_bins = map(a -> quadgk(cosθ -> dNdcosθ(cosθ, a), 0, 1)[1], plotdata.amps)
data_backward_in_bins = map(a -> quadgk(cosθ -> dNdcosθ(cosθ, a), -1, 0)[1], plotdata.amps)
# 

# # projections
const LMs = getindex.(description_ηπ, 1)
pw_projections = map(m -> pw_project_fixed_model(m, LMs), data.x[plotmap])
pw_intensities = map(x -> abs2.(x), pw_projections)
#
tolab(LM) = "$(LM.L)$(LM.M)"

# Odd and Even waves
fHeigher = (intensity_in_bins .- sum.(pw_intensities)) ./ intensity_in_bins
filtodd = isodd.(getindex.(LMs, 1))
fOdd = map(x -> sum(x .* filtodd), pw_intensities) ./ intensity_in_bins
fEven = map(x -> sum(x .* iszero.(filtodd)), pw_intensities) ./ intensity_in_bins
# 

tot_compass = map(plotdata.Iϕ) do _Iϕ
    sum(x -> x.I, _Iϕ.PWs)
end
fOdd_compass = map(plotdata.Iϕ) do _Iϕ
    sum(getproperty.(_Iϕ.PWs, :I) .* filtodd)
end ./ tot_compass
fEven_compass = map(plotdata.Iϕ) do _Iϕ
    sum(getproperty.(_Iϕ.PWs, :I) .* iszero.(filtodd))
end ./ tot_compass

data_pw_sums = [sum(getindex.(fitdata.I, i)) for (i, LM) in enumerate(LMs)]
data_pw_sums = [data_pw_sums..., 0]

b1 = bar(data_pw_sums, yaxis = nothing,
    xticks = (1:8, [tolab.(LMs)..., "higher"]), xlab = "LM", ylab = "intensity", lab = "")
# 
model_pw_sums = [sum(getindex.(pw_intensities[6:end], i)) for (i, LM) in enumerate(LMs)]
model_pw_sums = [model_pw_sums..., sum(intensity_in_bins .- sum.(pw_intensities))]

b2 = bar(model_pw_sums, yaxis = nothing,
    xticks = (1:8, [tolab.(LMs)..., "higher"]), xlab = "LM", ylab = "intensity", lab = "")
#
plot(b1, b2, layout = grid(2, 1), size = (400, 600), link = :x)
savefig(joinpath("plots", "pws_symmetric_model.pdf"))

let
    plot(ylab = "fraction", xlab = "m(ηπ) (GeV)", size = (500, 350), title = settings["tag"])
    plot!(plotdata.x, fHeigher, lab = "Higher waves L > 6", lw = 2)
    plot!(plotdata.x, fEven, lab = "Even waves L ≤ 6", lw = 2)
    plot!(plotdata.x, fOdd, lab = "Odd waves L ≤ 6", lw = 2)
    scatter!(plotdata.x, xerr = (plotdata.x[2] - plotdata.x[1]) / 2, fEven_compass, lab = "", mc = 2, ms = 3)
    scatter!(plotdata.x, xerr = (plotdata.x[2] - plotdata.x[1]) / 2, fOdd_compass, lab = "", mc = 3, ms = 3)
    plot!(ylims = (0, 1), leg = :left)
    vspan!(fitdata.x[[1, end]], lab = "", α = 0.1, seriescolor = 7)
end
savefig(
    joinpath("data", "exp_pro", settings["tag"],
        "odd-and-even_$(settings["tag"])_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf"))
