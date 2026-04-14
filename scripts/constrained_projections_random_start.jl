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

# 
using Plots
import Plots.PlotMeasures.mm
theme(:wong2; size = (500, 350), bottom_margin = 5mm)
# 

using DoubleRegge
const reaction_system = compass_ηπ


# settings_file = joinpath("data", "exp_pro","fit-results_bottom-Po_Np=3.toml")
settings_file = joinpath(
    "data",
    "exp_pro",
    "a2Po-f2f2-PoPo_opposite-sign",
    "fit-results_a2Po-f2f2-PoPo_opposite-sign_Np=3_alpha=0.8.toml",
)
! isfile(settings_file) && error("no file")
# 
parsed = TOML.parsefile(settings_file)
settings = parsed["settings"]
fit_results = parsed["fit_results"]
# 
constrained_pw_projection_fixed_model(m, init_pars) = constrained_pw_projection(
    (cosθ, ϕ)->intensity(m, cosθ, ϕ),
    init_pars,
    reaction_system.LMs,
)

# fit
const exchanges = six_exchanges[settings["exchanges"]]
const model = DoubleReggeModel(
    exchanges,
    settings["t2"],
    settings["scale_α"],
    reaction_system,
    fit_results["fit_minimizer"];
    s2shift = get(settings, "s2_shift", 0.0),
)
const fixed_pars = model.pars
fixed_model(m, cosθ, ϕ) = amplitude(model, m, cosθ, ϕ)
intensity(m, cosθ, ϕ) = abs2(fixed_model(m, cosθ, ϕ))*q(m, reaction_system)

function pw_project_fixed_model(m)
    amplitude(cosθ, ϕ) = fixed_model(m, cosθ, ϕ)*sqrt(q(m, reaction_system))
    pws = [pw_project(amplitude, L, M) for (L, M) in LMs]
    return pws
end
#

# data
const LMs = reaction_system.LMs
data = read_data(settings["pathtodata"], reaction_system)
# fit range
fitrangemap = map(x->inlims(x.x, settings["fitrange"]), data)
fitdata = data[fitrangemap]
# plot 
plotmap = map(x->inlims(x.x, (2.4, 3.0)), data)
plotdata = data[plotmap]
# 

struct samplePWA{N}
    amplitude::Function
    LMs::SVector{N}
    # 
    partialwaves::SVector{N,Complex{Float64}}
    samples::Vector
    # 
    function samplePWA(amplitude, LMs)
        pws = [pw_project(amplitude, L, M) for (L, M) in LMs]
        return new{length(LMs)}(amplitude, LMs, pws, [])
    end
end
Npars(s::samplePWA{N}) where {N} = N

function randcomplexwithnorm(Npars, integral)
    vals = 2*rand(Complex{Float64}, Npars) .- (1+1im)
    vals .*= integral / sum(abs2, vals)
end

function constrained_pw_projection_fixed_model!(s::samplePWA, N)
    integral = integrate_dcosθdϕ(abs2 ∘ s.amplitude)
    #
    for _ = 1:N
        initial_pars = randcomplexwithnorm(Npars(s), integral)
        fr = constrained_pw_projection(abs2 ∘ s.amplitude, initial_pars, s.LMs)
        push!(s.samples, fr)
    end
end

# create the structures
getsamplePWA(m) =
    samplePWA((cosθ, ϕ)->amplitude(model, m, cosθ, ϕ)*sqrt(q(m, reaction_system)), LMs)

# calculate: long ~ 10h
structures = getsamplePWA.(plotdata.x)
constrained_pw_projection_fixed_model!.(structures, 10)

# plot
let
    ps = [
        scatter(map(x->abs2(x.pars[2]), comb), getproperty.(comb, :min)) for
        comb in getproperty.(structures, :samples)
    ]
    plot(ps..., size = (1500, 1000))
end


# pack to a flat tables
summaryline(t::NamedTuple) = [t.min, real.(t.pars)..., imag.(t.pars)...]
function summaryline(s::samplePWA, m)
    ls = summaryline.(s.samples)
    [fill(m, length(ls)) hcat(ls...)']
end

# write to a file
packed_fitresutls = vcat(summaryline.(structures, plotdata.x)...)
writedlm(
    joinpath("data", "exp_pro", "sample_pws_$(settings["tag"]).txt"),
    packed_fitresutls,
)

# read from the file
packed_fitresutls =
    readdlm(joinpath("data", "exp_pro", "sample_pws_$(settings["tag"]).txt"))

# unpacked fitresutls
# Table(packed_fitresutls[:,1:2], (:a,:b))
