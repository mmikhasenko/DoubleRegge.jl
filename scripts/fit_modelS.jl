using DrWatson
@quickactivate "DoubleRegge"

using Plots
theme(:wong; size=(500,350))

using QuadGK
using Cuba
# 
using Optim
using TypedTables
using ForwardDiff
#
using TOML

# 
using DoubleRegge

# bottom-Po model same sign
# settings = Dict(
#     "system" => "compass_ηπ",
#     "pathtodata" => joinpath("data","exp_raw","PLB_shifted"),
#     "fitrange" => [2.4, 3.0],
#     "t2" => -0.2,
#     "tag" => "a2Po-f2Po-PoPo_oppoite-sign",
#     "exchanges" => [1,3,5],
#     "initial_pars" => [0.7, 0.7, 0.0 ],
#     "scale_α" => 0.8,
# )

# # bottom-Po model opposite sign
# settings = Dict(
#     "system" => "compass_ηπ",
#     "pathtodata" => joinpath("data","exp_raw","PLB_shifted"),
#     "fitrange" => [2.4, 3.0],
#     "t2" => -0.2,
#     "tag" => "a2Po-f2Po-PoPo_opposite-sign",
#     "exchanges" => [1,3,5],
#     "initial_pars" => [0.7, -0.7, 0.0 ],
#     "scale_α" => 0.8,
# )

# cesar model with f2/f2: same sign
# settings = Dict(
#     "system" => "compass_ηπ",
#     "pathtodata" => joinpath("data","exp_raw","PLB_shifted"),
#     "fitrange" => [2.4, 3.0],
#     "t2" => -0.2,
#     "tag" => "a2Po-f2f2-PoPo_same-sign",
#     "exchanges" => [1,4,5],
#     "initial_pars" => [0.7, 11.0, 0.0],
#     "scale_α" => 0.8,
# )

# cesar model with f2/f2: opposite sign
settings = Dict(
    "system" => "compass_ηπ",
    "pathtodata" => joinpath("data","exp_raw","PLB_shifted"),
    "fitrange" => [2.4, 3.0],
    "t2" => -0.2,
    "tag" => "a2Po-f2f2-PoPo_opposite-sign",
    "exchanges" => [1,4,5],
    "initial_pars" => [0.7, -11.0, 0.0],
    "scale_α" => 0.8,
)


setsystem!(Symbol(settings["system"]))

# data
const LMs = compass_ηπ_LMs
data = Table(x_IδI_ϕδϕ_compass_ηπ(settings["pathtodata"]))
amplitudes = [sqrt.(is) .* cis.(ϕs) for (is,ϕs) in zip(data.I, data.ϕ)]

# range
fitrangemap = map(x->inlims(x.x, settings["fitrange"]), data)
fitdata = Table(data[fitrangemap], amps = amplitudes[fitrangemap])

# fit
const exchanges = sixexchages[settings["exchanges"]]
const model = build_model(exchanges,settings["t2"], settings["scale_α"])
intensity(m, cosθ, ϕ; pars) = abs2(model(m, cosθ, ϕ; pars=pars))*q(m)

function integrand(cosθ,ϕ,pars)
    Id = abs2.(recamp.(cosθ, ϕ, fitdata.amps, Ref(LMs)))
    Am = model.(fitdata.x, cosθ, ϕ; pars=pars)
    Im = abs2.(Am) .* q.(fitdata.x)
    return sum(Im .- Id .* log.(Im))
end
ellh(pars) = integrate_dcosθdϕ((cosθ,ϕ)->integrand(cosθ,ϕ,pars))[1]

# 
integrand′(cosθ,ϕ,pars) = ForwardDiff.gradient(p->integrand(cosθ,ϕ,p), pars)
ellh′(pars) = integrate_dcosθdϕ((cosθ,ϕ)->integrand′(cosθ,ϕ,pars), dims=3)
ellh′!(stor,pars) = (stor .= ellh′(pars)) 

#
@time ellh([1.0,0,0])
# [ 0.62, -0.8, +0.0057]
ft = Optim.optimize(ellh, ellh′!, settings["initial_pars"], BFGS(),
               Optim.Options(show_trace = true, g_tol=1e-4, iterations=15))
#
fit_results = Dict(
    "fit_converged" => ft.x_converged,
    "fit_minimizer" => ft.minimizer,
    "fit_minimum" => ft.minimum)

output_name = joinpath("data", "exp_pro", settings["tag"],
    "fit-results_$(settings["tag"])_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).toml");

open(output_name,"w") do io
    TOML.print(io, Dict(
            "settings"=>settings,
            "fit_results"=>fit_results))
end
