using DrWatson
@quickactivate "DoubleRegge"

using Plots
theme(:wong; size=(500,350))

using QuadGK
using Cuba
using ForwardDiff
# 
using Optim
using TypedTables
#
using TOML

# 
using DoubleRegge

# current cesar model
settings = Dict(
    "system" => "compass_ηπ",
    "pathtodata" => joinpath("data","exp_raw","PLB_shifted"),
    "fitrange" => [2.4, 3.0],
    "t2" => -0.2,
    "tag" => "bottom-Po",
    "exchanges" => [1,3,5],
    "initial_pars" => [0.7, -0.7, 0.0 ],
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

# a2/P = 0.55770;  f2/P=0.72755 ; P/P=-0.00225
# ellh = -1245.33276

getindex.(exchanges,4)

function integrand(cosθ,ϕ,pars)
    Id = abs2.(recamp.(cosθ, ϕ, fitdata.amps, Ref(LMs)))
    Am = model.(fitdata.x, cosθ, ϕ; pars=pars)
    Im = abs2.(Am) .* q.(fitdata.x)
    return sum(Im .- Id .* log.(Im))
end

ellh(pars) = integrate_dcosθdϕ((cosθ,ϕ)->integrand(cosθ,ϕ,pars))[1]
#
# @time 
integrand′(cosθ,ϕ,pars) = ForwardDiff.gradient(p->integrand(cosθ,ϕ,p), pars)
ellh′(pars) = integrate_dcosθdϕ((cosθ,ϕ)->integrand′(cosθ,ϕ,pars), dims=3)
ellh′!(stor,pars) = (stor .= ellh′(pars))

# ellh(ft.minimizer)
# ellh′(ft.minimizer)
# ellh([0.55770, 0.72755, -0.00225])
# ellh([0.55659, -0.71550, 0.00384])
collect(zip(getindex.(exchanges,4),ft_deriv.minimizer))

# @show ft_deriv
ellh([0.55770, 0.72755,-0.00225])
ft_deriv.minimum
ellh(ft_deriv.minimizer)

ft_deriv = Optim.optimize(ellh, ellh′!, [0.6085504973680421, 0.7751736383298719, -0.004648902500409368], BFGS(),
               Optim.Options(show_trace = true))

ft = Optim.optimize(ellh, [0.6085504973680421, 0.7751736383298719, -0.004648902500409368], BFGS(),
               Optim.Options(show_trace = true))
#
fit_results = Dict(
    "fit_converged" => ft.x_converged,
    "fit_minimizer" => ft.minimizer,
    "fit_minimum" => ft.minimum)

output_name = joinpath("data", "exp_pro","fit-results_$(settings["tag"])_Np=$(length(settings["exchanges"]))");

open(output_name,"w") do io
    TOML.print(io, Dict(
            "settings"=>settings,
            "fit_results"=>fit_results))
end
