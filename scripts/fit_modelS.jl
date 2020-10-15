using DrWatson
@quickactivate "DoubleRegge"

using DoubleRegge
using Plots
theme(:wong; size=(500,350))

using QuadGK
using Cuba
# 
using Optim
using TypedTables
#
#


inlims(x,lims) = lims[1] ≤ x ≤ lims[2]


settings = Dict(
    "system" => :compass_ηπ,
    "pathtodata" => joinpath("data","exp_raw","PLB_shifted"),
    "system" => :compass_ηπ,
    "fitrange" => (2.5, 3.0),
    "t2" => -0.2,
    "tag" => "bottom-Po",
    "exchanges" => [1,3,5],
    "initial_pars" => [0.7, -0.7, 0.0 ],
    "scale_α" => 0.9,
)

setsystem!(settings["system"])

# data


LMs = compass_ηπ_LMs
data = Table(x_IδI_ϕδϕ_compass_ηπ(settings["pathtodata"]))
amplitudes = [sqrt.(is) .* cis.(ϕs) for (is,ϕs) in zip(data.I, data.ϕ)]

# range
fitrangemap = map(x->inlims(x.x, settings["fitrange"]), data)
fitdata = Table(data[fitrangemap], amps = amplitudes[fitrangemap])

# fit
exchanges = sixexchages[settings["exchanges"]]
model = build_model(exchanges, G.s0, settings["t2"], settings["scale_α"])
intensity(m, cosθ, ϕ; pars) = abs2(model(m, cosθ, ϕ; pars=pars))*q(m)

function integrand(cosθ,ϕ,pars)
    Id = abs2.(recamp.(cosθ, ϕ, fitdata.amps, Ref(LMs)))
    Am = model.(fitdata.x, cosθ, ϕ; pars=pars)
    Im = abs2.(Am) .* q.(fitdata.x)
    return sum(Im .- Id .* log.(Im))
end

ellh(pars) = integrate_dcosθdϕ(x->integrand(2x[1]-1,0.3+π*(2x[2]-1),pars))
#
@time ellh([1.0,0,0])
# [ 0.62, -0.8, +0.0057]
ft = Optim.optimize(ellh, settings["initial_pars"], BFGS(),
               Optim.Options(show_trace = true))
#
ft.minimizer

ft.status

# 0.8: [0.5658807102768846, -0.7473538045573105, 0.00403669556892034]

let
    cosθv = range(-1,1, length=101)
    function make_plot(bin)
        calv = dNdcosθ.(cosθv; amps=fitdata.amps[bin], LMs=LMs)
        plot(cosθv, calv, lab="")
        projection(cosθ) = quadgk(ϕ->intensity(fitdata.x[bin], cosθ, ϕ; pars=ft.minimizer), -π, π)[1]
        plot!(cosθv, projection.(cosθv), lab="")
    end
    ps = make_plot.(1:length(fitdata.x))
    plot(ps..., size=(1000,600))
end
