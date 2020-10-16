using DrWatson
@quickactivate "DoubleRegge"

using Test

using QuadGK
using Cuba

using Plots
theme(:wong; size=(500,350))

using DoubleRegge
setsystem!(:compass_ηπ)

using TOML

# settings_file = joinpath("data", "exp_pro","fit-results_bottom-Po_Np=3.toml")
settings_file = joinpath("data", "exp_pro","fit-results_a2Po-f2f2-PoPo_Np=3.toml")
! isfile(settings_file) && error("no file")
# 
parsed = TOML.parsefile(settings_file)
settings = parsed["settings"]
fit_results = parsed["fit_results"]

# 
setsystem!(Symbol(settings["system"]))


# fit
exchanges = sixexchages[settings["exchanges"]]
model = build_model(exchanges, G.s0, settings["t2"], settings["scale_α"])
fixed_model(m,cosθ,ϕ) = model(m,cosθ,ϕ; pars=fit_results["fit_minimizer"])


# data
LMs = compass_ηπ_LMs
data = Table(x_IδI_ϕδϕ_compass_ηπ(settings["pathtodata"]))
amplitudes = [sqrt.(is) .* cis.(ϕs) for (is,ϕs) in zip(data.I, data.ϕ)]
# fit range
fitrangemap = map(x->inlims(x.x, settings["fitrange"]), data)
fitdata = Table(data[fitrangemap], amps = amplitudes[fitrangemap])
# plot 
plotmap = map(x->inlims(x.x, (2.4,3.5)), data)
plotdata = data[plotmap]


function pw_iϕ(m,L,M)
    amplitude(cosθ,ϕ) = fixed_model(m,cosθ,ϕ)*sqrt(q(m))
    return pw_project(amplitude,L,M)
end
#
model_integral(m) = integrate_dcosθdϕ((cosθ,ϕ)->abs2(fixed_model(m,cosθ,ϕ))*q(m))
#
sum(sum, fitdata.I)
sum(model_integral.(fitdata.x))

# 
calv = [pw_iϕ.(data.x[plotmap],L,M) for (L,M) in LMs]

fit_results["fit_minimum"]

let
    plot(layout=grid(3,3), size=(900,900))
    for (i,(m,(L,M))) in enumerate(zip(calv,LMs))
        scatter!(sp=i, plotdata.x, getindex.(plotdata.I, i),
            yerr=getindex.(plotdata.δI, i), xerr=(plotdata.x[2]-plotdata.x[1])/2,
            c=:black, title="LM=$L$M", ms=3,
            lab=i!=1 ? "" : "data",)
        plot!(sp=i, plotdata.x, getindex.(m,1), lab=i!=1 ? "" : "a2Po-f2Po-PoPo", l=(2))
        vspan!(sp=i, fitdata.x[[1,end]], lab="", α=0.2, seriescolor=1)
    end
    plot!(xlab="m(ηπ) (GeV)")
end

savefig(
    joinpath("data", "exp_pro",
        "pw-projections_$(settings["tag"])_Np=$(length(settings["exchanges"])).pdf"))
#
