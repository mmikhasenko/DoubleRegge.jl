using DrWatson
@quickactivate "DoubleRegge"
using TOML
using TypedTables
using QuadGK
using Cuba
# 
using Plots
import Plots.PlotMeasures.mm
theme(:wong2; size=(500, 350), bottom_margin=5mm)
# 
using LaTeXStrings
# 
using DoubleRegge
# 
using Statistics
using LinearAlgebra

parse_bool_env(name, default) = lowercase(get(ENV, name, string(default))) in ("1", "true", "yes", "on")
parse_int_env(name, default) = parse(Int, get(ENV, name, string(default)))

function select_bin_indices(nbins; stride=1, max_bins=0)
    indices = collect(1:stride:nbins)
    max_bins <= 0 && return indices
    return first(indices, min(length(indices), max_bins))
end

function finite_span(xs)
    length(xs) >= 2 || return nothing
    return xs[[1, end]]
end

function xerr_width(xs)
    length(xs) >= 2 || return nothing
    return (xs[2] - xs[1]) / 2
end

function maybe_vspan!(; sp=nothing, xs)
    span = finite_span(xs)
    span === nothing && return nothing
    kwargs = sp === nothing ? (; lab="", α=0.1, seriescolor=7) : (; sp=sp, lab="", α=0.1, seriescolor=7)
    vspan!(span; kwargs...)
end

function maybe_scatter!(xs, ys; sp=nothing, kwargs...)
    args = sp === nothing ? (; kwargs...) : (; sp=sp, kwargs...)
    width = xerr_width(xs)
    width === nothing ? scatter!(xs, ys; args...) : scatter!(xs, ys; xerr=width, args...)
end

function bootstrap_stdev(f, sampler; nsamples)
    nsamples <= 1 && return fill(0.0, length(f(sampler())))
    stat = hcat([f(sampler()) for _ in 1:nsamples]...)
    return sqrt.(diag(cov(stat; dims=2)))
end

has_table_column(table, name::Symbol) = name in propertynames(table)

function intensity_vectors(table)
    has_table_column(table, :I) && return table.I
    return map(table.Iϕ) do x
        getproperty.(getproperty.(x.PWs, :I), :val)
    end
end

function intensity_error_vectors(table)
    has_table_column(table, :δI) && return table.δI
    return map(table.Iϕ) do x
        getproperty.(getproperty.(x.PWs, :I), :err)
    end
end

function rand_amp(expansion)
    sampled = Vector{NamedTuple{(:I, :ϕ), Tuple{Float64, Float64}}}(collect(
        NamedTuple{(:I, :ϕ)}((max(0.0, rand(pw.I)), rand(pw.ϕ))) for pw in expansion.PWs
    ))
    lms = Vector{Tuple{Int, Int}}(collect(expansion.LMs))
    return changerepresentation(TwoBodyPartialWaves(lms, sampled))
end

rand_amps(table) = rand_amp.(table.Iϕ)

function safe_dNdcosθ(cosθ, expansion)
    list_of_all = collect(zip(expansion.PWs, expansion.LMs))
    l1 = filter(x -> x[2][2] == 1, list_of_all)
    l2 = filter(x -> x[2][2] == 2, list_of_all)
    a1 = sum((a * Psi(L, M, cosθ, π / 2) for (a, (L, M)) in l1); init=0.0 + 0.0im)
    a2 = sum((a * Psi(L, M, cosθ, π / 4) for (a, (L, M)) in l2); init=0.0 + 0.0im)
    return π * (abs2(a1) + abs2(a2))
end


# # # # # # # # # # # # # # # # # # # # 
# 
tag = get(ENV, "DR_TAG", "etappi_a2Po-f2Po-a2f2-f2f2_opposite-sign")
const BIN_STRIDE = parse_int_env("DR_BIN_STRIDE", 2)
const MAX_FIT_BINS = parse_int_env("DR_MAX_FIT_BINS", 6)
const MAX_PLOT_BINS = parse_int_env("DR_MAX_PLOT_BINS", 6)
const PHI_BOOTSTRAP_SAMPLES = parse_int_env("DR_PHI_BOOTSTRAP_SAMPLES", 30)
const COSTHETA_BOOTSTRAP_SAMPLES = parse_int_env("DR_COSTHETA_BOOTSTRAP_SAMPLES", 100)
const COSTHETA_POINTS = parse_int_env("DR_COSTHETA_POINTS", 61)
const SKIP_CONTRIBUTIONS = parse_bool_env("DR_SKIP_CONTRIBUTIONS", false)
const MERGE_PDFS = parse_bool_env("DR_MERGE_PDFS", false)
# 
# # # # # # # # # # # # # # # # # # # # 


# 
settings_file = fitsfolder(tag, "fit-results.toml")
!isfile(settings_file) && error("no file")
parsed = TOML.parsefile(settings_file)
@unpack settings, fit_results = parsed

const reaction_system = getproperty(DoubleRegge, Symbol(settings["system"]))


# build model
const model = DoubleReggeModel(
    sixexchages[settings["exchanges"]],
    settings["t2"],
    settings["scale_α"],
    reaction_system,
    fit_results["fit_minimizer"];
    s2shift=get(settings, "s2_shift", 0.0))
const pfr = model.pars
fixed_model(m, cosθ, ϕ; pars=pfr) = amplitude(with_parameters(model, pars), m, cosθ, ϕ)
fixed_model_sqrtq(m, cosθ, ϕ; pars=pfr) = fixed_model(m, cosθ, ϕ; pars=pars) * sqrt(q(m, reaction_system))
intensity(m, cosθ, ϕ) = abs2(fixed_model(m, cosθ, ϕ)) * q(m, reaction_system)

# get data
data = read_data(settings["pathtodata"], reaction_system)
# fit range
fitdata = filter(data) do x
    inlims(x.x, settings["fitrange"])
end
fitdata = fitdata[select_bin_indices(length(fitdata); stride=BIN_STRIDE, max_bins=MAX_FIT_BINS)]
# plot 
plotdata = filter(data) do x
    inlims(x.x, (2.4, 3.0))
end
plotdata = plotdata[select_bin_indices(length(plotdata); stride=BIN_STRIDE, max_bins=MAX_PLOT_BINS)]

function pw_project_fixed(m::Float64, L, M)
    amplitude(cosθ, ϕ) = fixed_model_sqrtq(m, cosθ, ϕ)
    return pw_project(amplitude, L, M)
end
# 
model_integral(m; pars=pfr) = integrate_dcosθdϕ((cosθ, ϕ) -> abs2(fixed_model_sqrtq(m, cosθ, ϕ; pars=pars)))[1]
model_integral_forward(m; pars=pfr) = integrate_dcosθdϕ((cosθ, ϕ) -> abs2(fixed_model_sqrtq(m, cosθ, ϕ; pars=pars)), (0, 1))[1]
model_integral_backward(m; pars=pfr) = integrate_dcosθdϕ((cosθ, ϕ) -> abs2(fixed_model_sqrtq(m, cosθ, ϕ; pars=pars)), (-1, 0))[1]

const Np = length(pfr)
δ(i; n=Np) = (1:n .== i)
integral_interf(m, i, j) = model_integral(m; pars=pfr .* δ(i) + pfr .* δ(j)) -
                           model_integral(m; pars=pfr .* δ(i)) -
                           model_integral(m; pars=pfr .* δ(j));
# 
contribution_matrix(m) = [(i > j ? 0 : integral_interf(m, i, j) / (i == j ? 2 : 1)) for i in 1:Np, j in 1:Np]
if !SKIP_CONTRIBUTIONS
    let
        m = sum(contribution_matrix.(fitdata.x))
        heatmap([mi == 0.0 ? NaN : mi for mi in m],
            xticks=(1:Np, getproperty.(model.exchanges, :label)),
            yticks=(1:Np, getproperty.(model.exchanges, :label)), colorbar=false)
        # 
        mn = m ./ sum(m)
        for i in CartesianIndices(m)
            (i[1] > i[2]) && continue
            annotate!([(i[2], i[1],
                text("$(Int(round(m[i], digits=0))) / $(Int(round(100*mn[i], digits=0)))%", 10, :red))])
        end
        plot!(size=(400, 350), title="contributions of different diagrams")
        ellh = fit_results["fit_minimum"]
        s = prod("$l: $v,\n" for (l, v) in zip(getproperty.(model.exchanges, :label), round.(pfr, digits=2)))
        s *= "extLLH: " * string(round(ellh, digits=1))
        annotate!([(0.6, 3, text(s, 10, :left))])
    end
    savefig(
        joinpath("data", "exp_pro", tag,
            "contributions_$(tag)_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf"))
end
# 

# data
const LMs = reaction_system.LMs
data = read_data(settings["pathtodata"], reaction_system)
# fit range
fitrangemap = map(x -> inlims(x.x, settings["fitrange"]), data)
fitdata = data[fitrangemap]
fitdata = fitdata[select_bin_indices(length(fitdata); stride=BIN_STRIDE, max_bins=MAX_FIT_BINS)]
# plot 
plotmap = map(x -> inlims(x.x, (2.2, 3.0)), data)
plotdata = data[plotmap]
plotdata = plotdata[select_bin_indices(length(plotdata); stride=BIN_STRIDE, max_bins=MAX_PLOT_BINS)]
fit_I = intensity_vectors(fitdata)
plot_I = intensity_vectors(plotdata)
plot_δI = intensity_error_vectors(plotdata)


# ellh
fit_results["fit_minimum"]

# constraint
sum(sum, fit_I), sum(model_integral.(fitdata.x))

# intensity
intensity_in_bins = model_integral.(plotdata.x)

# asymmetry
intensity_forward_in_bins = model_integral_forward.(plotdata.x)
intensity_backward_in_bins = model_integral_backward.(plotdata.x)
data_forward_in_bins = [quadgk(cosθ -> safe_dNdcosθ(cosθ, a), 0, 1)[1] for a in plotdata.amps]
data_backward_in_bins = [quadgk(cosθ -> safe_dNdcosθ(cosθ, a), -1, 0)[1] for a in plotdata.amps]
# 
asymmetry(f, b) = (f - b) / (f + b)
δasymmetry(f, b, δ) = 2δ * sqrt(f^2 + b^2) / (f + b)^2

asymm_data = [asymmetry(f, b) for (f, b) in zip(data_forward_in_bins, data_backward_in_bins)]
δasymm_data = [δasymmetry(f, b, δ) for (f, b, δ) in zip(
    data_forward_in_bins,
    data_backward_in_bins,
    sqrt.(sum.(x -> x^2, plot_δI) / 2))]
#
asymm_model = [asymmetry(f, b) for (f, b) in zip(intensity_forward_in_bins, intensity_backward_in_bins)]
# 
let
    plot(size=(900, 350), layout=grid(1, 2),
        xlab="m(ηπ) (GeV)",
        ylab=["intensity" "(F-B) / (F+B)" "intensity" "intensity"],
        title=["number of events" "asymmetry"])
    #
    common_options = (m=(3,),)
    # 
    maybe_scatter!(plotdata.x, data_forward_in_bins; sp=1, lab="forward",
        common_options..., yerr=sqrt.(sum.(x -> x^2, plot_δI) / 2), seriescolor=3)
    maybe_scatter!(plotdata.x, data_backward_in_bins; sp=1, lab="backward",
        common_options..., yerr=sqrt.(sum.(x -> x^2, plot_δI) / 2), seriescolor=4)
    plot!(sp=1, plotdata.x, intensity_forward_in_bins, lab="", lw=2, seriescolor=3)
    plot!(sp=1, plotdata.x, intensity_backward_in_bins, lab="", lw=2, seriescolor=4)
    # 
    maybe_scatter!(plotdata.x, asymm_data; sp=2, lab="", common_options..., yerr=δasymm_data)
    plot!(sp=2, plotdata.x, asymm_model, lab="", ylim=(-1, 1), lw=2)
    #
    maybe_scatter!(plotdata.x, sum.(plot_I); sp=1,
        common_options..., lab="intensity",
        yerr=sqrt.(sum.(x -> x^2, plot_δI)), seriescolor=2)
    plot!(sp=1, plotdata.x, intensity_in_bins, lw=2, lab="", seriescolor=2)
    # 
    maybe_vspan!(sp=1; xs=fitdata.x)
    maybe_vspan!(sp=2; xs=fitdata.x)
end
savefig(
    joinpath("data", "exp_pro", tag,
        "intensity-assymetry_$(tag)_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf"))
#
# let bin = 1
#     mηπ = fitdata.x[bin]
#     cosθv = range(-1,1, length=100)
#     ϕv = range(-π, π, length=95)
#     calv = intensity.(mηπ,cosθv',ϕv)
#     heatmap(cosθv, ϕv, calv, colorbar=false,
#         xlab=L"\cos\theta", ylab=L"\phi", title="$(round(mηπ, digits=2)) GeV")
# end

setfirstargument(x1, f) = (x2, x3) -> f(x1, x2, x3)
setthirdfourtharguments(x3, x4, f) = (x1, x2) -> f(x1, x2, x3, x4)
# 
function phi_moment(f, l; forward=true)
    cosrange = forward ? (0, 1) : (-1, 0)
    return integrate_dcosθdϕ((cosθ, ϕ) -> cos(l * ϕ) * f(cosθ, ϕ), cosrange)[1] /
           integrate_dcosθdϕ(f, cosrange)[1]
end

phi_moment_data(amps; l, forward) = [
    phi_moment((cosθ, ϕ) -> abs2(recamp(cosθ, ϕ, a)), l; forward=forward) for a in amps
]
#

function take_stdiv_of_vector(vector_of_vector)
    mat = hcat(vector_of_vector...)
    err = sqrt.(diag(cov(mat; dims=2)))
    return err
end

# phi asymmetry
function phiasymmplot(l)
    plot(xlab=L"m_{\eta\pi}\,(\textrm{GeV})", ylab=(l == 1 ? L"<\cos\,\phi>" : L"<\cos\,2\phi>"), size=(500, 350))
    plot!(plotdata.x, phi_moment.(setfirstargument.(plotdata.x, intensity), l; forward=true), lw=2, c=2, lab="forward")
    plot!(plotdata.x, phi_moment.(setfirstargument.(plotdata.x, intensity), l; forward=false), lw=2, c=3, lab="backward")
    # let
    vals = phi_moment_data(plotdata.amps; l=l, forward=true)
    err = bootstrap_stdev(amps -> phi_moment_data(amps; l=l, forward=true), () -> rand_amps(plotdata); nsamples=PHI_BOOTSTRAP_SAMPLES)
    maybe_scatter!(plotdata.x, vals; yerr=err, c=2, lab="")
    #
    vals = phi_moment_data(plotdata.amps; l=l, forward=false)
    err = bootstrap_stdev(amps -> phi_moment_data(amps; l=l, forward=false), () -> rand_amps(plotdata); nsamples=PHI_BOOTSTRAP_SAMPLES)
    maybe_scatter!(plotdata.x, vals; yerr=err, c=3, lab="")
end
plot(phiasymmplot(1), phiasymmplot(2), size=(900, 350), layout=grid(1, 2))
maybe_vspan!(sp=1; xs=fitdata.x)
maybe_vspan!(sp=2; xs=fitdata.x)
savefig(
    joinpath("data", "exp_pro", tag,
        "intensity-cosphi_$(tag)_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf"))

# cosθ distributions
let
    cosθv = range(-1, 1, length=COSTHETA_POINTS)
    function make_plot(bin)
        calv = safe_dNdcosθ.(cosθv, Ref(fitdata.amps[bin]))
        #
        err = bootstrap_stdev(amps -> safe_dNdcosθ.(cosθv, Ref(amps[bin])), () -> rand_amps(fitdata); nsamples=COSTHETA_BOOTSTRAP_SAMPLES)
        #
        plot(xlab="cosθ", title="$(round(fitdata.x[bin]; digits=2)) GeV")
        plot!(cosθv, calv, ribbon=err, lab="")
        # 
        projection(cosθ) = quadgk(ϕ -> intensity(fitdata.x[bin], cosθ, ϕ), -π, π)[1]
        plot!(cosθv, projection.(cosθv), lab="", lw=2, yaxis=nothing)
    end
    ps = make_plot.(1:length(fitdata.x))
    plot(ps..., size=(1100, 500), layout=grid(3, 5))
end
savefig(
    joinpath("data", "exp_pro", tag,
        "cos-distributions_$(tag)_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf"))
#


# # projections
pw_projections = [map(LM -> pw_project_fixed(m, LM...), LMs) for m in plotdata.x]
pw_intensities = map(x -> abs2.(x), pw_projections)
let
    plot(layout=grid(3, 3), size=(900, 900))
    for (i, (L, M)) in enumerate(LMs)
        maybe_scatter!(plotdata.x, getindex.(plot_I, i); sp=i,
            yerr=getindex.(plot_δI, i),
            c=:black, title="LM=$L$M", ms=3,
            lab=i != 1 ? "" : "data",)
        #
        plot!(sp=i, plotdata.x, getindex.(pw_intensities, i), lab=i != 1 ? "" : "PW projection", l=(2,))
        maybe_vspan!(sp=i; xs=fitdata.x)
    end
    plot!(xlab="m(ηπ) (GeV)")
end
savefig(
    joinpath("data", "exp_pro", tag,
        "pw-projections_$(tag)_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf"))
#
tolab(LM) = "$(LM[1])$(LM[2])"

bar([sum(getindex.(fit_I, i)) for (i, LM) in enumerate(LMs)],
    xticks=(1:7, tolab.(LMs)), xlab="LM", ylab="intensity", lab="")

bar([sum(getindex.(pw_intensities[6:end], i)) for (i, LM) in enumerate(LMs)],
    xticks=(1:7, tolab.(LMs)), xlab="LM", ylab="intensity", lab="")
#

# Odd and Even waves
fHeigher = (intensity_in_bins .- sum.(pw_intensities)) ./ intensity_in_bins
filtodd = isodd.(getindex.(LMs, 1))
fOdd = map(x -> sum(x .* filtodd), pw_intensities) ./ intensity_in_bins
fEven = map(x -> sum(x .* iszero.(filtodd)), pw_intensities) ./ intensity_in_bins
# 
fOdd_compass = map(x -> sum(x .* filtodd), plot_I) ./ sum.(plot_I)
fEven_compass = map(x -> sum(x .* iszero.(filtodd)), plot_I) ./ sum.(plot_I)

let
    plot(ylab="fraction", xlab="m(ηπ) (GeV)", size=(500, 350), title=tag)
    plot!(plotdata.x, fHeigher, lab="Higher waves L > 6", lw=2)
    plot!(plotdata.x, fEven, lab="Even waves L ≤ 6", lw=2)
    plot!(plotdata.x, fOdd, lab="Odd waves L ≤ 6", lw=2)
    maybe_scatter!(plotdata.x, fEven_compass; lab="", mc=2, ms=3)
    maybe_scatter!(plotdata.x, fOdd_compass; lab="", mc=3, ms=3)
    plot!(ylims=(0, 1), leg=:left)
    maybe_vspan!(xs=fitdata.x)
end
savefig(
    joinpath("data", "exp_pro", tag,
        "odd-and-even_$(tag)_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf"))
#


let
    pathtofolder = joinpath("data", "exp_pro", tag)
    inputfiles = readdir(pathtofolder)
    outputfile = joinpath(pathtofolder, "_combined.pdf")
    inputfiles = filter(f -> splitext(f)[2] == ".pdf" && f != outputfile, inputfiles)
end


produced_files = filter(!isnothing, [
    "intensity-assymetry_$(tag)_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf",
    "cos-distributions_$(tag)_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf",
    SKIP_CONTRIBUTIONS ? nothing : "contributions_$(tag)_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf",
    "odd-and-even_$(tag)_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf",
    "intensity-cosphi_$(tag)_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf",
    "pw-projections_$(tag)_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf"
])

if MERGE_PDFS
    let
        pathtofolder = joinpath("data", "exp_pro", tag)
        outputfile = joinpath(pathtofolder, "combined_$(tag)_Np=$(length(settings["exchanges"]))_alpha=$(settings["scale_α"]).pdf")
        inputfiles = joinpath.(Ref(pathtofolder), produced_files)
        run(`pdftk $inputfiles cat output $outputfile`)
    end
end
