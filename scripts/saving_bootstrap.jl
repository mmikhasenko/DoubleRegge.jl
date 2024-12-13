using DrWatson
@quickactivate "DoubleRegge"

using DoubleRegge
using QuadGK
using Plots
using DelimitedFiles
# 
theme(:wong2;
    frame = :box, grid = false, lab = "",
    xlims = (:auto, :auto), ylims = (:auto, :auto))
#

data = read_data(settings["pathtodata"], description_ηπ);
const xdata = data.x;
const Nbins = length(xdata)
#
const amplitudes = data.amps;

let nPoints = 100
    main_point = Matrix{Float64}(undef, Nbins, nPoints)
    cosθv = range(-1, 1, nPoints)
    for bin in 1:Nbins
        calv = dNdcosθ.(cosθv, Ref(amplitudes[bin]))
        main_point[bin, :] .= calv
    end
    writedlm(datadir("exp_pro", "main_point.txt"), main_point)
end

let nPoints = 100, nSamples = 1000
    lower_bands = Matrix{Float64}(undef, Nbins, nPoints)
    upper_bands = Matrix{Float64}(undef, Nbins, nPoints)
    for bin in 1:Nbins
        v1, v2 = bootstrap_band(data.Iϕ[bin]; nPoints, nSamples)
        lower_bands[bin, :] .= v1
        upper_bands[bin, :] .= v2
    end
    writedlm(datadir("exp_pro", "bootstrap_lower.txt"), lower_bands)
    writedlm(datadir("exp_pro", "bootstrap_upper.txt"), upper_bands)
end

let
    lower_bands = readdlm(datadir("exp_pro", "bootstrap_lower.txt"))
    upper_bands = readdlm(datadir("exp_pro", "bootstrap_upper.txt"))
    Np = size(lower_bands, 2)
    ps = map(1:Nbins) do bin
        plot(range(-1, 1, length = Np), lower_bands[bin, :],
            fill_between = upper_bands[bin, :], lab = "", xaxis = nothing, yaxis = nothing)
    end
    plot(ps..., size = (1000, 1000))
end
savefig(plotsdir("bootstrap_interval_costheta.pdf"))
