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

const data_folder = "data/exp_raw/PLB_shifted"
const data = read_data(data_folder, description_ηπ);
const nBins = length(data)

let nPoints = 100
    main_point = Matrix{Float64}(undef, nBins, nPoints)
    cosθv = range(-1, 1, nPoints)
    for bin in 1:nBins
        calv = dNdcosθ.(cosθv, Ref(data.amps[bin]))
        main_point[bin, :] .= calv
    end
    writedlm(datadir("exp_pro", "main_point.txt"), main_point)
end

let nPoints = 100, nSamples = 1000
    lower_bands = Matrix{Float64}(undef, nBins, nPoints)
    upper_bands = Matrix{Float64}(undef, nBins, nPoints)
    for bin in 1:nBins
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
    ps = map(1:nBins) do bin
        plot(range(-1, 1, length = Np), lower_bands[bin, :],
            fill_between = upper_bands[bin, :], lab = "", xaxis = nothing, yaxis = nothing)
    end
    plot(ps..., size = (1000, 1000))
end
savefig(plotsdir("bootstrap_interval_costheta.pdf"))
