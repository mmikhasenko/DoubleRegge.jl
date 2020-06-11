using DrWatson
@quickactivate "DoubleRegge"

using DoubleRegge
using QuadGK
using Plots
using DelimitedFiles
theme(:wong; size=(500,350))
#
data = [(inten = get_intesity(L,M),
         phase = get_phase(L,M)) for (L,M) in LMs];
const xdata = data[1].inten[:,1];
const Nbins = length(xdata)
#
const amplitudes = constructamps(
    slice(data, :inten, 1:Nbins),
    slice(data, :phase, 1:Nbins)
);

@time bootstrap_band(10, data)

let Npoints = 100, Nsamples=1000
    lower_bands = Matrix{Float64}(undef,Nbins,Npoints)
    upper_bands = Matrix{Float64}(undef,Nbins,Npoints)
    for bin in 1:Nbins
        v1, v2 = bootstrap_band(bin, data; Npoints=Npoints, Nsamples=Nsamples)
        lower_bands[bin,:] .= v1
        upper_bands[bin,:] .= v2
    end
    writedlm(datadir("exp_pro","bootstrap_lower.txt"), lower_bands)
    writedlm(datadir("exp_pro","bootstrap_upper.txt"), upper_bands)
end

let
    lower_bands = readdlm(datadir("exp_pro","bootstrap_lower.txt"))
    upper_bands = readdlm(datadir("exp_pro","bootstrap_upper.txt"))
    Np = size(lower_bands,2)
    ps = [plot(range(-1,1,length=Np), lower_bands[bin,:], fill_between=upper_bands[bin,:], lab="") for bin in 1:Nbins]
    plot(ps..., size=(1000,1000))
end
