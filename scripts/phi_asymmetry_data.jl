using DrWatson
@quickactivate "DoubleRegge"

using DoubleRegge
using QuadGK
using Plots
using Statistics
theme(:wong; size=(500,350))
#
data = [(inten = get_intesity(L,M),
         phase = get_phase(L,M)) for (L,M) in LMs];
const xdata = data[1].inten[:,1];
const Nbins = length(xdata)
#
slice(data, property, bin; ind=2) = [v[bin, ind] for v in getproperty.(data, property)]
slice(data, bin) = getindex.(data, bin)
#
amplitudes = constructamps(
    slice(data, :inten, 1:Nbins),
    slice(data, :phase, 1:Nbins)
);

let bin=40
    cosθv = range(-1,1, length=101)
    calv = dNdcosθ.(cosθv; amps=slice(amplitudes, bin))
    plot(cosθv, calv)
end

phi_asymmetry_data(bin, cosθ) = 
    phi_asymmetry(ϕ->abs2(recamp(cosθ,ϕ,slice(amplitudes,bin))))
#

let 
    plot()
    for bin in 40:3:55
        plot!(cosθ->phi_asymmetry_data(bin, cosθ), -1:0.01:1, lab="m(ηπ)=$(round(xdata[bin], digits=2))")
    end
    plot!(xlab="cosθ", ylab="Assymetry in position of ϕ peak",)
end