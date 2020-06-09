using DrWatson
@quickactivate "DoubleRegge"

using DoubleRegge
using QuadGK
using Plots
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

# test
let bin=40
    i1 = quadgk(cosθ->dNdcosθ(cosθ; amps=slice(amplitudes, bin)),-1, 1)[1]
    i2 = sum(slice(data, :inten, bin))
    i1 ≈ i2
end

using Statistics

let bin=40
    i = slice(data, :inten, bin)
    δi = slice(data, :inten, bin; ind=3)
    bar(i, yerr = δi, linecolor=:black, markerstrokecolor=:black,
        xticks=(1:length(i), ["($L,$M)" for (L,M) in LMs]), xlab="waves", lab="")
end
savefig(plotsdir("intensities.with.errors.2.32.pdf"))

let bin=40
    i = slice(data, :inten, bin)
    # errors
    δi = slice(data, :inten, bin; ind=3)
    # computation
    function randsumi()
        i .+ randn(length(i)) .* δi
        i′ = i .+ randn(length(i)) .* δi
        i′ = map(v -> v < 0 ? 0.0 : v, i′)
        # i′ .*= sum(i)/sum(i′)
        return sum(i′)
    end
    stephist([randsumi() for _ in 1:1000], xlab="sum(I)", lab="bootstrap")
    vline!([sum(i)], l=(3,:black), lab="N")
    vline!(sum(i) .+ [-1,1] .* sqrt(sum(i)), lab="±√N", l=(2,:red))
end
savefig(plotsdir("sumI_sqrtN.pdf"))

let bin=40, bin1 = 10, bin2=90
    i = slice(data, :inten, bin)
    p = slice(data, :phase, bin)
    # errors
    δi = slice(data, :inten, bin; ind=3)
    δp = slice(data, :phase, bin; ind=3)
    # computation
    amps = constructamps(i,p)
    cosθv = range(-1,1, length=101)
    calv = dNdcosθ.(cosθv; amps=amps)
    # bootstrap
    bstp = [let
                i′ = i .+ randn(length(i)) .* δi
                i′ = map(v -> v < 0 ? 0.0 : v, i′)
                i′ .*= sum(i) / sum(i′)
                p′ = p+randn()*δp
                dNdcosθ.(cosθv; amps=constructamps(i′,p′))
        end for _ in 1:1000]
    plot(layout=grid(2,2), size=(900,700), title=["examples" "ribbon" "slice1" "slice2"])
    [plot!(sp=1, cosθv, b, lab="") for b in bstp[1:20]]
    #
    qls = [quantile(getindex.(bstp,i), [0.16, 0.84]) for i in 1:length(bstp[1])]
    plot!(sp=2, cosθv, calv, ribbon=(calv-getindex.(qls,1), getindex.(qls,2)-calv), lab="")
    vline!(sp=2, [cosθv[bin1], cosθv[bin2]], l=(:gray,:dash), lab="")
    #
    #
    s = getindex.(bstp,bin1)
    stephist!(sp=3, s, bins=50, lab="")
    vline!(sp=3, quantile(s,[0.16, 0.84]), lab="", l=(1,:black))
    vline!(sp=3, [calv[bin1]], lab="", l=(2,:red))
    # #
    s = getindex.(bstp,bin2)
    stephist!(sp=4, s, bins=50, lab="")
    vline!(sp=4, quantile(s,[0.16, 0.84]), lab="", l=(1,:black))
    vline!(sp=4, [calv[bin2]], lab="", l=(2,:red))
end
savefig(plotsdir("bootstrap.2.32.pdf"))
# savefig(plotsdir("bootstrap.2.32.norm.pdf"))
