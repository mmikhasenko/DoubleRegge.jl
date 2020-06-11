#
slice(data, property, bin; ind=2) = [v[bin, ind] for v in getproperty.(data, property)]
slice(data, bin) = getindex.(data, bin)
#

function bootstrap_band(bin, data; Npoints::Int=101, Nsamples::Int=1000)
    i = slice(data, :inten, bin)
    p = slice(data, :phase, bin)
    # errors
    δi = slice(data, :inten, bin; ind=3)
    δp = slice(data, :phase, bin; ind=3)
    # computation
    amps = constructamps(i,p)
    cosθv = range(-1,1, length=Npoints)
    calv = dNdcosθ.(cosθv; amps=amps)
    # bootstrap
    bstp = [let
                # i′ = [μ < σ ? μ : μ+randn()*σ for (μ,σ) in zip(i,δi)]
                i′ = i .+ randn(length(i)) .* δi
                i′ = map(v -> v < 0 ? 0.0 : v, i′)
                # i′ .*= sum(i) / sum(i′)
                p′ = p+randn()*δp
                dNdcosθ.(cosθv; amps=constructamps(i′,p′))
        end for _ in 1:Nsamples]
    qls = [quantile(getindex.(bstp,i), [0.16, 0.84]) for i in 1:length(bstp[1])]
    return getindex.(qls,1), getindex.(qls,2)
end
