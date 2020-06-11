using DrWatson
@quickactivate "DoubleRegge"

using DoubleRegge
using QuadGK
using Plots
using Cuba
using Optim
theme(:wong; size=(500,350))
#
#
#
# data
const data = [(inten = get_intesity(L,M),
         phase = get_phase(L,M)) for (L,M) in LMs];
const xdata = data[1].inten[:,1];
const Nbins = length(xdata)
#
slice(data, property, bin; ind=2) = [v[bin, ind] for v in getproperty.(data, property)]
slice(data, bin) = getindex.(data, bin)
#
const amplitudes = constructamps(
    slice(data, :inten, 1:Nbins),
    slice(data, :phase, 1:Nbins)
);

# range
fit_range = (2.3, 3.0)
filt_fr = fit_range[1] .< xdata .< fit_range[2]
const amplitudes_fr = map(as->as[filt_fr], amplitudes)
const xdata_fr = xdata[filt_fr];

# fit
model(mηπ,cosθ,ϕ; pars) = pars[1]*modelS(pars[2:end],
    (s = DoubleRegge.s0, s1 = mηπ^2,
        cosθ = cosθ, ϕ = ϕ,
        t2 = -0.45))
#
function χ2(cosθ,ϕ,c)
    Id = abs2.(recamp(cosθ,ϕ, amplitudes_fr))
    Am = model.(xdata_fr,cosθ,ϕ; pars=c)
    Im = abs2.(Am)
    return sum(abs2, Id - Im)
end
#
# @time χ2(0.3,0.2,[1,1.1]) # test
#
integrate_dcosθdϕ(g) = cuhre((x,f)->f[1]=g(x),2,1)[1][1]*(4π)
χ2(pars) = integrate_dcosθdϕ(x->χ2(2x[1]-1,π*(2x[1]-1),pars))
#
# @time χ2([1,1.1]) # test
#
ft = Optim.minimizer(Optim.optimize(χ2, [1.1e-3,1.1], BFGS(),
               Optim.Options(show_trace = true)))
#
# @time quadgk(ϕ->model(xdata_fr[1], 0.2,ϕ; pars=ft), -π, π)[1] # test
# plotting
let
    cosθv = range(-1,1, length=101)
    ps = [let
        calv = dNdcosθ.(cosθv; amps=slice(amplitudes_fr, bin))
        plot(cosθv, calv, lab="")
        calv = [quadgk(ϕ->abs2(model(xdata_fr[bin], cosθ, ϕ; pars=ft)),-π, π)[1] for cosθ in cosθv]
        plot!(cosθv, calv, lab="")
    end for bin in 1:length(xdata_fr)]
    plot(ps..., size=(1500,1000))
end
