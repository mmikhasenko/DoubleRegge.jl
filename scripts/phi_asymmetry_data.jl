using DrWatson
@quickactivate "DoubleRegge"

using DoubleRegge
using QuadGK
using Statistics
# 
using Plots
theme(:wong2;
    frame = :box, grid = false, lab = "",
    xlims = (:auto, :auto), ylims = (:auto, :auto))
##
const data_folder = "data/exp_raw/PLB_shifted"
const data = read_data(data_folder, description_ηπ);
const nBins = length(data)

let bin = 40
    cosθv = range(-1, 1, length = 101)
    calv = dNdcosθ.(cosθv, Ref(data.amps[bin]))
    plot(cosθv, calv)
end

phi_asymmetry_data(bin, cosθ; start = -π / 2) =
    phi_asymmetry(ϕ -> abs2(recamp(cosθ, ϕ, data.amps[bin])); start)
#
let
    plot()
    for bin in 40:3:55
        plot!(cosθ -> phi_asymmetry_data(bin, cosθ), -1:0.01:1, lab = "m(ηπ)=$(round(data.x[bin], digits=2))")
    end
    plot!(xlab = "cosθ", ylab = "Assymetry in position of ϕ peak")
end

phi_asymmetry_2d_data(bin) = phi_asymmetry_2d((cosθ, ϕ) -> abs2(recamp(cosθ, ϕ, data.amps[bin])))
@assert phi_asymmetry_2d_data(52) isa Float64

let
    scatter(bin -> phi_asymmetry_2d_data(bin), 1:5:55, lab = "")
    plot!(xlab = "bin: mηπ", ylab = "Assymetry in position of ϕ peak")
end
