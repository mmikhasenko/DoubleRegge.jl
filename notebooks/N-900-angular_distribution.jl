### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 0b87f2e0-b865-11ed-11dc-efc0e7cfec5f
# ╠═╡ show_logs = false
begin
	cd(joinpath(@__DIR__, ".."))
	using Pkg
	Pkg.activate(".")
	Pkg.instantiate()
	# 
	using DoubleRegge
	using QuadGK
	using Plots
	using DelimitedFiles
end

# ╔═╡ ff467cd0-a60f-44ce-886f-9494a73298b4
md"""
# Angular distribution

The partial-wave analysis of the $\eta^{(\prime)}\pi$ systems from the 2015 [COMPASS analysis](https://inspirehep.net/literature/1311486) parametrize the angular distributions for the system.

In this notebook, the distributions are computed for a selected bin.
"""

# ╔═╡ 4a22d70a-21dd-48f8-a039-ccc72969e587
theme(:vibrant, size=(500,350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    legendfontsize=9, legend =:topright,
    xlim=(:auto,:auto), ylim=(:auto,:auto))

# ╔═╡ 41338900-025d-457a-97ee-c65140bbf17f
setsystem!(:compass_ηπ)

# ╔═╡ 5acc90e9-443b-4777-8fcd-cf81130ac147
data = read_data(joinpath("data", "exp_raw", "PLB_shifted"), description_ηπ);

# ╔═╡ cfa8bfe0-2380-480f-9710-f6ebbb7ec429
bin = 5

# ╔═╡ fd1b5012-a2b3-4bd0-a447-9617f69a22d0
md"""
The asymmetry in the cosθ distribution is due to the interference of the odd and even waves.

 - It is strong in low energy, peaking forward, cosθ ~ 1 at `bin ~ 5`, and
 - peaking backward, cosθ ~ -1 for high energy, `bin ~ 50`
"""

# ╔═╡ c8697621-2454-4a09-b451-3b24d4b26bf3
let
	cosθv = range(-1,1, length=100)
	yv = dNdcosθ.(cosθv, Ref(data.amps[bin]))
	xtitle = round(data.x[bin], digits=2)
	plot(cosθv, yv, title="bin at $(xtitle) MeV",
		xlab="cosθ", ylab = "dN/dcosθ", lab="")
end

# ╔═╡ 27f953e6-7fa4-4331-9e27-dfc5327d6df9
let
	cosθv = range(-1,1, 200)
	ϕv = range(-π,π, 200)
	# 
	Av = recamp.(cosθv', ϕv, Ref(data.amps[bin]))
	Iv = abs2.(Av)
	xtitle = round(data.x[bin], digits=2)
	# 
	heatmap(cosθv, ϕv, Iv, title="bin at $(xtitle) MeV",
		xlab="cosθ", ylab = "ϕ", lab="")
end

# ╔═╡ bfc17bc8-971d-41db-843d-746ce898b109
md"""
The M=2 projection make ϕ distribution deviate from sin²(ϕ). Particularly, the peaks are slightly asymmetric, higher in the center for bin ~ 5.
"""

# ╔═╡ 3101cb9c-5ee8-4117-8a5e-dc4d025f8ef7
let 
	ϕv = range(-π,π, 100)
	yv = dNdϕ.(ϕv, Ref(data.amps[bin]); LMs=nothing)
	xtitle = round(data.x[bin], digits=2)
	plot(ϕv, yv, title="bin at $(xtitle) MeV",
		xlab="ϕ", ylab = "dN/dcosθ", lab="")
end

# ╔═╡ Cell order:
# ╟─ff467cd0-a60f-44ce-886f-9494a73298b4
# ╠═0b87f2e0-b865-11ed-11dc-efc0e7cfec5f
# ╠═4a22d70a-21dd-48f8-a039-ccc72969e587
# ╠═41338900-025d-457a-97ee-c65140bbf17f
# ╠═5acc90e9-443b-4777-8fcd-cf81130ac147
# ╠═cfa8bfe0-2380-480f-9710-f6ebbb7ec429
# ╠═27f953e6-7fa4-4331-9e27-dfc5327d6df9
# ╟─fd1b5012-a2b3-4bd0-a447-9617f69a22d0
# ╠═c8697621-2454-4a09-b451-3b24d4b26bf3
# ╟─bfc17bc8-971d-41db-843d-746ce898b109
# ╠═3101cb9c-5ee8-4117-8a5e-dc4d025f8ef7
