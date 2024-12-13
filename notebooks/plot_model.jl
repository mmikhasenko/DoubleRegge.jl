### A Pluto.jl notebook ###
# v0.19.47

using Markdown
using InteractiveUtils

# ╔═╡ 64d5bf96-d096-4864-9fbd-67df6d7aeb5e
# ╠═╡ show_logs = false
begin
	using Pkg
	cd(joinpath(@__DIR__, ".."))
	Pkg.activate(".")
	# 
	using TOML
	# 
	using Plots
	import Plots.PlotMeasures.mm
	using LaTeXStrings
	# 
	using DoubleRegge
	using Parameters
end

# ╔═╡ d9173b8c-2863-4726-820f-01970bba84f0
theme(:wong2;
	size = (500, 350), xlab = L"m_{\eta'\pi}\,\,(\mathrm{GeV})",
	frame=:box, grid=false, lab="",
	xlims=(:auto, :auto), ylims=(:auto, :auto))

# ╔═╡ a93d70c8-50bb-4836-8a13-5d69b824945c
md"""
## Model to pick
"""

# ╔═╡ 42ede67c-de09-4edc-b05f-c8bdc6205b74
@unpack settings, fit_results = let
	tag = "etappi_a2Po-a2f2-f2f2-PoPo_opposite-sign"
	settings_file = fitsfolder(tag, "fit-results.toml")
	!isfile(settings_file) && error("no file")
	TOML.parsefile(settings_file)
end

# ╔═╡ f0a7b1b4-b971-11ef-0c7f-3b497c108dc6
setsystem!(Symbol(settings["system"]))

# ╔═╡ d900c75f-23be-4000-a365-a4bd7706ef7b
const model = build_model(
    sixexchages[settings["exchanges"]],
    settings["t2"],
	settings["scale_α"]);

# ╔═╡ cb933425-f0e9-49ae-8aac-00d17f236525
const fixed_pars = fit_results["fit_minimizer"];

# ╔═╡ 6ad94d49-2e9f-43a6-ae50-c3f06809310b
fixed_model(m, cosθ, ϕ) = model(m, cosθ, ϕ; pars = fixed_pars)

# ╔═╡ bf32445c-fbfd-4809-8452-f1e17394d4a3
intensity(m, cosθ, ϕ) = abs2(fixed_model(m, cosθ, ϕ)) * q(m)

# ╔═╡ 9cf4002c-c678-40e2-9520-46cb5443e2fd
plot(cosθ -> intensity(2.2, cosθ, 0.3), -1, 1)

# ╔═╡ Cell order:
# ╠═64d5bf96-d096-4864-9fbd-67df6d7aeb5e
# ╠═d9173b8c-2863-4726-820f-01970bba84f0
# ╟─a93d70c8-50bb-4836-8a13-5d69b824945c
# ╠═42ede67c-de09-4edc-b05f-c8bdc6205b74
# ╠═f0a7b1b4-b971-11ef-0c7f-3b497c108dc6
# ╠═d900c75f-23be-4000-a365-a4bd7706ef7b
# ╠═cb933425-f0e9-49ae-8aac-00d17f236525
# ╠═6ad94d49-2e9f-43a6-ae50-c3f06809310b
# ╠═bf32445c-fbfd-4809-8452-f1e17394d4a3
# ╠═9cf4002c-c678-40e2-9520-46cb5443e2fd
