### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ 35bd41b0-9467-11eb-2e9d-f170d3275f02
begin
	using SymPy
	
	import PyCall
	PyCall.pyimport_conda("sympy.physics.wigner",       "sympy")
	PyCall.pyimport_conda("sympy.physics.quantum.spin", "sympy")
	
	import_from(sympy.physics.wigner)
    import_from(sympy.physics.quantum.spin, (:WignerD,), typ=:Any)
end

# ╔═╡ d226598b-7293-496d-a5b6-4d3987668ed8
begin
	using LinearAlgebra
	using Polynomials
	using Plots
end

# ╔═╡ 62ce2d64-b2eb-4c9d-912e-bdf0b2dbdd72
const ampls = [1,-10+1im,2.2,3.3im,-3.2,3.3]

# ╔═╡ b1100e02-cca9-4d30-aed9-5971cd3da8aa
md"""
## Polynomial coefficients of the Wigner d
figure out once and tabulate
"""

# ╔═╡ 621aeccb-91ba-4a35-b94b-694ed3184b65
θ, t = @vars θ t positive=true

# ╔═╡ 57185174-71e4-4164-bb65-46edc35b7681
dl(L) = sqrt(2L+Sym(1))*WignerD(L,1,0,0,θ,0).doit()

# ╔═╡ dbda6de0-bd51-426b-8c9a-5402510fb5e3
function half_angle_substitution(L; Lmax=6)
	expr = sympy.expand_trig(dl(L))
	poly = simplify(expr.subs(sin(θ), 2t/(1+t^2)).subs(cos(θ), (1-t^2)/(1+t^2)) * (1+t^2)^Lmax)/t
	poly.subs(t^2,t)
end

# ╔═╡ 764ec4b8-2a39-4f5c-968f-9e5df718384c
half_angle_substitution(2)

# ╔═╡ c65e6e8c-77f2-4b44-9ac1-881f8cfa9776
function tanθh_coeff_wignerd(L; Nmax=5)
	poly = expand(half_angle_substitution(L; Lmax=Nmax+1))
	[poly.coeff(t, n) for n in 0:Nmax]
end

# ╔═╡ 79e93062-ac1c-4704-99e0-3302263301fd
cs = [tanθh_coeff_wignerd(L; Nmax=5) for L in 1:6];

# ╔═╡ 10f1a80e-d1fb-4de3-b6fa-ed5823d1f24d
hcat(cs...)

# ╔═╡ 483b0a8c-fdcb-4e1b-a1fb-f9552262f6de
md"""
Use sympy command to print the strings to be directly copied to the code
```julia
sympy.ccode.(cs)
``` 
"""

# ╔═╡ 0963c797-2ca8-4da7-8090-0a7831fb8cd2
const fcs = map(x->Float64.(x), cs)

# ╔═╡ 0094a836-a656-435b-bbbc-93d143a7b33b
md"""
## Amplitude to polynomial and roots
"""

# ╔═╡ a84d534f-d096-4c06-8666-57b3a028210b
function conjugate_subvector(v, indices)
	vnew = copy(v)
	vnew[indices] .= conj(vnew[indices])
	return vnew
end

# ╔═╡ 38879339-9302-4813-af1a-6caefffc571f
conjindices(i, N) = digits(i, base=2, pad=N) .!= 0

# ╔═╡ 4aa5c9e6-81cb-45eb-b1a2-bbf824a7ab42
get_coeff_from_roots(rs) = coeffs(prod(Polynomial([-r, 1]) for r in rs))

# ╔═╡ 815bea3a-2725-496e-b53d-e5c86ae54629
const coeff0 = sum(ampls .* fcs)

# ╔═╡ b1a557f2-099c-4a1c-badf-ba9285850984
const roots0 = Polynomials.roots(Polynomial(coeff0))

# ╔═╡ 71cd9f0c-2e37-40a3-9f42-43c6a251be4b
function amplitude_from_roots(rs)
	coeff = get_coeff_from_roots(rs)
	return vcat(hcat(fcs...) \ coeff)
end

# ╔═╡ fac5cc90-1a60-461d-8a7f-fd3c24628be1
rootsi(i) = conjugate_subvector(roots0, conjindices(i,length(roots0)))

# ╔═╡ c4dd64be-7bf8-41c1-9a79-176bbb88c8d5
all_amplitudes = [amplitude_from_roots(rootsi(i))*coeff0[end] for i in 0:2^5-1]

# ╔═╡ 78f87724-e3c0-40f2-a8b8-82c9071b896b
sum.(abs2,all_amplitudes)

# ╔═╡ 2d714e45-ecda-440d-8e87-b90083b0d83b
plot(map(x->abs2.(x),all_amplitudes), lab="", )

# ╔═╡ 323159f7-7a12-442c-b18b-a065d88a8340
plot(layout=grid(6,1), size=(500,1000),
	histogram(abs2.(getindex.(all_amplitudes,1)), bins=1000),
	histogram(abs2.(getindex.(all_amplitudes,2)), bins=1000),
	histogram(abs2.(getindex.(all_amplitudes,3)), bins=1000),
	histogram(abs2.(getindex.(all_amplitudes,4)), bins=1000),
	histogram(abs2.(getindex.(all_amplitudes,5)), bins=1000),
	histogram(abs2.(getindex.(all_amplitudes,6)), bins=1000))

# ╔═╡ 9a686395-a2ba-4be9-8e82-58ccb9bb4f97
md"""
#### Tests
"""

# ╔═╡ c029b3c3-1399-47bb-836e-da1d057a76d3
testA(th,PWs) = abs2(sum(dl(L)*a for (L,a) in enumerate(PWs))(th))

# ╔═╡ 13450685-1581-41d7-9eb9-a57c5cc7dedd
function poldl(θ, l, coeff; Lmax=6)
	t = tan(θ/2)
	pol = Polynomial(vcat([[f,0] for f in coeff]...))*Polynomial([0,1])
	pol(t)/(1+t^2)^Lmax
end

# ╔═╡ b04400d9-0521-40f9-8a40-85ef9729000c
let L=3, th = 1.1
 	poldl(th, L, fcs[L]), N(dl(L)(th))
end

# ╔═╡ 9247ef50-71cf-423c-8ce5-3c47f51e1746
poltext(t,rs) = abs2(prod(t-r for r in rs))

# ╔═╡ 56bafe07-f4a7-4cbf-b3d5-9e9a6a1279d8
poltext(1.3,roots0)

# ╔═╡ e93722d3-db82-4565-8949-ae9eebc805f8
poltext(1.3,rootsi(31))

# ╔═╡ 43493eda-0e10-49a0-ac83-853fe46ecd3f
bar([real.(roots0) imag.(roots0)], lab=["real" "imag"])

# ╔═╡ 90262948-4417-4c3d-8091-5b8752610020
amplitude_from_roots(roots0)

# ╔═╡ 291c8818-802c-48f9-b1a4-7c6dec7524b4
amplitude_from_roots(
	conjugate_subvector(roots0, conjindices(0,length(roots0)).!= 0))

# ╔═╡ 7f111ea7-7889-4046-a14b-856eabdc2ed8
begin
	new_ampl = vcat(hcat(fcs...) \ coeff0)
	new_ampl .≈ ampls
end

# ╔═╡ 522e379d-f735-4314-8fce-e745afd9af73
Hij = [dot(v./norm(v), h./norm(h)) for v in fcs, h in fcs]

# ╔═╡ c1b4732a-568b-4508-baf2-e1ee98d21522
det(Hij)

# ╔═╡ d176cf44-b259-46a1-8f12-6415251148e6
prod(abs.(get_coeff_from_roots(roots0)*coeff0[end] - coeff0) .< 1e13)

# ╔═╡ 88578939-6efc-4908-88e3-07269c5ad86a
[sympy.integrate(dl(L)^2*sin(θ), (θ, 0, PI)) for L in 1:6]

# ╔═╡ ea1a41a2-1210-4f9e-8692-cf5c1ecb33b4
md"""
## debugging area
""" 

# ╔═╡ 28e4f5e6-a1cb-469c-a26d-f9724bec37e4
scatter([Polynomial(sum(amplitude_from_roots(rootsi(i)) .* fcs))(1.1) for i in 0:31])

# ╔═╡ 80e236c1-c255-41cb-8729-8e16ae3f0675
Float64(testA(1.1,ampls))

# ╔═╡ 88d3aa1d-86bf-4d23-bb4c-8744cb95ce8c
Float64(testA(1.1,amplitude_from_roots(rootsi(2))))

# ╔═╡ 2634b211-0569-4d34-a6c0-992297a387ee
abs2(poldl(1.1, 6, sum(amplitude_from_roots(rootsi(0)) .* fcs)))

# ╔═╡ 60f23099-4adb-4aac-9b8d-75a194b13e30
sum(abs2,amplitude_from_roots(rootsi(2)))

# ╔═╡ c58bfb8a-03c0-4387-843f-3376f5cb1cf0
amplitude_from_roots(rootsi(2))

# ╔═╡ Cell order:
# ╠═35bd41b0-9467-11eb-2e9d-f170d3275f02
# ╠═d226598b-7293-496d-a5b6-4d3987668ed8
# ╠═62ce2d64-b2eb-4c9d-912e-bdf0b2dbdd72
# ╟─b1100e02-cca9-4d30-aed9-5971cd3da8aa
# ╠═621aeccb-91ba-4a35-b94b-694ed3184b65
# ╠═57185174-71e4-4164-bb65-46edc35b7681
# ╠═dbda6de0-bd51-426b-8c9a-5402510fb5e3
# ╠═764ec4b8-2a39-4f5c-968f-9e5df718384c
# ╠═c65e6e8c-77f2-4b44-9ac1-881f8cfa9776
# ╠═79e93062-ac1c-4704-99e0-3302263301fd
# ╠═10f1a80e-d1fb-4de3-b6fa-ed5823d1f24d
# ╟─483b0a8c-fdcb-4e1b-a1fb-f9552262f6de
# ╠═0963c797-2ca8-4da7-8090-0a7831fb8cd2
# ╟─0094a836-a656-435b-bbbc-93d143a7b33b
# ╠═a84d534f-d096-4c06-8666-57b3a028210b
# ╠═38879339-9302-4813-af1a-6caefffc571f
# ╠═4aa5c9e6-81cb-45eb-b1a2-bbf824a7ab42
# ╠═815bea3a-2725-496e-b53d-e5c86ae54629
# ╠═b1a557f2-099c-4a1c-badf-ba9285850984
# ╠═71cd9f0c-2e37-40a3-9f42-43c6a251be4b
# ╠═fac5cc90-1a60-461d-8a7f-fd3c24628be1
# ╠═c4dd64be-7bf8-41c1-9a79-176bbb88c8d5
# ╠═78f87724-e3c0-40f2-a8b8-82c9071b896b
# ╠═2d714e45-ecda-440d-8e87-b90083b0d83b
# ╠═323159f7-7a12-442c-b18b-a065d88a8340
# ╟─9a686395-a2ba-4be9-8e82-58ccb9bb4f97
# ╠═c029b3c3-1399-47bb-836e-da1d057a76d3
# ╠═13450685-1581-41d7-9eb9-a57c5cc7dedd
# ╠═b04400d9-0521-40f9-8a40-85ef9729000c
# ╠═9247ef50-71cf-423c-8ce5-3c47f51e1746
# ╠═56bafe07-f4a7-4cbf-b3d5-9e9a6a1279d8
# ╠═e93722d3-db82-4565-8949-ae9eebc805f8
# ╠═43493eda-0e10-49a0-ac83-853fe46ecd3f
# ╠═90262948-4417-4c3d-8091-5b8752610020
# ╠═291c8818-802c-48f9-b1a4-7c6dec7524b4
# ╠═7f111ea7-7889-4046-a14b-856eabdc2ed8
# ╠═522e379d-f735-4314-8fce-e745afd9af73
# ╠═c1b4732a-568b-4508-baf2-e1ee98d21522
# ╠═d176cf44-b259-46a1-8f12-6415251148e6
# ╠═88578939-6efc-4908-88e3-07269c5ad86a
# ╟─ea1a41a2-1210-4f9e-8692-cf5c1ecb33b4
# ╠═28e4f5e6-a1cb-469c-a26d-f9724bec37e4
# ╠═80e236c1-c255-41cb-8729-8e16ae3f0675
# ╠═88d3aa1d-86bf-4d23-bb4c-8744cb95ce8c
# ╠═2634b211-0569-4d34-a6c0-992297a387ee
# ╠═60f23099-4adb-4aac-9b8d-75a194b13e30
# ╠═c58bfb8a-03c0-4387-843f-3376f5cb1cf0
