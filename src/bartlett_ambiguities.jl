const Lmax = 6

const basis = [
	[   -sqrt(6),  -5*sqrt(6), -10*sqrt(6), -10*sqrt(6),  -5*sqrt(6),    -sqrt(6)],
	[  -sqrt(30), -3*sqrt(30), -2*sqrt(30),  2*sqrt(30),  3*sqrt(30),    sqrt(30)],
	[-2*sqrt(21),           0, 10*sqrt(21), 10*sqrt(21),           0, -2*sqrt(21)],
	[ -6*sqrt(5),  24*sqrt(5),  30*sqrt(5), -30*sqrt(5), -24*sqrt(5),   6*sqrt(5)],
	[   -sqrt(330),   9*sqrt(330), -10*sqrt(330), -10*sqrt(330),   9*sqrt(330),    -sqrt(330)],
	[   -sqrt(546),  15*sqrt(546), -50*sqrt(546),  50*sqrt(546), -15*sqrt(546),     sqrt(546)]]


function conjugate_subvector(v, indices::BitArray)
	vnew = copy(v)
	vnew[indices] .= conj(vnew[indices])
	return vnew
end
conjindices(i, N) = digits(i, base=2, pad=N) .!= 0
conjugate_subvector(v, i::Int) =
	conjugate_subvector(v, conjindices(i,length(v)))

get_coeff_from_roots(rs) = coeffs(prod(Polynomial([-r, 1]) for r in rs))

function amplitude_from_roots(rs)
	coeff = get_coeff_from_roots(rs)
	return vcat(hcat(basis...) \ coeff)
end

function bartlettambiguities(PWs6)
	# 
	coeff0 = sum(PWs6 .* basis)
	roots0 = Polynomials.roots(Polynomial(coeff0))
	#
	all_amplitudes = [amplitude_from_roots(
		conjugate_subvector(roots0,i)) * coeff0[end]
			for i in 0:2^(Lmax-1)-1]
	# 
	return all_amplitudes
end