function phi_asymmetry(f)
    L = quadgk(f, π/2, 3π/2)[1]
    R = quadgk(f, -π/2, π/2)[1]
    return (L - R) / (L + R)
end
