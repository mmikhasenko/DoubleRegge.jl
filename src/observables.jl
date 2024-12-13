function phi_asymmetry(f)
    L = quadgk(f, π/2, 3π/2)[1]
    R = quadgk(f, -π/2, π/2)[1]
    return (L - R) / (L + R)
end

function phi_asymmetry_2d(f2d)
    L = cuhre((x,f)->f[1]=f2d(2x[1]-1, -π/2+π*x[2]), 2, 1)[1][1]
    R = cuhre((x,f)->f[1]=f2d(2x[1]-1,  π/2+π*x[2]), 2, 1)[1][1]
    return (L - R) / (L + R)
end
