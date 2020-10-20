
arg(z) = atan(imag(z), real(z))
Psi(L,M,cosθ,ϕ) = sqrt((2L+1)/(2π))*wignerd(L,M,0,cosθ)*sin(M*ϕ)

function pw_project(amplitude_cosθϕ, L, M)
    integrator_output = cuhre((x,f)->(f[1],f[2])=reim(amplitude_cosθϕ(2*x[1]-1,π*(2*x[2]-1))*Psi(L,M,2*x[1]-1,π*(2*x[2]-1))), 2, 2)
    pw = 4π*complex(integrator_output[1]...)
    return pw
end
