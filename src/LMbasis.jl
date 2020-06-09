
arg(z) = atan(imag(z), real(z))
Psi(L,M,cosθ,ϕ) = sqrt((2L+1)/(2π))*wignerd(L,M,0,cosθ)*sin(M*ϕ)

LMs = [
    (L = 1, M = 1), (L = 2, M = 1), (L = 2, M = 2),
    (L = 3, M = 1), (L = 4, M = 1),
    (L = 5, M = 1), (L = 6, M = 1)];
