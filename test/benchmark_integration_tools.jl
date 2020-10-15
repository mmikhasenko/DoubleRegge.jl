using BenchmarkTools

using HCubature
using NIntegration
using QuadGK
using Cuba

function model(cosθ,ϕ; αp1 = 2.3)
    α1 = α2 = αp1-1.0
    v = (exp(-α1*(cosθ+1)) + exp(-α2*(1-cosθ)))*(1+sin(ϕ))^(1.2)
    return abs2(v)
end
#
# using Plots
# plot(cosθ -> model(cosθ,π/4; αp1 = 2.3), -1, 1)
#
# QuadGK
I_quadgk(αp1) = quadgk(cosθ->quadgk(ϕ->model(cosθ,ϕ;αp1=αp1),-π,π)[1],-1.0, 1.0)[1]
#
# NIntegration
I_nintegration(αp1) = nintegrate((cosθ,ϕ,z)->model(cosθ,ϕ;αp1=αp1), (-1.0, -π, 0.0), (1.0, π, 1.0))[1]
#
# HCubature
I_hcubature(αp1) =  hcubature(x->model(x[1],x[2]; αp1 = αp1), [-1.0, -π], [1.0, π])[1]
#
# Cuba
I_cuba(αp1) = cuhre((x,f)->f[1]=model(2x[1]-1, π*(2x[2]-1);αp1=αp1),2,1).integral[1]*(4π)
#
@btime I_cuba(2.3) # test
@btime I_hcubature(2.3) # test
@btime I_nintegration(2.3) # test
@btime I_quadgk(2.3) # test
