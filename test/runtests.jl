using DoubleRegge
using Test

@testset "ϕ asymmetry" begin
    @test phi_asymmetry(ϕ->sin(ϕ)^2) < 1e-12
    @test phi_asymmetry(ϕ->sin(2ϕ)^2) < 1e-12
end