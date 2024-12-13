using DoubleRegge
using Test

t0 = TwoBodyPartialWaves([(1,1),(2,1),(3,1)], [1, 2-3im, -4-4im])

@testset "TwoBodyPartialWaves: changerepresentation" begin
    @test prod(changerepresentation(changerepresentation(t0; iref=1)).PWs .≈ t0.PWs)
end

@testset "ϕ asymmetry" begin
    @test phi_asymmetry(ϕ->sin(ϕ)^2) < 1e-12
    @test phi_asymmetry(ϕ->sin(2ϕ)^2) < 1e-12
end

@testset "alignment" begin
    shiftbyperiod(0.999, 0; period=2) == 0
    shiftbyperiod(1.999, 0; period=2) == -2
    shiftbyperiod(2.9,0; period=2) == -2
    shiftbyperiod(3.1,0; period=2) == -4
    # 
    testseq = range(0.1, 0.8, length=10)
    spoiled = testseq .+ 2*[0, rand(-1:1, 9)...]
    cured = alignperiodicsequence(spoiled; period=2)
    prod(testseq .≈ cured)
end

