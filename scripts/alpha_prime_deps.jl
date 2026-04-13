using DrWatson
@quickactivate "DoubleRegge"

using DoubleRegge
using QuadGK
using Plots
using Cuba # HCubature # NIntegration
using Optim
using DelimitedFiles
theme(:wong; size=(500,350))
#
#
const reaction_system = compass_ηπ

function Model_of_cosθ(cosθ; mηπ = 2.2, α′ = 0.9)
    f(ϕ) = abs2(modelDR(α_ℙ, α_ℙ,
        (s = reaction_system.s0, s1 = mηπ^2, cosθ = cosθ, ϕ = ϕ, t2 = -0.2), reaction_system;
        η_forward = true, α′ = α′))
    return quadgk(f, -π, π)[1]
end


p1 = let
    plot(title = "varying α′")
    for α in 0.6:0.1:1.3
        cosθv = range(-1,1,length=100)
        calv = Model_of_cosθ.(cosθv; α′ = α)
        plot!(cosθv, calv/sum(calv), lab="α = $α")
    end
    plot!(leg=:left)
end
p2 = let
    plot(title = "same scan")
    for α in 0.6:0.1:1.3
        cosθv = range(-1,1,length=100)
        calv = Model_of_cosθ.(cosθv; α′ = α)
        plot!(cosθv,calv/sum(calv), lab="α = $α")
    end
    plot!(leg=:left)
end
p3 = let
    plot(title = "same scan")
    for α in 0.6:0.1:1.3
        cosθv = range(-1,1,length=100)
        calv = Model_of_cosθ.(cosθv; α′ = α)
        plot!(cosθv,calv/sum(calv), lab="α = $α")
    end
    plot!(leg=:left)
end
plot(p1,p2,p3, size=(1200,350), layout=grid(1,3))
