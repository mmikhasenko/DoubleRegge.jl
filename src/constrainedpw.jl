# 
fold(longp) = (Np = div(length(longp),2); longp[1:Np] .+ 1im .* longp[(Np+1):end])
unfold(p) = vcat(real.(p), imag.(p))
# test
(v = rand(Complex{Float64},10); prod(fold(unfold(v)) .== v))

function constrained_pw_projection(intensity_cosθϕ, init_pars, LMs)
    
    function integrand(cosθ,ϕ,pars)
        Id = intensity_cosθϕ(cosθ, ϕ)
        Im = abs2(sum(p*Psi(L,M,cosθ,ϕ) for (p,(L,M)) in zip(pars,LMs)))
        Im ≈ 0.0 && (Im=nextfloat(0.0))
        return -Id*log(Im)
    end
    f(pars) = sum(abs2, pars) + integrate_dcosθdϕ((cosθ,ϕ)->integrand(cosθ,ϕ,fold(pars)))[1]
    #
    f = Optim.optimize(f, unfold(init_pars), BFGS(), # f′!, 
                Optim.Options(show_trace = true, f_tol=1e-7, iterations=100))
    return (min = f.minimum, pars = fold(f.minimizer))
end 

function constrained_pw_projection_with_derivative(intensity_cosθϕ, init_pars, LMs)
    
    function integrand(cosθ,ϕ,pars)
        Id = intensity_cosθϕ(cosθ, ϕ)
        Im = abs2(sum(p*Psi(L,M,cosθ,ϕ) for (p,(L,M)) in zip(pars,LMs)))
        Im ≈ 0.0 && (Im=nextfloat(0.0))
        return -Id*log(Im)
    end
    f(pars) = sum(abs2, pars) + integrate_dcosθdϕ((cosθ,ϕ)->integrand(cosθ,ϕ,fold(pars)))[1]

    # authomatic derivative
    function f′(pars)
        dint(cosθ,ϕ) = ForwardDiff.gradient(
                p->sum(abs2, p) / (4π) + integrand(cosθ,ϕ, fold(p)),
            unfold(pars))
        integral = integrate_dcosθdϕ(dint; dims=2*length(pars))
        fold(integral)
    end
    f′!(stor,pars) = copyto!(stor, f′(pars))
    #
    f = Optim.optimize(f, f′!, unfold(init_pars), BFGS(), # f′!, 
                Optim.Options(show_trace = true, iterations=15))
    return (min = f.minimum, pars = fold(f.minimizer))
end 
