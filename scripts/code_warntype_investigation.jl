
arg(z) = atan(imag(z), real(z))
Psi(L, M, cosθ, ϕ) = sqrt((2L + 1) / (2π)) * wignerd(L, M, 0, cosθ) * sin(M * ϕ)

function pw_project(amplitude_cosθϕ, L, M)
    integrator_output = cuhre((x, f) -> (f[1], f[2]) = reim(amplitude_cosθϕ(2 * x[1] - 1, π * (2 * x[2] - 1)) * Psi(L, M, 2 * x[1] - 1, π * (2 * x[2] - 1))), 2, 2)
    pw = 4π * complex(integrator_output[1]...)
    return pw
end

@code_warntype model(2.1, 0.3, 0.3; pars = [1, 1, 1.0]) # fine, body::Complex{Float64}
#
x(a, b, c) = model(a, b, c; pars = [1, 1, 1.0])
@code_warntype x(2.1, 0.3, 0.3) # body::Any


#      _|_|  _|                        _|      
#    _|          _|  _|_|    _|_|_|  _|_|_|_|  
#  _|_|_|_|  _|  _|_|      _|_|        _|      
#    _|      _|  _|            _|_|    _|      
#    _|      _|  _|        _|_|_|        _|_|  

module MyTest
using PartialWaveFunctions
Psi(L, M, cosθ, ϕ) = sqrt((2L + 1) / (2π)) * wignerd(L, M, 0, cosθ) * sin(M * ϕ)
recamp(cosθ, ϕ, amps, LMs) = sum(a * Psi(L, M, cosθ, ϕ) for (a, (L, M)) in zip(amps, LMs))
export recamp, Psi
end
# 
using .MyTest

const testLM = [(L = 1, M = 1), (L = 3, M = 1)]
const testA = rand(Complex{Float64}, 2)
#
@code_warntype Psi(3, 1, 0.3, 0.3) # ::Float64
recamp(0.3, 0.3, testA, testLM)
@code_warntype recamp(0.3, 0.3, testA, testLM) # ::Any

function test()
    return recamp(0.3, 0.3, testA, testLM)
end
@code_warntype test() # Body::Float64


#                                                          _|  
#    _|_|_|    _|_|      _|_|_|    _|_|    _|_|_|      _|_|_|  
#  _|_|      _|_|_|_|  _|        _|    _|  _|    _|  _|    _|  
#      _|_|  _|        _|        _|    _|  _|    _|  _|    _|  
#  _|_|_|      _|_|_|    _|_|_|    _|_|    _|    _|    _|_|_|  



module MyTest2
using Parameters
function mysum(t::T where T <: Type, itrs::Base.Generator)
    v = zero(t)
    for it in itrs
        v += it
    end
    return v
end
function build_model(list_of_settings)
    function model(x; pars)
        return mysum(typeof(1im * pars[1]), (p * single_term_model(x; settings = s)
                                             for (p, s) in zip(pars, list_of_settings)))
        # return sum(p*single_term_model(x; settings=s)
        #     for (p,s) in zip(pars, list_of_settings))
    end
    return model
end
function single_term_model(x; settings)
    @unpack N, k, p = settings
    value = sum((k * x)^i * cis(i * p * x) for i in 1:N)
    return value
end
export build_model
end

using .MyTest2


# fine, body::Complex{Float64}
# 
const list_of_settings = [(N = 3, k = 0.1, p = 0.3),
    (N = 5, k = 0.5, p = 0.3),
    (N = 7, k = 0.2, p = 0.3)]
#
#!!!!!!!!!!
# MyTest2.single_term_model(2.2; settings=list_of_settings[1])
#!!!!!!!!!!
# 
const model = build_model(list_of_settings)
model(2.2; pars = [1.0, 1.0, 1.0]) # Body::Any
@code_warntype model(2.2; pars = [1.0, 1.0, 1.0]) # Body::Any


M(x) = model(x; pars = [1, 1, 1.0])
M(2.2)
@code_warntype M(2.2) # Body::Any

function test()
    return M(2.2)
end
@code_warntype test() # Body::Any





module MyTest2
using Parameters
function build_model(list_of_settings)
    function model(x; pars)
        return sum(p * modelA(x; settings = s)
                   for (p, s) in zip(pars, list_of_settings))
    end
    return model
end
function modelA(x; settings)
    @unpack N, k, p = settings
    value = sum((k * x)^i * cis(i * p * x) for i in 1:N)
    return value
end

export build_model
end

using .MyTest2

const list_of_settings = [(N = 3, k = 0.1, p = 0.3),
    (N = 5, k = 0.5, p = 0.3),
    (N = 7, k = 0.2, p = 0.3)]
# 
const model = build_model(list_of_settings) # added const
MyTest2.modelA(2.1; settings = list_of_settings[1])
# 
model(2.2; pars = [1, 1, 1.0])
@code_warntype model(2.2; pars = [1, 1, 1.0]) # fine, body::Complex{Float64}
#
M(x) = model(x; pars = [1, 1, 1.0])
@code_warntype M(2.2) # still Body::Any







M(x) = model(x; pars = [1, 1, 1.0])
M(2.2)
@code_warntype M(2.2) # Body::Any

function test()
    return M(2.2)
end
@code_warntype test() # Body::Any



const model = x -> 3x
@code_warntype model(6)  # Body::Int64
modelcall() = model(7)
@code_warntype modelcall()  # Body::Any
#
function test()
    modelcall()
end
@code_warntype test() # Body::Any
