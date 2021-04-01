inlims(x,lims) = lims[1] ≤ x ≤ lims[2]

(..)(x::AbstractArray,i...) = getindex.(x,i...)
(..)(x::AbstractArray,i::Symbol) = getfield.(x,i)


function mysum(t::T where T<:Type, itrs::Base.Generator)
    v = zero(t)
    for it in itrs
        v += it
    end
    return v
end 

plotsfolder(tag...) = joinpath("plots", tag...)
fitsfolder( tag...) = joinpath("data", "exp_pro", tag...)
