inlims(x,lims) = lims[1] ≤ x ≤ lims[2]

(..)(x::AbstractArray,i...) = getindex.(x,i...)
(..)(x::AbstractArray,i::Symbol) = getfield.(x,i)

reorder(vecofvec::Union{Vector{Vector{T}},Vector{SVector{N,T}}} where {T,N}) = [getindex.(vecofvec,i) for i in 1:length(vecofvec[1])]

function mysum(t::T where T<:Type, itrs::Base.Generator)
    v = zero(t)
    for it in itrs
        v += it
    end
    return v
end 

plotsfolder(tag...) = joinpath("plots", tag...)
fitsfolder( tag...) = joinpath("data", "exp_pro", tag...)


# phase alignment
function shiftbyperiod(Δ; period=2π)
    Nfull = div(Δ, period)
    Nhalfs = div(Δ-period*Nfull, period/2)
    # 
    -period * (Nfull + Nhalfs)
end
shiftbyperiod(el, pref; period=2π) = shiftbyperiod(el-pref; period=2π)
# 
meanshiftbyperiod(phases, ref=0) = phases .+ shiftbyperiod(mean(phases)-ref)


function alignperiodicsequence(v::Vector; period=2π)
    vc = copy(v)
    for i in 2:length(vc)
        vc[i] += shiftbyperiod(vc[i], vc[i-1]; period=period)
    end
    return vc
end
