inlims(x,lims) = lims[1] ≤ x ≤ lims[2]



function mysum(t::T where T<:Type, itrs::Base.Generator)
    v = zero(t)
    for it in itrs
        v += it
    end
    return v
end 