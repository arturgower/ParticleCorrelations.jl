function trapezoidal_scheme(x::AbstractVector{T}; x0::T = first(x), xn::T = last(x)) where T<:AbstractFloat

    inds = axes(x,1)

    h = (x[inds[2]]-x[inds[1]])
    σs = similar(x)
    σs .= h

    σs[inds[1]] -= h/T(2)

    σs[end] -= h/T(2)
    # accounts for the ends
    σs[inds[1]] +=   (x[inds[1]] - x0)*(x[inds[2]] - x0 + h)/(2.0h)
    σs[inds[2]] +=  -(x0 - x[inds[1]])^T(2)/(T(2)*h)
    σs[end] +=  (xn - x[end])*(xn - x[end-1] + h)/(2h)
    σs[end-1] += - (xn - x[end])^T(2)/(T(2)*h)

    return σs
end