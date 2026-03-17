## Quadrature ##
# this backend, unlike other backends, does not require an extension
# we use the generic conversions to deal with Infinity boxes, and then scale the domain
# to [-1,1] to use the quadrature weights

convert_domain(s::Domain.Box, ::Backend.Quadrature) =
    convert_domain_generic(s)

convert_integrand(i::Integral{<:Any,<:Backend.Quadrature}, domain, args; kw...) =
    convert_integrand_generic(i, domain, args; kw...)

(s::Backend.Quadrature)(f, domain, ::Nothing) =
    sum(((x,w),) -> f(x) * w, node_weight_iterator(s, domain...); init = 0.0)

function (s::Backend.Quadrature)(f!, domain, result; kw...)
    itr = node_weight_iterator(s, domain...)
    fill!(result, 0)
    out = similar(result)
    foreach(itr) do (x, w)
        f!(out, x)
        @. result += out * w
        return nothing
    end
    return result
end

function node_weight_iterator(s::Backend.Quadrature, xmin::NumberOrInfinity, xmax::NumberOrInfinity)
    Δx = 0.5 * (xmax - xmin)
    error_if_inf(Δx)
    itr = ((xmin + Δx * (t + 1), Δx * w) for (t, w) in zip(s.nodes, s.weights))
    return itr
end

function node_weight_iterator(s::Backend.Quadrature, mins::NTuple{N,NumberOrInfinity}, maxs::NTuple{N,NumberOrInfinity}) where {N}
    Δxs = 0.5 .* (mins .- maxs)
    error_if_inf.(Δxs)
    Δxprod = prod(Δxs)
    zs = Iterators.product(ntuple(Returns(zip(s.nodes, s.weights)), Val(N))...)
    itr = ((mins .+ Δxs .* (first.(nws) .+ 1), Δxprod * prod(last.(nws))) for nws in zs)
    return itr
end

error_if_inf(x) = isinf(x) && throw(ArgumentError("The Quadrature backend cannot deal with domains with `Inf`s. Use `Infinity` instead."))
