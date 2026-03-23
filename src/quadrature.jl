## Quadrature ##
# this backend, unlike other backends, does not require an extension
# we use the generic conversions to deal with Infinity boxes, and then scale the domain
# to [-1,1] to use the quadrature weights

convert_domain(s::Domain.Box, ::Backend.Quadrature) =
    convert_domain_generic(s)

convert_integrand(i::Integral{<:Any}, ::Backend.Quadrature, domain, args; kw...) =
    convert_integrand_generic(i, domain, args; kw...)

(s::Backend.Quadrature)(f, domain, ::Nothing) =
    sum(((x,w),) -> f(x) * w, node_weight_iterator(s, domain...))

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

function node_weight_iterator(s::Backend.Quadrature, xmin::T, xmax::T) where {T<:Number}
    Δx = T((xmax - xmin) / 2)
    error_if_inf_Quadrature(Δx)
    itr = ((xmin + Δx * T(t + 1), Δx * T(w)) for (t, w) in zip(s.nodes, s.weights))
    return itr
end

function node_weight_iterator(s::Backend.Quadrature, mins::NTuple{N,T}, maxs::NTuple{N,T}) where {N,T<:Number}
    Δxs = T.((maxs .- mins) ./ 2)
    error_if_inf_Quadrature.(Δxs)
    Δxprod = prod(Δxs)
    zs = Iterators.product(ntuple(Returns(zip(s.nodes, s.weights)), Val(N))...)
    itr = ((mins .+ Δxs .* T.(first.(nws) .+ 1), Δxprod * T(prod(last.(nws)))) for nws in zs)
    return itr
end

error_if_inf_Quadrature(x) = isinf(x) && throw(ArgumentError("The Quadrature backend cannot deal with domains with `Inf`s. Use `Infinity` instead."))
