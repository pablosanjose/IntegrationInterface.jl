## Quadrature ##
# this backend, unlike other backends, does not require an extension
# we use the generic conversions to deal with Infinity boxes, and then scale the domain
# to [-1,1] to use the quadrature weights

convert_domain(s::Domain.Box1D, ::Backend.Quadrature) =
    convert_domain_generic(s)

convert_integrand(i::Integral{<:Any,<:Backend.Quadrature}, domain, args; kw...) =
    convert_integrand_generic(i, domain, args; kw...)

function (s::Backend.Quadrature)(f, domain, ::Nothing)
    xmin, xmax = domain
    Δx = 0.5 * (xmax - xmin)
    result = sum(zip(s.nodes, s.weights); init = zero(float(Δx))) do (t, w)
        x = (xmin + Δx * (t + 1))
        return f(x) * w * Δx
    end
    return result
end

function (s::Backend.Quadrature)(f!, domain, result, args; kw...)
    xmin, xmax = domain
    Δx = 0.5 * (xmax - xmin)
    fill!(result, 0)
    out = similar(result)
    foreach(zip(s.nodes, s.weights)) do (t, w)
        x = (xmin + Δx * (t + 1))
        f!(out, x)
        @. result += out * w * Δx
        return nothing
    end
    return result
end
