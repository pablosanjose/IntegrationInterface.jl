
### Infinity API ###
# Deal with Infinity and complex domain bounds.
# See https://github.com/pablosanjose/IntegrationInterface.jl/issues/3 for details

point(d::Infinity) = d.point
point(x::Number) = x
point(x::Tuple) = x

Base.:+(d::Infinity) = d
Base.:-(d::Infinity) = Infinity(-point(d))

### Integral/integral API ###

(i::Integral)(args...; kw...) =
    integrate_or_zero(i, evaluate_domain(domain(i), args; kw...), args; kw...)

# fastpath for zero-size domains
function integrate_or_zero(i, domain::AbstractEvaluatedDomain{<:Any,T}, args; kw...) where {T}
    if isempty(domain)  # fastpath to zero
        f = integrand(i)
        x0 = first(domain)
        f0 = f(x0..., args...; kw...) * zero(T)
        return f0
    else
        return integrate(i, domain, args; kw...)
    end
end

# integrate can assume the domain is not a Domain.Functional.
# We just need backend-specific conversions
integrate(i::Integral, domain, args; kw...) =
    i.backend(convert_integrand(i, domain, args; kw...), convert_domain(domain, i.backend), i.result)

# Any Domain Sum is handled by summing over the domains ##
# non-mutating version
function integrate(i::Integral{Nothing}, domain::Domain.Sum, args; kw...)
    result = sum(ungroup(domain)) do subdomain
        integrate(i, evaluate_domain(subdomain, args; kw...), args; kw...)
    end
    return result
end

# mutating version
function integrate(i::Integral, domain::Domain.Sum, args; kw...)
    resultsum = zero(result(i))
    foreach(ungroup(domain)) do subdomain
        resultsum .+= integrate(i, evaluate_domain(subdomain, args; kw...), args; kw...)
    end
    return resultsum
end

evaluate_domain(d::Domain.Functional, args; kw...) = d(args...; kw...)
evaluate_domain(d::AbstractDomain, args; kw...) = d

## Translation interface: how to map our API to backends, including Infinity handling ##
#   convert_domain: domain passed to backend -- convert_domain_generic is a useful default
#       transform_domain: domain after change of variables. Called by convert_domain_generic
#   convert_integrand: integrand passed to backend
#       change_of_variables: map coordinates of transform_domain to original, with jacobian
#   backend(f, domain, result): implemented in extensions, where f and domain are converted
#   --
#   See changeofvariables.jl for transformations for different domains

convert_domain_generic(d::Domain.Box) = firstlast(transform_domain(d))

firstlast(d) = (first(d), last(d))

# generic integrand conversions (extensions must opt-in to these explicitly)
function convert_integrand_generic(i::Integral{Nothing}, domain, args; post = identity, kw...)
    f = integrand(i)
    function f´(t)
        x, dx = change_of_variables(t, domain)
        return post(dx * f(x..., args...; kw...))
    end
    return f´
end

function convert_integrand_generic(i::Integral, domain, args; post = identity, kw...)
    f! = integrand(i)
    function f!´(out, t)
        x, dx = change_of_variables(t, domain)
        f!(out, x..., args...; kw...)
        out .*= dx
        maybe_post!(out, post)
        return out
    end
    return f!´
end

maybe_post!(out, ::typeof(identity)) = out
maybe_post!(out, post) = (out .= post.(out))

# failures
convert_domain(d::AbstractDomain, s::AbstractBackend) =
    throw(ArgumentError("No conversion method for domain $(domainname(d)) defined for $(backendname(s)) backend, or backend package not loaded."))

convert_integrand(i::Integral, domain, args; kw...) =
    throw(ArgumentError("No conversion method for the integrand defined for the $(backendname(i)) backend, or backend package not loaded."))

(s::AbstractBackend)(f, domain, result) =
    error("The integration backend package for $(nameof(typeof(s))) is not loaded.")
