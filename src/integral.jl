## API ##
function integral(f::F, domain::AbstractDomain; result = nothing, backend::AbstractBackend = default_backend(domain)) where {F}
	return Integral(f, result, domain, backend)
end

# currying version
integral(domain; kw...) = f -> integral(f, domain; kw...)

# Can be overridden for user-defined f types and domains
default_backend(::Domain.Box1D) = Backend.QuadGK()
default_backend(::Domain.Sum{<:NTuple{<:Any,Domain.Box1D}}) = Backend.QuadGK()
default_backend(::Domain.Functional{<:Domain.Box1D}) = Backend.QuadGK()
default_backend(::Domain.Box) = Backend.HCubature()
default_backend(::Domain.Sum{<:NTuple{<:Any,Domain.Box}}) = Backend.HCubature()
default_backend(::Domain.Functional{<:Domain.Box}) = Backend.HCubature()
default_backend(d) = throw(ArgumentError("No default backend exists for domain $(domainname(d)), please specify one explicitly."))

ismutating(i::Integral) = !isnothing(i.result)

integrand(i::Integral) = i.integrand

backend(i::Integral) = i.backend

result(i::Integral) = i.result

backendname(i::Integral) = backendname(backend(i))
backendname(s::AbstractBackend) = nameof(typeof(s))
backendname(s::Type{<:AbstractBackend}) = nameof(s)

domain(i::Integral) = i.domain

domainname(i::Integral) = domainname(domain(i))

## call syntax (scalar and in-place) ##
(i::Integral)(args...; kw...) =
    integrate(i, evaluate_domain(domain(i), args; kw...), args; kw...)

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

## Extension fallbacks and generics ##
#   convert_domain, convert_integrand and backend(f, domain, result)

# generic domain conversions (extensions must opt-in to these explicitly)
#   transform_domain deals with Infinity and complex-to-real conversions - see changeofvariables.jl
convert_domain_generic(d::Domain.Box) = firstlast(transform_domain(d))

firstlast(d) = (first(d), last(d))

# generic integrand conversions (extensions must opt-in to these explicitly)
function convert_integrand_generic(i::Integral{Nothing}, domain, args; post = identity, kw...)
    f = integrand(i)
    f´(t) = post(jacobian(t, domain) * f(change_of_variables(t, domain)..., args...; kw...))
    return f´
end

function convert_integrand_generic(i::Integral, domain, args; post = identity, kw...)
    f! = integrand(i)
    function f!´(out, t)
        f!(out, change_of_variables(t, domain)..., args...; kw...)
        out .*= jacobian(t, domain)
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

## Show ##

Base.summary(i::Integral) = "Integral"

function Base.show(io::IO, J::Integral)
    i = get(io, :indent, "")
    ioindent = IOContext(io, :indent => i * "  ")
    print(io, summary(J), "\n",
"$i  Mutating   : $(ismutating(J))
$i  Domain     : $(domainname(J))
$i  Backend    : $(backendname(J))
$i  Integrand  : ")
  print(ioindent, integrand(J))
end
