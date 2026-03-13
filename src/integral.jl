## API ##
function integral(f::F, domain::AbstractDomain; result = missing, solver::AbstractBackend = default_solver(domain)) where {F}
	return Integral(f, result, domain, solver)
end

# currying version
integral(domain; kw...) = f -> integral(f, domain; kw...)

# Can be overridden for user-defined f types and domains
default_solver(::Domain.Segment) = Backend.QuadGK()
default_solver(::Domain.Sum{<:NTuple{<:Any,Domain.Segment}}) = Backend.QuadGK()
default_solver(::Domain.Functional{<:Domain.Segment}) = Backend.QuadGK()
default_solver(::Domain.Box) = Backend.HCubature()
default_solver(::Domain.Sum{<:NTuple{<:Any,Domain.Box}}) = Backend.HCubature()
default_solver(::Domain.Functional{<:Domain.Box}) = Backend.HCubature()
default_solver(d) = throw(ArgumentError("No default solver exists for domain $(domainname(d))"))

ismutating(i::Integral) = !ismissing(i.result)

integrand(i::Integral) = i.integrand

solver(i::Integral) = i.solver

result(i::Integral) = i.result

solvername(i::Integral) = solvername(solver(i))
solvername(s::AbstractBackend) = nameof(typeof(s))
solvername(s::Type{<:AbstractBackend}) = nameof(s)

domain(i::Integral) = i.domain

domainname(i::Integral) = domainname(domain(i))

## call syntax (scalar and in-place) ##
(i::Integral)(args...; params...) =
    integrate(i, evaluate_domain(domain(i), args; params...), args; params...)

# integrate can assume the domain is not a Domain.Functional.
# We just need backend-specific conversions
integrate(i::Integral, domain, args; params...) =
    i.solver(convert_integrand(i, domain, args; params...), convert_domain(domain, i.solver), i.result)

# Any Domain Sum is handled by summing over the domains ##
# non-mutating version
function integrate(i::Integral{Missing}, domain::Domain.Sum, args; params...)
    result = sum(ungroup(domain)) do subdomain
        subdomain´ = evaluate_domain(subdomain, args; params...)
        integrate(i, convert_domain(subdomain´, i.solver), args; params...)
    end
    return result
end

# mutating version
function integrate(i::Integral, domain::Domain.Sum, args; params...)
    resultsum = zero(result(i))
    foreach(ungroup(domain)) do subdomain
        subdomain´ = evaluate_domain(subdomain, args; params...)
        resultsum .+= integrate(i, convert_domain(subdomain´, i.solver), args; params...)
    end
    return resultsum
end

evaluate_domain(d::Domain.Functional, args; params...) = d(args...; params...)
evaluate_domain(d::AbstractDomain, args; params...) = d

## Extension fallbacks and generics ##
#   convert_domain, convert_integrand and solver(f, domain, result)

# generic domain conversions (extensions must opt-in to these explicitly)
convert_domain_generic(d::Domain.Segment{<:Real,<:Real}) =
    (d.x1, d.x2)
convert_domain_generic(d::Domain.Box{N,<:NTuple{N,Real},<:NTuple{N,Real}}) where {N} =
    (d.mins, d.maxs)

# Deal with Infinity and complex-to-real conversions - see changeofvariables.jl
convert_domain_generic(d::Domain.Segment) = convert_domain_generic(transform_domain(d))
convert_domain_generic(d::Domain.Box) = convert_domain_generic(transform_domain(d))

# generic integrand conversions (extensions must opt-in to these explicitly)
function convert_integrand_generic(i::Integral{Missing}, domain, args; post = identity, params...)
    f = integrand(i)
    f´(t) = post(jacobian(t, domain) * f(change_of_variables(t, domain)..., args...; params...))
    return f´
end

function convert_integrand_generic(i::Integral, domain, args; post = identity, params...)
    f! = integrand(i)
    function f!´(out, t)
        f!(out, change_of_variables(t, domain)..., args...; params...)
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
    throw(ArgumentError("No conversion method for domain $(domainname(d)) defined for this backend $(solvername(s)), or solver backend not loaded."))

convert_integrand(i::Integral, domain, args; params...) =
    throw(ArgumentError("No conversion method for the integrand defined for this backend $(solvername(i)), or solver backend not loaded."))

(s::AbstractBackend)(f, domain, result) =
    error("The integral solver solver for $(nameof(typeof(s))) is not loaded.")

## Show ##

Base.summary(i::Integral) = "Integral"

function Base.show(io::IO, J::Integral)
    i = get(io, :indent, "")
    ioindent = IOContext(io, :indent => i * "  ")
    print(io, summary(J), "\n",
"$i  Mutating   : $(ismutating(J))
$i  Domain     : $(domainname(J))
$i  Solver     : $(solvername(J))
$i  Integrand  : ")
  print(ioindent, integrand(J))
end
