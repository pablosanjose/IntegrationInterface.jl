# Integration backends are provided by external libraries and their extensions
# Extensions must provide:
#   - `convert_domain(domain, backend)`: a domain understood by backend, can fall back to
#     `convert_domain_generic(domain)`
#   - `convert_integrand(i::Integral, backend, domain, args; kw...)`: a function understood by
#      the backend, can fall back to `convert_domain_generic(domain)`
#   - `(::Backend)(integrand::Function, domain, result)`: call to actual implementation,
#     once `integrand` and `domain` have been converted.
#
# We bundle all AbstractBackends in a Backend submodule to avoid name clashes with their
# respective backend packages.


## Collection of backend types ##
module Backend

using IntegrationInterface: AbstractBackend, AbstractEvaluatedDomain, Domain, domainname

# singleton default backend
struct Default <: AbstractBackend end

# Requires the QuadGK package
struct QuadGK{O<:NamedTuple} <: AbstractBackend
    opts::O
end

QuadGK(; opts...) = QuadGK(NamedTuple(opts))

# Requires the Cubature package
struct Cubature{O<:NamedTuple} <: AbstractBackend
	opts::O
end

Cubature(; opts...) = Cubature(NamedTuple(opts))

# Requires the HCubature package
struct HCubature{O<:NamedTuple} <: AbstractBackend
	opts::O
end

HCubature(; opts...) = HCubature(NamedTuple(opts))

struct Quadrature{T} <: AbstractBackend
	nodes::Vector{T}
    weights::Vector{T}
end

Quadrature((nodes, weights)) = Quadrature(nodes, weights)

struct HAdaptiveIntegration{O<:NamedTuple} <: AbstractBackend
	opts::O
end

HAdaptiveIntegration(; opts...) = HAdaptiveIntegration(NamedTuple(opts))

# defaults for each domain type
default(::Domain.Box{1}) = QuadGK()
default(::Domain.Box) = HCubature()
default(::Domain.Simplex) = HAdaptiveIntegration()
default(_) = Default()

# Can be overridden for user-defined f types and domains
resolve(b::AbstractBackend, _) = b
resolve(::Backend.Default, d::Domain.AbstractEvaluatedDomain) = default(d)
resolve(::Backend.Default, d) =
    throw(ArgumentError("No default backend exists for domain $(domainname(d)), please specify one explicitly."))

end # module
