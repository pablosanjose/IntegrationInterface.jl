# Integration solver backends are provided by external libraries and their extensions
# Extensions must provide:
#   - `convert_domain(domain, backend)`: a domain understood by backend, can fall back to
#     `convert_domain_generic(domain)`
#   - `convert_integrand(i::Integral, domain, args; params...)`: a function understood by
#      the backend, can fall back to `convert_domain_generic(domain)`
#   - `(::Backend)(integrand::Function, domain, result)`: call to actual implementation,
#     once `integrand` and `domain` have been converted.
#
# We bundle all AbstractBackends in a Backend submodule to avoid name clashes with their
# respective solver packages.


## Collection of backend solver types ##
module Backend

using IntegrationInterface: AbstractBackend

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

struct Quadrature <: AbstractBackend
	nodes::Vector{Float64}
    weights::Vector{Float64}
end

Quadrature((nodes, weights)) = Quadrature(nodes, weights)

end # module

const IS = Backend
