# We bundle AbstractBackends in a submodule to avoid name clashes with their
# respective solver packages.

abstract type AbstractBackend end

# Support for Domain.Functional conversion for all backends. This dispatches to the
# domain-specific convert_domain method, which is defined by each backend extension
convert_domain(d::Domain.Functional, s::AbstractBackend, args) =
    convert_domain(d(args...), s)
convert_domain(d::AbstractDomain, s::AbstractBackend, _) =
    convert_domain(d, s)
convert_domain(::AbstractDomain, s::AbstractBackend) =
    throw(ArgumentError("No conversion method for domain $(domainname(d)) defined for this backend. Forgot `using` the package for backend $(solvername(s))?"))

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
