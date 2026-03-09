# We bundle AbstractBackends in a submodule to avoid name clashes with their
# respective solver packages. We define the IS alias for easy access to Backend.

abstract type AbstractBackend end

module Backend

using IntegrationInterface: AbstractBackend

## Solvers ##

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
