# We bundle AbstractIntegralSolvers in a submodule to avoid name clashes with their
# respective backend packages. We define the IS alias for easy access to IntegrationSolvers.

abstract type AbstractIntegralSolver end

module IntegrationSolvers

using IntegrationInterface: AbstractIntegralSolver

## Solvers ##

# Requires the QuadGK package
struct QuadGK{O<:NamedTuple} <: AbstractIntegralSolver
    opts::O
end

QuadGK(; opts...) = QuadGK(NamedTuple(opts))

# Requires the Cubature package
struct Cubature{O<:NamedTuple} <: AbstractIntegralSolver
	opts::O
end

Cubature(; opts...) = Cubature(NamedTuple(opts))

# Requires the HCubature package
struct HCubature{O<:NamedTuple} <: AbstractIntegralSolver
	opts::O
end

HCubature(; opts...) = HCubature(NamedTuple(opts))

struct Quadrature <: AbstractIntegralSolver
	nodes::Vector{Float64}
    weights::Vector{Float64}
end

Quadrature((nodes, weights)) = Quadrature(nodes, weights)

struct Multi{N,S<:NTuple{N,AbstractIntegralSolver}} <: AbstractIntegralSolver
  solvers::S
end

Multi(solvers::AbstractIntegralSolver...) = Multi(solvers)
Multi(solver::AbstractIntegralSolver, N) = Multi(ntuple(Returns(solver), N))

end # module

const IS = IntegrationSolvers

## Quadrature ##

function (s::IS.Quadrature)(f, domain, ::Missing, args; params...)
    xmin, xmax = domain
    Δx = 0.5 * (xmax - xmin)
    result = sum(zip(s.nodes, s.weights); init = zero(float(Δx))) do (x, w)
        f(xmin + Δx * (x + 1), args...; params...) * w * Δx
    end
    return result
end

function (s::IS.Quadrature)(f!, domain, result, args; params...)
    xmin, xmax = domain
    Δx = 0.5 * (xmax - xmin)
    fill!(result, 0)
    out = similar(result)
    result = foreach(zip(s.nodes, s.weights)) do (x, w)
        result .+= f!(out, xmin + Δx * (x + 1), args...; params...) .* w .* Δx
    end
    return result
end

## Multisolver for nested integrals ##

# add an extra argument with outer point variables
(s::IS.Multi)(f, domains, result, args; params...) = s(f, domains, result, args, (); params...)
# recursive solver call, starting from last. Only innermost is potentially mutating
(s::IS.Multi)(f, domains, result, args, point; params...) =
	last(s)(maybe_evaluate_domain(last(domains), point), missing, args; params...) do z
		Base.front(s)(f, Base.front(domains), result, args, (point..., z); params...)
	end
# Innermost integral
(s::IS.Multi{1})(f, domains, ::Missing, args, point; params...) =
	only(s)((x, as...; ps...) -> f((x, point...), as...; ps...), maybe_evaluate_domain(only(domains), point), missing, args; params...)
(s::IS.Multi{1})(f!, domains, result, args, point; params...) =
	only(s)((out, x, as...; ps...) -> f!(out, (x, point...), as...; ps...), maybe_evaluate_domain(only(domains), point), result, args; params...)

maybe_evaluate_domain(domain, _) = domain
maybe_evaluate_domain(domain::Function, point) = domain(point...)

# required for recursion
Base.front(s::IS.Multi) = IS.Multi(Base.front(s.solvers))
Base.last(s::IS.Multi) = last(s.solvers)
Base.only(s::IS.Multi) = only(s.solvers)
Base.length(s::IS.Multi) = length(s.solvers)

sanitize_domain(domain, s::IS.Multi) = length(domain) == length(s) ? domain :
    throw(ArgumentError("Wrong number of domains for Multi solver, expected $(length(s)), got $(length(domain))"))

solvername(s::IS.Multi) = string("Multi(", join(solvername.(s.solvers), ", "), ")")

## Serialize/Deserialize ##
# Needed to convert arrays to vectors of a given type and back, in-place (non-allocating)

serialize_array(::Type{T}, a::AbstractArray) where {T} = reinterpret(T, serialize_array(a))
serialize_array(a::AbstractArray) = vec(a)

deserialize_array(a::AbstractArray{T}, v::AbstractArray) where {T} =
    deserialize_array(a, reinterpret(T, v))
deserialize_array(a::AbstractArray{T}, v::AbstractArray{T}) where {T} =
    (check_deserializer(a, v); unsafe_deserialize_array(a, v))

# assumes equal eltype T and compatible size
unsafe_deserialize_array(::AbstractArray{T,N}, v::AbstractArray{T,N}) where {T,N} = v
unsafe_deserialize_array(a::AbstractArray, v::AbstractArray) = reshape(v, size(a))

check_deserializer(a, v) =
    size(serialize(a)) == size(v) ||
        argerror("Wrong size of serialized array, expected $(size(serialize(a))), got $(size(v))")

check_deserializer(a, v::AbstractVector) =
    length(serialize(a)) == length(v) ||
        argerror("Wrong length of serialized array, expected $(length(serialize(a))), got $(length(v))")
