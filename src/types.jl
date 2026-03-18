# All types not belonging in submodules are defined up-front to avoid include ordering problems

abstract type AbstractDomain end

abstract type AbstractBackend end

# general callable object representing an integral over a domain using a backend
struct Integral{R,S<:AbstractBackend,D<:AbstractDomain,F}
    integrand::F	# usually a function, but can be other kind of object
    result::R		# will be `nothing` if not in-place
    domain::D		# Tuple of AbstractDomains or a single AbstractDomain
    backend::S		# object representing the integration  backend(s)
end

# represents an infinite ray passing through a point (direction fixed by another point´)
struct Infinity{T}
    point::T
    function Infinity{T}(point) where {T}
        isinf(point) && throw(ArgumentError("Ininity(Inf) not allowed"))
        return new(point)
    end
end

Infinity(point::T) where {T} = Infinity{T}(point)

const NumberOrInfinity = Union{Number,Infinity}
