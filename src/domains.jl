module Domain

using IntegrationInterface: AbstractDomain, Infinity
import IntegrationInterface as II

const MaybeInfinity{T} = Union{T,Infinity{T}}
const NumberOrInfinity = MaybeInfinity{<:Number}            # a number or a scalar ray
const NumberOrNTuple = Union{Number,NTuple{<:Any,Number}}   # a number or a point of numbers


struct Sum{T} <: AbstractDomain
    subdomains::T   # subdomains should be an iterator over AbstractDomains
end

struct Functional{D<:AbstractDomain,F} <: AbstractDomain
    type::Type{D}
    f::F
end

struct Box{N,T<:Number,P1<:NTuple{N,MaybeInfinity{T}},P2<:NTuple{N,MaybeInfinity{T}}} <: AbstractDomain
    mins::P1
    maxs::P2
    function Box{N,T,P1,P2}(mins, maxs) where {N,T<:Number,P1<:NTuple{N,MaybeInfinity{T}},P2<:NTuple{N,MaybeInfinity{T}}}
        error_if_degenerate.(mins, maxs)
        error_if_mixed_infs.(mins, maxs)
        return new(mins, maxs)
    end
end

# type promotion
function Box(mins::NTuple{N,Any}, maxs::NTuple{N,Any}) where {N}
    T = promote_type_infs(mins..., maxs...)
    mins´, maxs´ = promote_inf.(T, mins), promote_inf.(T, maxs)
    P1, P2 = typeof(mins´), typeof(maxs´)
    return Box{N,T,P1,P2}(mins´, maxs´)
end

## Sanitization ##

# We cannot have integers as Domain types
function promote_type_infs(xs...)
    T = promote_type(typeof.(II.point.(xs))...)
    isconcretetype(T) || throw(ArgumentError("Couldn't find a concrete promotion type for domain"))
    return float(T)
struct Simplex{N,P<:SimplexVertices{N}} <: AbstractDomain
    vertices::P  # collection of N+1 vertices in N-dimensions
    function Simplex{N,P}(vertices) where {N,P<:SimplexVertices{N}}
        for j in 1:N+1, i in j+1:N+1
            error_if_degenerate(vertices[i], vertices[j])
            error_if_mixed_infs(vertices[i], vertices[j])
        end
        return new(vertices)
    end
end

Simplex(vertices::P) where {N,P<:SimplexVertices{N}} = Simplex{N,P}(vertices)

const RealFiniteSimplex{N} = Simplex{N,<:SimplexVerticesRealFinite{N}}

## Sanitization ##

promote_inf(T::Type, x::Number) = convert(T, x)
promote_inf(T::Type, x::Infinity) = Infinity(convert(T, II.point(x)))

error_if_degenerate(::NumberOrNTuple, ::NumberOrNTuple) = nothing
error_if_degenerate(x1, x2) = II.point(x1) == II.point(x2) &&
    throw(ArgumentError("Got a domain corresponding to an unbounded ray with an ill-defined direction. "))

# Should not mix Infinity with Inf. Also, Complex infs are not allowed
error_if_mixed_infs(x1, x2) =
    _error_if_mixed_infs(_sanitize_complex_inf(x1), _sanitize_complex_inf(x2))

_error_if_mixed_infs(::Infinity, ::Infinity) = nothing
_error_if_mixed_infs(_, _) = nothing
_error_if_mixed_infs(::Infinity, x) = isinf_any(x) && error_mixed_infs()
_error_if_mixed_infs(x, ::Infinity) = isinf_any(x) && error_mixed_infs()
_sanitize_complex_inf(x::Complex) =
    isinf(x) ? throw(ArgumentError("Complex Inf not allowed in domains")) : x

_sanitize_complex_inf(p::Tuple) = _sanitize_complex_inf.(p)
_sanitize_complex_inf(x) = x

isinf_any(x::Number) = isinf(x)
isinf_any(x::Tuple) = any(isinf, x)

error_mixed_infs() = throw(ArgumentError("Mixing `Inf` with `Infinity` is ambiguous."))

## Show ##
II.domainname(d::Sum) = string("Sum(", join(II.domainname.(d.subdomains), ", "), ")")
II.domainname(d::Box{N}) where {N} = "Box{$N}(($(short_show(d.mins...))), ($(short_show(d.maxs...)))))"
II.domainname(d::Box{1}) = "Box{1}($(short_show(only(d.mins))), $(short_show(only(d.maxs))))"
II.domainname(d::Functional) = "Functional{$(II.domainname(d.type))}"
II.domainname(::Type{<:Box{N}}) where {N} = "Box{$N}"
II.domainname(d::Type) = nameof(d)

short_show(d::Infinity) = "Infinity($(II.point(d)))"
short_show(d::Number) = "$d"
short_show(xs...) = join(short_show.(xs), ", ")

## API ##

# 1D box
Box{1}(a::NumberOrInfinity, b::NumberOrInfinity) = Box((a,), (b,))
Box(a::NumberOrInfinity, b::NumberOrInfinity) = Box((a,), (b,))

# Sum of consecutive 1D boxes
function Box{1}(node1::NumberOrInfinity, node2::NumberOrInfinity, node3::NumberOrInfinity, nodes::NumberOrInfinity...)
    nodes´ = (node1, node2, node3, nodes...)
    return Sum(Box{1}.(Base.front(nodes´), Base.tail(nodes´)))
end

# As above, but with an AbstractVector
function Box{1}(nodes::AbstractVector)
    return Sum(Box{1}(nodes[i], nodes[i+1]) for i in eachindex(nodes)[1:end-1])
end

# Simplex

Simplex(xs::NTuple...) = Simplex(xs)

# unit simplex

Base.oneunit(d::Simplex) = oneunit(typeof(d))

function Base.oneunit(::Type{<:Simplex{N}}) where {N}
    z = ntuple(Returns(0), Val(N))
    u = ntuple(j ->ntuple(i->ifelse(i==j, 1, 0), Val(N)), Val(N))
    return Simplex(z, u...)
end

# Functional domain
Box{N}(f::Function) where {N} = Functional(Box{N}, f)
Simplex{N}(f::Function) where {N} = Functional(Simplex{N}, f)

(::Type{D})(::Function) where {D<:AbstractDomain} =
    throw(ArgumentError("`$(nameof(D))(::Function)` not supported, use `$(nameof(D)){N}(::Function)` instead."))

# call
(f::Functional{Box{N}})(args...; kw...) where {N} = Box(f.f(args...; kw...)...)::Box{N}
(f::Functional{Simplex{N}})(args...; kw...) where {N} = Simplex(f.f(args...; kw...)...)::Simplex{N}

# Sum of domains
Sum(xs::AbstractDomain...) = Sum(xs)

# Domain.interval helper

interval(a::NumberOrInfinity, b::NumberOrInfinity) = Box{1}(a, b)
interval(is::NTuple{2,NumberOrInfinity}...) = Box(first.(is), last.(is))

# accessors #

Base.first(d::Box{1}) = only(d.mins)
Base.last(d::Box{1}) = only(d.maxs)

Base.first(d::Box) = d.mins
Base.last(d::Box) = d.maxs

vertices(s::Simplex) = s.vertices

# ungroup domain sums

II.ungroup(ss::Sum) = ss.subdomains

# conversion

to_1D_boxes(d::Box{N}) where {N} = ntuple(i -> Box{1}(d.mins[i], d.maxs[i]), Val(N))

to_box(d::NTuple{N,Box{1}}) where {N} = Box(first.(d), last.(d))

end
