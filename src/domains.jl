module Domain

using IntegrationInterface: AbstractDomain, Infinity
import IntegrationInterface as II

const MaybeInfinity{T} = Union{T,Infinity{<:T}}
const NumberOrInfinity{T<:Number} = MaybeInfinity{T}            # a number or a scalar ray
const NumberOrNTuple = Union{Number,NTuple{<:Any,Number}}   # a number or a point of numbers

const TypedVertices{N,T} = Tuple{T,Vararg{T,N}}
const SimplexVertices{N,T} = TypedVertices{N,MaybeInfinity{NTuple{N,T}}}

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

function Box(mins::NTuple{N,Any}, maxs::NTuple{N,Any}) where {N}
    T = promote_type_infs(mins..., maxs...)
    mins´, maxs´ = promote_inf.(T, mins), promote_inf.(T, maxs)
    P1, P2 = typeof(mins´), typeof(maxs´)
    return Box{N,T,P1,P2}(mins´, maxs´)
end

const FiniteBox{N,T} = Box{N,T,NTuple{N,T},NTuple{N,T}}
const FiniteRealBox{N} = FiniteBox{N,<:Real}

struct Simplex{N,T<:Number,P<:SimplexVertices{N,T},D} <: AbstractDomain
    vertices::P          # collection of N+1 vertices in N-dimensions
    basisdata::D         # will be nothing (not needed) if the Simplex is real and finite
    function Simplex{N,T,P,D}(vertices, basisdata) where {N,T<:Number,P<:SimplexVertices{N,T},D}
        onpairs(error_if_degenerate, vertices...)
        onpairs(error_if_mixed_infs, vertices...)
        error_if_inf_Simplex(vertices...)
        return new(vertices, basisdata)
    end
end

function Simplex(vertices::MaybeInfinity{NTuple{N,Any}}...) where {N}
    length(vertices) == N+1 ||
        throw(ArgumentError("Expected $(N+1) vertices for an $N-simplex, got $(length(vertices))."))
    T = promote_type_infs(vertices...)
    vertices´ = promote_inf.(T, vertices)
    basisdata´ = basisdata(vertices´)   # see section below on unbounded simplices
    P, D = typeof(vertices´), typeof(basisdata´)
    return Simplex{N,T,P,D}(vertices´, basisdata´)
end

basisdata(d::Simplex) = d.basisdata

const FiniteSimplex{N,T} = Simplex{N,T,<:TypedVertices{N,NTuple{N,T}}}
const FiniteRealSimplex{N} = FiniteSimplex{N,<:Real}
const InfiniteSimplex{N,T} = Simplex{N,T,<:TypedVertices{N,Infinity{NTuple{N,T}}}}

## Sanitization ##

# We avoid Integers as Domain types
function promote_type_infs(xs::NumberOrInfinity...)
    T = promote_type(typeof.(II.point.(xs))...)
    isconcretetype(T) || throw(ArgumentError("Couldn't find a concrete promotion type for domain, got $T"))
    return float(T)
end

promote_type_infs(xs::MaybeInfinity{NTuple{N,Any}}...) where {N} =
    promote_type_infs(tuplejoin(II.point.(xs)...)...)

promote_inf(T::Type, x::Number) = convert(T, x)
promote_inf(T::Type, x::Tuple) = promote_inf.(T, x)
promote_inf(T::Type, x::Infinity) = Infinity(promote_inf(T, II.point(x)))

# detect v=Infinity(x) such that x === point(v) == point(v´) for another vertex v´
error_if_degenerate(::NumberOrNTuple, ::NumberOrNTuple) = nothing
error_if_degenerate(x1, x2) = II.point(x1) == II.point(x2) &&
    throw(ArgumentError("Got a domain corresponding to an unbounded ray with an ill-defined direction. "))

# Should not mix Infinity with Inf. Also, Complex Infs are not allowed
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

error_mixed_infs() = throw(ArgumentError("Mixing `Inf` with `Infinity` is not allowed."))

error_if_inf_Simplex(::Infinity...) = throw(ArgumentError("Simplex must have at least one finite vertex."))
error_if_inf_Simplex(_...) = nothing

## Util ##

# joins several tuples into one

tuplejoin() = ()
tuplejoin(x) = x
tuplejoin(x, y) = (x..., y...)
tuplejoin(x, y, z...) = (x..., tuplejoin(y, z...)...)

# operate f on all pairs of arguments, with no regard for order

onpairs(f::Function, a, b) = f(a, b)

function onpairs(f::Function, a, b, c, d...)
    foreach(x -> f(a, x), (b, c, d...))
    return onpairs(f, b, c, d...)
end

## Show ##
II.domainname(d::Box{1,T}) where {T}= "Box{1,$T}($(short_show(only(d.mins))), $(short_show(only(d.maxs))))"
II.domainname(d::Box{N,T}) where {N,T} = "Box{$N, $T}(($(short_show(d.mins...))), ($(short_show(d.maxs...)))))"
II.domainname(d::Simplex{N,T}) where {N,T} = "Simplex{$N, $T}($(short_show(d.vertices...))))"
II.domainname(::Type{<:Box{N}}) where {N} = "Box{$N}"
II.domainname(::Type{<:Simplex{N}}) where {N} = "Simplex{$N}"
II.domainname(d::Functional) = "Functional{$(II.domainname(d.type))}"
II.domainname(d::Sum) = string("Sum(", join(II.domainname.(d.subdomains), ", "), ")")
II.domainname(d::Type) = nameof(d)

short_show(d::Infinity) = "Infinity($(II.point(d)))"
short_show(d::Number) = "$d"
short_show(xs::Tuple) = join(short_show.(xs), ", ")
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

## Unbounded simplices ##

# Basis data: for non real-finite-simplices,
#   - order vertices with one finte, then infinite, then rest. Best for change of variables.
# and then collect:
#   - origin: point(vertices´[1])
#   - basis: vectors between origin and the rest.
#   - boxes: normalized Box{1} for each of the simplex axes
#   - signature: sign of the basis permutation

basisdata(::TypedVertices{N,NTuple{N,Real}}) where {N} = nothing

function basisdata(vertices::SimplexVertices{N,T}) where {N,T}
    signature, vertices´ = infinity_to_front_but_one(vertices...)
    origin, vs... = II.point.(vertices´)
    basis = ntuple(i -> vs[i] .- origin, Val(N))
    boxes = ntuple(i -> Box{1}(zero(T), maybe_Infinity(T(1), vertices´[1+i])), Val(N))
    return origin, basis, boxes, signature
end

# returns reorder and signature of permutation
function infinity_to_front_but_one(vs...)
    ref = Ref(1)
    vs´ = infinity_to_front(ref, vs...)
    vs´´ = (last(vs´), Base.front(vs´)...)
    signature = ref[] * ifelse(isodd(length(vs´)-1), -1, 1)
    return signature, vs´´
end

function infinity_to_front(ref, v, vs...)
    isodd(length(vs)) && (ref[] *= -1)
    return (infinity_to_front(ref, vs...)..., v)
end
infinity_to_front(ref, v::Infinity, vs...) = (v, infinity_to_front(ref, vs...)...)
infinity_to_front(ref) = ()

# orthogonal: build an "orthogonal domain" corresponding to interval (a,b) along all axes,
# where a, b can be numbers or Infinity

orthogonal(::D, args...) where {D<:AbstractDomain} = orthogonal(D, args...)

function orthogonal(::Type{<:Simplex{N}}, a::NumberOrInfinity{T}, b::NumberOrInfinity{T}) where {N,T}
    Ta, Tb = T(II.point(a)), T(II.point(b))
    z = maybe_Infinity(ntuple(Returns(Ta), Val(N)), a)
    u = ntuple(Val(N)) do j
            p = ntuple(Val(N)) do i
                ifelse(i==j, Tb, Ta)
            end
            maybe_Infinity(p, b)
        end
    return Simplex(z, u...)
end

orthogonal(::Type{<:Box{N}}, a::NumberOrInfinity{T}, b::NumberOrInfinity{T}) where {N,T} =
    Box(ntuple(Returns(cast(T, a)), Val(N)), ntuple(Returns(cast(T, b)), Val(N)))

cast(T::Type, a::Infinity) = Infinity(cast(T, point(a)))
cast(T::Type, a::Number) = convert(T, a)

maybe_Infinity(x, ::Infinity) = Infinity(x)
maybe_Infinity(x, _) = x


end # module
