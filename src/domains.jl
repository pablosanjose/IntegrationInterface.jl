module Domain

using IntegrationInterface: AbstractDomain, NumberOrInfinity
import IntegrationInterface as II

struct Box{N,P1<:NTuple{N,NumberOrInfinity},P2<:NTuple{N,NumberOrInfinity}} <: AbstractDomain
    mins::P1
    maxs::P2
    function Box{N,P1,P2}(mins, maxs) where {N,P1<:NTuple{N,NumberOrInfinity},P2<:NTuple{N,NumberOrInfinity}}
        error_if_degenerate.(mins, maxs)
        mins´, maxs´ = sanitize_infs(mins, maxs)
        return new(mins´, maxs´)
    end
end

const Box1D{T1,T2} = Box{1,Tuple{T1},Tuple{T2}}

Box(mins::P1, maxs::P2) where {N,P1<:NTuple{N,NumberOrInfinity},P2<:NTuple{N,NumberOrInfinity}} =
    Box{N,P1,P2}(mins, maxs)
struct Functional{D<:AbstractDomain,F} <: AbstractDomain
    type::Type{D}
    f::F
end

struct Sum{T} <: AbstractDomain
    subdomains::T   # subdomains should be an iterator over AbstractDomains
end

## Sanitization ##
error_if_degenerate(::Number, ::Number) = nothing
error_if_degenerate(x1::NumberOrInfinity, x2::NumberOrInfinity) = II.point(x1) == II.point(x2) &&
    throw(ArgumentError("Got a domain corresponding to an unbounded ray with an ill-defined direction. "))

# Should not mix Infinity with Inf. Also, Complex infs are not allowed
function sanitize_infs(x1::Tuple, x2::Tuple)
    s = sanitize_infs.(x1, x2)
    return first.(s), last.(s)
end

sanitize_infs(x1::NumberOrInfinity, x2::NumberOrInfinity) =
    _sanitize_infs(sanitize_complex_inf(x1), sanitize_complex_inf(x2))
_sanitize_infs(x1::NumberOrInfinity, x2::NumberOrInfinity) = (x1, x2)
_sanitize_infs(x1::II.Infinity, x2::Number) = reverse(_sanitize_infs(x2, x1))
_sanitize_infs(x1::Number, x2::II.Infinity) =
    isinf(x1) ? throw(ArgumentError("Mixing `Inf` with `Infinity` is ambiguous.")) : (x1, x2)
sanitize_complex_inf(x::Complex) =
    isinf(x) ? throw(ArgumentError("Complex Inf not allowed in domains")) : x
sanitize_complex_inf(x) = x

## Show ##
II.domainname(d::Sum) = string("Sum(", join(II.domainname.(d.subdomains), ", "), ")")
II.domainname(d::Box{N}) where {N} = "Box{$N}(($(short_show(d.mins...))), ($(short_show(d.maxs...)))))"
II.domainname(d::Box{1}) = "Box{1}($(short_show(only(d.mins))), $(short_show(only(d.maxs))))"
II.domainname(d::Functional) = "Functional{$(nameof(d.type))}"

short_show(d::II.Infinity) = "Infinity($(II.point(d)))"
short_show(d::Number) = "$d"
short_show(xs...) = join(short_show.(xs), ", ")

## API ##

# 1D box
Box1D(a::NumberOrInfinity, b::NumberOrInfinity) = Box((a,), (b,))

# Sum of consecutive 1D boxes
function Box1D(node1::NumberOrInfinity, node2::NumberOrInfinity, node3::NumberOrInfinity, nodes::NumberOrInfinity...)
    nodes´ = (node1, node2, node3, nodes...)
    return Sum(Box1D.(Base.front(nodes´), Base.tail(nodes´)))
end

# As above, but with an AbstractVector
function Box1D(nodes::AbstractVector)
    return Sum(Box1D(nodes[i], nodes[i+1]) for i in eachindex(nodes)[1:end-1])
end

# Functional domain
Box{N}(f::Function) where {N} = Functional(Box{N}, f)

(f::Functional{Box{N}})(args...; kw...) where {N} = Box(f.f(args...; kw...)...)::Box{N}

# Sum of domains
Sum(xs::AbstractDomain...) = Sum(xs)

# Domain.interval helper

interval(a::NumberOrInfinity, b::NumberOrInfinity) = Box1D(a, b)
interval(is::NTuple{2,NumberOrInfinity}...) = Box(first.(is), last.(is))

# accessors #

Base.first(d::Box1D) = only(d.mins)
Base.last(d::Box1D) = only(d.maxs)

Base.first(d::Box) = d.mins
Base.last(d::Box) = d.maxs

# call #

(f::Functional{D})(args...; kw...) where {D} = D(f.f(args...; kw...)...)

# ungroup domain sums

II.ungroup(ss::Sum) = ss.subdomains

# conversion

to_1D_boxes(d::Box{N}) where {N} = ntuple(i -> Box1D(d.mins[i], d.maxs[i]), Val(N))

to_box(d::NTuple{N,Box1D}) where {N} = Box(first.(d), last.(d))

end
