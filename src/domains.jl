module Domain

using IntegrationInterface: AbstractDomain, NumberOrInfinity
import IntegrationInterface as II

struct Line{T1<:NumberOrInfinity,T2<:NumberOrInfinity} <: AbstractDomain
    x1::T1
    x2::T2
    function Line{T1,T2}(x1, x2) where {T1<:NumberOrInfinity,T2<:NumberOrInfinity}
        error_if_degenerate(x1, x2)
        x1´, x2´ = sanitize_infs(x1, x2)
        return new(x1´, x2´)
    end
end

Line(x1::T1, x2::T2) where {T1<:NumberOrInfinity,T2<:NumberOrInfinity} =
    Line{T1,T2}(x1, x2)

struct Box{N,T1<:NTuple{N,NumberOrInfinity},T2<:NTuple{N,NumberOrInfinity}} <: AbstractDomain
    mins::T1
    maxs::T2
    function Box{N,T1,T2}(mins, maxs) where {N,T1<:NTuple{N,NumberOrInfinity},T2<:NTuple{N,NumberOrInfinity}}
        error_if_degenerate.(mins, maxs)
        mins´, maxs´ = sanitize_infs(mins, maxs)
        return new(mins´, maxs´)
    end
end

Box(x1::T1, x2::T2) where {N,T1<:NTuple{N,NumberOrInfinity},T2<:NTuple{N,NumberOrInfinity}} =
    Box{N,T1,T2}(x1, x2)

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
    throw(ArgumentError("Got a Domain.Line corresponding to an unbounded ray with an ill-defined direction. "))

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
II.domainname(d::Line) = "Line($(short_show(d.x1, d.x2)))"
II.domainname(d::Sum) = string("Sum(", join(II.domainname.(d.subdomains), ", "), ")")
II.domainname(d::Box) = "Box(($(short_show(d.mins...))), ($(short_show(d.maxs...)))))"
II.domainname(d::Functional) = "Functional{$(nameof(d.type))}"

short_show(d::II.Infinity) = "Infinity($(II.point(d)))"
short_show(d::Number) = "$d"
short_show(xs...) = join(short_show.(xs), ", ")

## API ##

Line(s::NTuple{2,NumberOrInfinity}...) = Sum(Line.(s))

function Line(node1::NumberOrInfinity, node2::NumberOrInfinity, node3::NumberOrInfinity, nodes::NumberOrInfinity...)
    nodes´ = (node1, node2, node3, nodes...)
    return Sum(Line.(Base.front(nodes´), Base.tail(nodes´)))
end

function Line(nodes::AbstractVector)
    return Sum(Line(nodes[i], nodes[i+1]) for i in eachindex(nodes)[1:end-1])
end

(::Type{D})(f::Function) where {D<:AbstractDomain} = Functional(D, f)

Sum(xs::AbstractDomain...) = Sum(xs)

# accessors #

Base.first(d::Line) = d.x1
Base.last(d::Line) = d.x2

Base.first(d::Box) = d.mins
Base.last(d::Box) = d.maxs

# call #

(f::Functional{D})(args...; kw...) where {D} = D(f.f(args...; kw...)...)

# ungroup domain sums

II.ungroup(ss::Sum) = ss.subdomains

# conversion

to_lines(d::Box{N}) where {N} = ntuple(i -> Line(d.mins[i], d.maxs[i]), Val(N))

to_box(d::NTuple{N,Line}) where {N} = Box(first.(d), last.(d))

end
