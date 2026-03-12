module Domain

using IntegrationInterface: AbstractDomain, NumberOrInfinity
import IntegrationInterface as II

struct Segment{T1<:NumberOrInfinity,T2<:NumberOrInfinity} <: AbstractDomain
    x1::T1
    x2::T2
    function Segment{T1,T2}(x1, x2) where {T1<:NumberOrInfinity,T2<:NumberOrInfinity}
        error_if_degenerate(x1, x2)
        return new(x1, x2)
    end
end

Segment(x1::T1, x2::T2) where {T1<:NumberOrInfinity,T2<:NumberOrInfinity} =
    Segment{T1,T2}(x1, x2)

struct Box{N,T1<:NTuple{N,NumberOrInfinity},T2<:NTuple{N,NumberOrInfinity}} <: AbstractDomain
    mins::T1
    maxs::T2
    function Box{N,T1,T2}(mins, maxs) where {N,T1<:NTuple{N,NumberOrInfinity},T2<:NTuple{N,NumberOrInfinity}}
        error_if_degenerate.(mins, maxs)
        return new(mins, maxs)
    end
end

Box(x1::T1, x2::T2) where {N,T1<:NTuple{N,NumberOrInfinity},T2<:NTuple{N,NumberOrInfinity}} =
    Box{N,T1,T2}(x1, x2)

struct Functional{D<:AbstractDomain,F} <: AbstractDomain
    type::Type{D}
    f::F
end

struct Sum{T<:NTuple{<:Any,AbstractDomain}} <: AbstractDomain
    subdomains::T
end

error_if_degenerate(::Number, ::Number) = nothing
error_if_degenerate(x1::NumberOrInfinity, x2::NumberOrInfinity) = II.point(x1) == II.point(x2) &&
    throw(ArgumentError("Got a Domain.Segment corresponding to an unbounded ray with an ill-defined direction. "))

II.domainname(d::Segment) = "Segment($(short_show(d.x1, d.x2)))"
II.domainname(d::Sum) = string("Sum(", join(II.domainname.(d.subdomains), ", "), ")")
II.domainname(d::Box) = "Box(($(short_show(d.mins...))), ($(short_show(d.maxs...)))))"
II.domainname(d::Functional) = "Functional{$(nameof(d.type))}"

short_show(d::II.Infinity) = "Infinity($(II.point(d)))"
short_show(d::Number) = "$d"
short_show(xs...) = join(short_show.(xs), ", ")

## API ##

# constructors #

Segment(s::NTuple{2,NumberOrInfinity}...) = Sum(Segment.(s))

function Segment(node1::NumberOrInfinity, node2::NumberOrInfinity, node3::NumberOrInfinity, nodes::NumberOrInfinity...)
    nodes´ = (node1, node2, node3, nodes...)
    return Sum(Segment.(Base.front(nodes´), Base.tail(nodes´)))
end

(::Type{D})(f::Function) where {D<:AbstractDomain} = Functional(D, f)

Sum(xs::AbstractDomain...) = Sum(xs)

# Other

Base.first(d::Segment) = d.x1
Base.last(d::Segment) = d.x2

# call #

(f::Functional{D})(args...; params...) where {D} = D(f.f(args...; params...)...)

# ungroup domain sums

II.ungroup(ss::Sum) = ss.subdomains

# conversion

to_segments(d::Box{N}) where {N} = ntuple(i -> Segment(d.mins[i], d.maxs[i]), Val(N))

to_box(d::NTuple{N,Segment}) where {N} = Box(first.(d), last.(d))

end
