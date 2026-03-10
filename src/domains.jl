abstract type AbstractDomain end

module Domain

using IntegrationInterface: AbstractDomain
import IntegrationInterface as II

struct Segment{T<:Number} <: AbstractDomain
    x1::T
    x2::T
end

struct SegmentGroup{T<:Number} <: AbstractDomain
    segments::Vector{NTuple{2,T}}
end

struct Box{N,T<:Number} <: AbstractDomain
    mins::NTuple{N,T}
    maxs::NTuple{N,T}
end

struct Functional{F} <: AbstractDomain
    f::F
end

II.domainname(d::Segment) = "Segment($(d.x1), $(d.x2))"
II.domainname(d::SegmentGroup) = "SegmentGroup($(d.segments))"
II.domainname(d::Box) = "Box($(d.mins), $(d.maxs))"

## API ##

# constructors #

Segment(x1::Number, x2::Number) = Segment(promote(x1, x2)...)

Segment(s::NTuple{2,Number}...) = SegmentGroup(collect(s))

function Segment(node1::Number, node2::Number, node3::Number, nodes::Number...)
    nodes´ = promote(node1, node2, node3, nodes...)
    return SegmentGroup(collect(zip(Base.front(nodes´), Base.tail(nodes´))))
end

# conversions #
ungroup(ss::SegmentGroup) = (Segment(s) for s in ss.segments)

end
