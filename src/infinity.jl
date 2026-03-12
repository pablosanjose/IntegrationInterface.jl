## Infinity ##
# Deal with Infinity domain bounds.
# See https://github.com/pablosanjose/IntegrationInterface.jl/issues/3 for details

## API ##

point(d::Infinity) = d.point
point(x::Number) = x

Base.:+(d::Infinity) = d
Base.:-(d::Infinity) = Infinity(-point(d))

## Tools for convert_integrand in the presence of Infinity ##

# fallbacks
change_of_variables(t, ::AbstractDomain) = t

transform_domain(d::AbstractDomain) = d

jacobian(t, d::AbstractDomain) = 1

# Domain.Segment
change_of_variables(t, d::Domain.Segment{<:Number,<:Infinity}) = d.x1 + delta(d) * t/(1-t)
change_of_variables(t, d::Domain.Segment{<:Infinity, <:Number}) = d.x2 - delta(d) * t/(1-t)
change_of_variables(t, d::Domain.Segment{<:Infinity, <:Infinity}) =
    0.5*(point(d.x1)+point(d.x2)) + 0.75 * delta(d) * t/(1-t^2)

jacobian(t, d::Domain.Segment{<:Number,<:Infinity}) = delta(d)/(1-t)^2
jacobian(t, d::Domain.Segment{<:Infinity,<:Number}) = delta(d)/(1-t)^2
jacobian(t, d::Domain.Segment{<:Infinity,<:Infinity}) = 0.75*delta(d)*(1+t^2)/(1-t^2)^2

transform_domain(::Domain.Segment{<:Number,<:Infinity}) = Domain.Segment(0.0, 1.0)
transform_domain(::Domain.Segment{<:Infinity,<:Number}) = Domain.Segment(0.0, 1.0)
transform_domain(::Domain.Segment{<:Infinity,<:Infinity}) = Domain.Segment(-1.0, 1.0)

delta(d::Domain.Segment) = point(d.x2) - point(d.x1)

# Domain.Box - this is finicky, can easily produce allocations
change_of_variables(t, d::Domain.Box) = change_of_variables.(t, Domain.to_segments(d))

jacobian(t, d::Domain.Box) = prod(jacobian.(t, Domain.to_segments(d)))

transform_domain(d::Domain.Box) = Domain.to_box(transform_domain.(Domain.to_segments(d)))
