### Change of variables and corresponding domain transformations ###


### Infinity ###
# Deal with Infinity and complex domain bounds.
# See https://github.com/pablosanjose/IntegrationInterface.jl/issues/3 for details

## API ##

point(d::Infinity) = d.point
point(x::Number) = x

Base.:+(d::Infinity) = d
Base.:-(d::Infinity) = Infinity(-point(d))

## Tools for convert_integrand in the presence of Infinity and complex domains ##

# fallbacks
change_of_variables(t, ::AbstractDomain) = t

transform_domain(d::AbstractDomain) = d

jacobian(t, d::AbstractDomain) = 1

# Domain.Line
change_of_variables(t, d::Domain.Line{<:Number,<:Infinity}) = d.x1 + delta(d) * t/(1-t)
change_of_variables(t, d::Domain.Line{<:Infinity, <:Number}) = d.x2 - delta(d) * t/(1-t)
change_of_variables(t, d::Domain.Line{<:Infinity, <:Infinity}) =
    0.5*(point(d.x1)+point(d.x2)) + 0.75 * delta(d) * t/(1-t^2)
# complex case
change_of_variables(t, d::Domain.Line{<:ComplexOrReal,<:ComplexOrReal}) = d.x1 + delta(d)*t
change_of_variables(t, ::Domain.Line{<:Real,<:Real}) = t             # no transformation

jacobian(t, d::Domain.Line{<:Number,<:Infinity}) = delta(d)/(1-t)^2
jacobian(t, d::Domain.Line{<:Infinity,<:Number}) = delta(d)/(1-t)^2
jacobian(t, d::Domain.Line{<:Infinity,<:Infinity}) = 0.75*delta(d)*(1+t^2)/(1-t^2)^2
jacobian(t, d::Domain.Line{<:ComplexOrReal,<:ComplexOrReal}) = delta(d)
jacobian(t, ::Domain.Line{<:Real,<:Real}) = 1                        # no transformation

transform_domain(::Domain.Line{<:Number,<:Infinity}) = Domain.Line(0.0, 1.0)
transform_domain(::Domain.Line{<:Infinity,<:Number}) = Domain.Line(0.0, 1.0)
transform_domain(::Domain.Line{<:Infinity,<:Infinity}) = Domain.Line(-1.0, 1.0)
transform_domain(::Domain.Line{<:ComplexOrReal,<:ComplexOrReal}) = Domain.Line(0.0, 1.0)
transform_domain(d::Domain.Line{<:Real,<:Real}) = d                  # no transformation

delta(d::Domain.Line) = point(d.x2) - point(d.x1)

# Domain.Box - this is finicky, can easily produce allocations
change_of_variables(t, d::Domain.Box) = change_of_variables.(t, Domain.to_lines(d))

jacobian(t, d::Domain.Box) = prod(jacobian.(t, Domain.to_lines(d)))

transform_domain(d::Domain.Box) = Domain.to_box(transform_domain.(Domain.to_lines(d)))
