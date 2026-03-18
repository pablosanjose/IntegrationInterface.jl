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

# Domain.Box1D
change_of_variables(t, d::Domain.Box1D{<:Number,<:Infinity}) = first(d) + delta(d) * t/(1-t)
change_of_variables(t, d::Domain.Box1D{<:Infinity, <:Number}) = last(d) - delta(d) * t/(1-t)
change_of_variables(t, d::Domain.Box1D{<:Infinity, <:Infinity}) =
    0.5*(point(first(d))+point(last(d))) + 0.75 * delta(d) * t/(1-t^2)
# Anything number that is not real (e.g. complex, unitful) is converted for compatibility
change_of_variables(t, d::Domain.Box1D{<:Number,<:Number}) = first(d) + delta(d)*t
change_of_variables(t, ::Domain.Box1D{<:Real,<:Real}) = t             # no transformation

jacobian(t, d::Domain.Box1D{<:Number,<:Infinity}) = delta(d)/(1-t)^2
jacobian(t, d::Domain.Box1D{<:Infinity,<:Number}) = delta(d)/(1-t)^2
jacobian(t, d::Domain.Box1D{<:Infinity,<:Infinity}) = 0.75*delta(d)*(1+t^2)/(1-t^2)^2
jacobian(t, d::Domain.Box1D{<:Number,<:Number}) = delta(d)
jacobian(t, ::Domain.Box1D{<:Real,<:Real}) = 1                        # no transformation

transform_domain(::Domain.Box1D{<:Number,<:Infinity}) = Domain.Box{1}(0.0, 1.0)
transform_domain(::Domain.Box1D{<:Infinity,<:Number}) = Domain.Box{1}(0.0, 1.0)
transform_domain(::Domain.Box1D{<:Infinity,<:Infinity}) = Domain.Box{1}(-1.0, 1.0)
transform_domain(::Domain.Box1D{<:Number,<:Number}) = Domain.Box{1}(0.0, 1.0)
transform_domain(d::Domain.Box1D{<:Real,<:Real}) = d                  # no transformation

delta(d::Domain.Box1D) = point(last(d)) - point(first(d))

# Domain.Box - this is finicky, can easily produce allocations
change_of_variables(t, d::Domain.Box) = change_of_variables.(t, Domain.to_1D_boxes(d))

jacobian(t, d::Domain.Box) = prod(jacobian.(t, Domain.to_1D_boxes(d)))

transform_domain(d::Domain.Box) = Domain.to_box(transform_domain.(Domain.to_1D_boxes(d)))
