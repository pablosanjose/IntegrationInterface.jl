module IntegrationInterfaceQuadGKExt

using QuadGK
using IntegrationInterface
import IntegrationInterface as II

# supported domains (no need to deal with Domain.Functionals here)
II.check_domain_solver(::Type{<:Union{Domain.Segment, Domain.SegmentGroup}}, ::Backend.QuadGK) = nothing

II.convert_domain(s::Domain.Segment, ::Backend.QuadGK) = (s.x1, s.x2)
II.convert_domain(s::Domain.SegmentGroup, ::Backend.QuadGK) = (s.segments,)

## Call ##

(s::Backend.QuadGK)(f, domain, ::Missing, args; params...) =
	quadgk(x -> f(x..., args...; params...), II.convert_domain(domain, s, args)...; s.opts...) |> first
(s::Backend.QuadGK)(f!, domain, result, args; params...) =
	quadgk!((out, x) -> f!(out, x..., args...; params...), result, II.convert_domain(domain, s, args)...; s.opts...) |> first

end # module
