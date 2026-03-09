module IntegrationInterfaceQuadGKExt

using QuadGK
using IntegrationInterface
import IntegrationInterface as II

# supported domains (no need to deal with Domains.Functionals here)
II.check_domain_solver(::Union{Domains.Segment, Domains.SegmentGroup}, ::Backends.QuadGK) = nothing

II.convert_domain(s::Domains.Segment, ::Backends.QuadGK) = (s.x1, s.x2)
II.convert_domain(s::Domains.SegmentGroup, ::Backends.QuadGK) = (s.segments,)

## Call ##

(s::Backends.QuadGK)(f, domain, ::Missing, args; params...) =
	quadgk(x -> f(x, args...; params...), II.convert_domain(domain, s)...; s.opts...) |> first
(s::Backends.QuadGK)(f!, domain, result, args; params...) =
	quadgk!((out, x) -> f!(out, x, args...; params...), result, II.convert_domain(domain, s)...; s.opts...) |> first

end # module
