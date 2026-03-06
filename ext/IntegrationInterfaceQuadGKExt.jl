module IntegrationInterfaceQuadGKExt

using QuadGK
using IntegrationInterface

(s::IS.QuadGK)(f, domain, ::Missing, args; params...) =
	quadgk(x -> f(x, args...; params...), domain; s.opts...) |> first
(s::IS.QuadGK)(f!, domain, result, args; params...) =
	quadgk!((out, x) -> f!(out, x, args...; params...), result, domain; s.opts...) |> first

end # module
