module IntegrationInterfaceHCubatureExt

using HCubature
using IntegrationInterface
import IntegrationInterface as II

# supported domains (no need to deal with Domains.Functionals here)
II.check_domain_solver(::Domains.Box, ::Backends.HCubature) = nothing

II.convert_domain(s::Domains.Box, ::Backends.HCubature) = (s.mins, s.maxs)

(s::Backends.HCubature)(f, domain, ::Missing, args; params...) =
	hcubature(point -> f(point, args...; params...), II.convert_domain(domain, s)...; s.opts...) |> first

(::Backends.HCubature)(f!, domain, result, args; params...) =
    throw(ArgumentError("HCubature does not support in-place integration. Use StaticArrays for array-valued integrands."))

end # module
