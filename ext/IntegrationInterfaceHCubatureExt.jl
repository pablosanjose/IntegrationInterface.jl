module IntegrationInterfaceHCubatureExt

using HCubature
using IntegrationInterface
import IntegrationInterface as II

# supported domains (no need to deal with Domain.Functionals here)
II.check_domain_solver(::Domain.Box, ::Backend.HCubature) = nothing

II.convert_domain(s::Domain.Box, ::Backend.HCubature) = (s.mins, s.maxs)

(s::Backend.HCubature)(f, domain, ::Missing, args; params...) =
	hcubature(point -> f(point, args...; params...), II.convert_domain(domain, s)...; s.opts...) |> first

(::Backend.HCubature)(f!, domain, result, args; params...) =
    throw(ArgumentError("HCubature does not support in-place integration. Use StaticArrays for array-valued integrands."))

end # module
