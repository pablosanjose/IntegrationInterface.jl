module IntegrationInterfaceHCubatureExt

using HCubature
using IntegrationInterface
import IntegrationInterface as II

II.convert_domain(s::Domain.Box, ::Backend.HCubature) = II.convert_domain_generic(s)

II.convert_integrand(i::II.Integral{Missing,<:Backend.HCubature}, domain, args; params...) =
    II.convert_integrand_generic(i, domain, args; params...)

II.convert_integrand(::II.Integral{<:Any,<:Backend.Cubature}, domain, args; params...) =
    throw(ArgumentError("HCubature does not support in-place integration. Use StaticArrays for array-valued integrands."))

(s::Backend.HCubature)(f, domain, ::Missing) = hcubature(f, domain...; s.opts...) |> first

end # module
