module IntegrationInterfaceHCubatureExt

using HCubature
using IntegrationInterface
import IntegrationInterface as II

# HCubature assumes a real domain
II.convert_domain(s::Domain.Box, ::Backend.HCubature) = II.convert_domain_generic(error_if_Inf(s))

II.convert_integrand(i::II.Integral{Missing,<:Backend.HCubature}, domain, args; params...) =
    II.convert_integrand_generic(i, domain, args; params...)

II.convert_integrand(::II.Integral{<:Any,<:Backend.HCubature}, domain, args; params...) =
    throw(ArgumentError("HCubature does not support in-place integration. Use StaticArrays for array-valued integrands."))

(s::Backend.HCubature)(f, domain, ::Missing) = hcubature(f, domain...; s.opts...) |> first

error_if_Inf(s::Domain.Box) = any(error_if_Inf, first(s)) || any(error_if_Inf, last(s)) || s
error_if_Inf(x::Number) = isinf(x) &&
    throw(ArgumentError("HCubature doesn't understand domains with `Inf`s. Use `Infinity(point)` instead."))
error_if_Inf(x::Infinity) = false

end # module
