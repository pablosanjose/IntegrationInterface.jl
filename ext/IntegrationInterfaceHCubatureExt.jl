module IntegrationInterfaceHCubatureExt

using HCubature
using IntegrationInterface
import IntegrationInterface as II

# HCubature assumes a real domain
II.convert_domain(s::Domain.Box, ::Backend.HCubature) = II.convert_domain_generic(error_if_Inf(s))

II.convert_integrand(i::II.Integral{Nothing}, ::Backend.HCubature, domain, args; kw...) =
    II.convert_integrand_generic(i, domain, args; post = ensure_real_or_complex, kw...)

II.convert_integrand(::II.Integral{<:Any}, ::Backend.HCubature, domain, args; kw...) =
    throw(ArgumentError("HCubature does not support in-place integration. Use StaticArrays for array-valued integrands."))

(s::Backend.HCubature)(f, domain, ::Nothing, witherror) = hquadrature_or_hcubature(f, domain...; s.opts...)
(s::Backend.HCubature)(f, domain, result) = s(f, domain, result, true) |> first

hquadrature_or_hcubature(f, min::Number, max::Number; kw...) = hquadrature(f, min, max; kw...)
hquadrature_or_hcubature(args...; kw...) = hcubature(args...; kw...)

error_if_Inf(s::Domain.Box) = any(error_if_Inf, first(s)) || any(error_if_Inf, last(s)) || s
error_if_Inf(x::Number) = isinf(x) &&
    throw(ArgumentError("HCubature doesn't understand domains with `Inf`s. Use `Infinity(point)` instead."))
error_if_Inf(x::Infinity) = false

ensure_real_or_complex(x::Union{Real,Complex}) = x
ensure_real_or_complex(x::AbstractArray{<:Union{Real,Complex}}) = x
ensure_real_or_complex(_) = throw(ArgumentError("HCubature only understand functions with Real or Complex values or eltypes, but not other Numbers like Unitful."))

end # module
