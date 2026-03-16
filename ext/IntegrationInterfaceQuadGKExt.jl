module IntegrationInterfaceQuadGKExt

using QuadGK
using IntegrationInterface
import IntegrationInterface as II

II.convert_domain(s::Domain.Line, ::Backend.QuadGK) =
    II.convert_domain_generic(s)

II.convert_integrand(i::II.Integral{<:Any,<:Backend.QuadGK}, domain, args; kw...) =
    II.convert_integrand_generic(i, domain, args; kw...)

## Call ##

(s::Backend.QuadGK)(f, domain, ::Nothing) = quadgk(f, domain...; s.opts...) |> first
(s::Backend.QuadGK)(f!, domain, result) = quadgk!(f!, result, domain...; s.opts...) |> first

end # module
