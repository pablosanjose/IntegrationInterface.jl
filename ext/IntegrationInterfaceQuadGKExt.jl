module IntegrationInterfaceQuadGKExt

using QuadGK
using IntegrationInterface
import IntegrationInterface as II

II.convert_domain(s::Domain.Box{1}, ::Backend.QuadGK) =
    II.convert_domain_generic(s)

II.convert_integrand(i::II.Integral{<:Any}, ::Backend.QuadGK, domain, args; kw...) =
    II.convert_integrand_generic(i, domain, args; kw...)

## Call ##

(s::Backend.QuadGK)(f, domain, ::Nothing, witherror) = quadgk(f, domain...; s.opts...)
(s::Backend.QuadGK)(f!, domain, result, witherror) = quadgk!(f!, result, domain...; s.opts...)
(s::Backend.QuadGK)(f, domain, result) = s(f, domain, result, true) |> first

end # module
