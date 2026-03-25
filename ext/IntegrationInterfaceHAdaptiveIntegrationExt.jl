module IntegrationInterfaceHAdaptiveIntegrationExt

using HAdaptiveIntegration
using HAdaptiveIntegration.StaticArrays
using HAdaptiveIntegration.LinearAlgebra
using IntegrationInterface
import IntegrationInterface as II

function II.convert_domain(d::Domain.Box, ::Backend.HAdaptiveIntegration)
    d´ = II.transform_domain(d)
    x0, x1 = first(d´), last(d´)
    foreach(error_if_inf_HAdaptiveIntegration, (x0..., x1...))
    return Orthotope(SVector(x0), SVector(x1))
end

function II.convert_domain(d::Domain.Simplex, ::Backend.HAdaptiveIntegration)
    vertices = II.Domain.vertices(II.transform_domain(d))
    foreach(error_if_inf_HAdaptiveIntegration, vertices)
    return Simplex(map(SVector, II.Domain.vertices(II.transform_domain(d)))...)
end

II.convert_integrand(i::II.Integral{Nothing}, ::Backend.HAdaptiveIntegration, domain::II.Domain.Box, args; kw...) =
    II.convert_integrand_generic(i, domain, args; post = ensure_real_or_complex, kw...)

II.convert_integrand(i::II.Integral{Nothing}, ::Backend.HAdaptiveIntegration, domain::II.Domain.FiniteRealSimplex, args; kw...) =
    II.convert_integrand_generic(i, domain, args; post = ensure_real_or_complex, kw...)

function II.convert_integrand(i::II.Integral{Nothing}, ::Backend.HAdaptiveIntegration, domain::II.Domain.Simplex{N}, args; kw...) where {N}
    _, basis, _, _ = Domain.basisdata(domain)
    volume = det(hcat(SVector.(basis)...))
    return II.convert_integrand_generic(i, domain, args; post = x -> ensure_real_or_complex(x) * volume, kw...)
end

II.convert_integrand(::II.Integral{<:Any}, ::Backend.HAdaptiveIntegration, domain, args; kw...) =
    throw(ArgumentError("HAdaptiveIntegration does not support in-place integration. Use StaticArrays for array-valued integrands."))

(s::Backend.HAdaptiveIntegration)(f, domain, ::Nothing, witherror) = integrate(f, domain; s.opts...)
(s::Backend.HAdaptiveIntegration)(f, domain, result) = s(f, domain, result, true) |> first

error_if_inf_HAdaptiveIntegration(xs) = foreach(error_if_inf_HAdaptiveIntegration, xs)
error_if_inf_HAdaptiveIntegration(x::Number) = isinf(x) &&
    throw(ArgumentError("The HAdaptiveIntegration backend cannot deal with domains with `Inf`s. Use `Infinity` instead."))

ensure_real_or_complex(x::Union{Real,Complex}) = x
ensure_real_or_complex(x::AbstractArray{<:Union{Real,Complex}}) = x
ensure_real_or_complex(_) = throw(ArgumentError("HAdaptiveIntegration only understand functions with Real or Complex values or eltypes, but not other Numbers like Unitful."))

end # module
