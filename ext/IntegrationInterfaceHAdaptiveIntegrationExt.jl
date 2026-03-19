module IntegrationInterfaceHAdaptiveIntegrationExt

using HAdaptiveIntegration
using HAdaptiveIntegration.StaticArrays
using IntegrationInterface
import IntegrationInterface as II

II.convert_domain(s::Domain.Simplex, ::Backend.HAdaptiveIntegration) =
    Simplex(map(SVector, II.vertices(II.transform_domain(s))))

II.transform_domain(s::II.Domain.RealFiniteSimplex{N}) where {N} = s
# normalize to unit simplex
II.transform_domain(s::II.Domain.Simplex) = II.Domain.oneunit(s)

# simplex basis
function basis(s::II.Domain.Simplex) where {N}
    z, us... = II.point.(II.Domain.vertices(s))
    return ntuple(i-> SVector(us[i] .- z), Val(N))
end

volume(basis) = det(hcat(basis...))

function II.convert_integrand(i::II.Integral{<:Any,<:Backend.HAdaptiveIntegration}, domain::II.Domain.Simplex, args; kw...)
    b = basis(domain)
    ortho = Box(b...)
    v = volume(b)
    return II.convert_integrand_generic(i, ortho, args; post = x->x*v, kw...)
end

(s::Backend.HAdaptiveIntegration)(f, domain, ::Nothing) = integrate(f, domain...; s.opts...) |> first

end # module
