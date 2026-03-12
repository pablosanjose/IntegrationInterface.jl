module IntegrationInterfaceCubatureExt

using Cubature
using IntegrationInterface
import IntegrationInterface as II

II.convert_domain(s::Domain.Box, ::Backend.Cubature) = II.convert_domain_generic(s)

# Scalar case is easy
II.convert_integrand(i::II.Integral{Missing,<:Backend.Cubature}, domain, args; params...) =
    II.convert_integrand_generic(i, domain, args; params...)

# Cubature only understands Scalar or Vector{Float32} arguments. We must serialize/deserialize
function II.convert_integrand(i::II.Integral{<:Any,<:Backend.Cubature}, d, args; params...)
    f! = II.integrand(i)
    fd!(t, out) = jacobian(t, d) * f!(II.deserialize_array(t, out), change_of_variables(t, d)..., args...; params...)
    return fd!
end

(s::Backend.Cubature)(f, domain, ::Missing) = hcubature(f, domain...; s.opts...) |> first

# Cubature requires vectors of Float64, so we serialize/deserialize
function (s::Backend.Cubature)(f!, domain, result)
    v = II.serialize_array(Float64, result)
    v .= first(hcubature(length(v), f!, domain...; s.opts...))
	return result
end

end # module
