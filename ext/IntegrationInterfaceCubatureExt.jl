module IntegrationInterfaceCubatureExt

using Cubature
using IntegrationInterface
import IntegrationInterface as II

II.convert_domain(s::Domain.Box, ::Backend.Cubature) = II.convert_domain_generic(error_if_Inf(s))

# Scalar case is easy
II.convert_integrand(i::II.Integral{Nothing,<:Backend.Cubature}, domain, args; params...) =
    II.convert_integrand_generic(i, domain, args; post = ensure_real, params...)

# Cubature only understands Scalar or Vector{Float32} arguments. We must serialize/deserialize
function II.convert_integrand(i::II.Integral{<:Any,<:Backend.Cubature}, d, args; params...)
    f! = II.integrand(i)
    result = II.result(i)
    function fd!(t, out)
        f!(II.deserialize_array(out, result), II.change_of_variables(t, d)..., args...; params...)
        II.deserialize_array(out, result) .*= II.jacobian(t, d)
        return out
    end
    return fd!
end

(s::Backend.Cubature)(f, domain, ::Nothing) = hcubature(f, domain...; s.opts...) |> first

# Cubature requires vectors of Float64, so we serialize/deserialize
function (s::Backend.Cubature)(f!, domain, result)
    v = II.serialize_array(Float64, result)
    v .= first(hcubature(length(v), f!, domain...; s.opts...))
	return result
end

ensure_real(x::Real) = x
ensure_real(_) = throw(ArgumentError("Cubature doesn't understand complex-valued functions. You can try with a mutating complex-vector-valued function."))

error_if_Inf(s::Domain.Box) = any(error_if_Inf, first(s)) || any(error_if_Inf, last(s)) || s
error_if_Inf(x::Number) = isinf(x) &&
    throw(ArgumentError("Cubature doesn't understand domains with `Inf`s. Use `Infinity(point)` instead."))
error_if_Inf(x::Infinity) = false

end # module
