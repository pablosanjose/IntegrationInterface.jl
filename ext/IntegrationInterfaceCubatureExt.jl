module IntegrationInterfaceCubatureExt

using Cubature
using IntegrationInterface
import IntegrationInterface as II

II.convert_domain(s::Domain.Box, ::Backend.Cubature) = II.convert_domain_generic(error_if_Inf(s))

# Scalar case is easy
II.convert_integrand(i::II.Integral{Nothing}, ::Backend.Cubature, domain, args; kw...) =
    II.convert_integrand_generic(i, domain, args; post = ensure_real, kw...)

# Cubature only understands Scalar or Vector{Float32} arguments. We must serialize/deserialize
function II.convert_integrand(i::II.Integral{<:Any}, ::Backend.Cubature, d, args; kw...)
    f! = II.integrand(i)
    result = II.result(i)
    function fd!(t, out)
        x, dx = II.change_of_variables(t, d)
        f!(II.deserialize_array(out, result), x..., args...; kw...)
        II.deserialize_array(out, result) .*= dx
        return out
    end
    return fd!
end

(s::Backend.Cubature)(f, domain, ::Nothing, witherror) = hquadrature_or_hcubature(f, domain...; s.opts...)

# Cubature requires vectors of Float64, so we serialize/deserialize
function (s::Backend.Cubature)(f!, domain, result, witherror)
    v = II.serialize_array(Float64, result)
    value, error = hquadrature_or_hcubature(length(v), f!, domain...; s.opts...)
    v .= value
	return result, error
end

hquadrature_or_hcubature(f, min::Number, max::Number; kw...) = hquadrature(f, min, max; kw...)
hquadrature_or_hcubature(n, f!, min::Number, max::Number; kw...) = hquadrature(n, f!, min, max; kw...)
hquadrature_or_hcubature(args...; kw...) = hcubature(args...; kw...)

(s::Backend.Cubature)(f, domain, result) = first(s(f, domain, result, true))

ensure_real(x::Real) = x
ensure_real(_) = throw(ArgumentError("Cubature doesn't understand non-real functions. You can try with a mutating complex-vector-valued function."))

error_if_Inf(s::Domain.Box{1}) = error_if_Inf(first(s)) || error_if_Inf(last(s)) || s
error_if_Inf(s::Domain.Box) = any(error_if_Inf, first(s)) || any(error_if_Inf, last(s)) || s
error_if_Inf(x::Number) = isinf(x) &&
    throw(ArgumentError("Cubature doesn't understand domains with `Inf`s. Use `Infinity(point)` instead."))
error_if_Inf(x::Infinity) = false

end # module
