module IntegrationInterfaceCubatureExt

using Cubature
using IntegrationInterface
import IntegrationInterface as II

# supported domains (no need to deal with Domains.Functionals here)
II.check_domain_solver(::Union{Domains.Box}, ::Backends.Cubature) = nothing

II.convert_domain(s::Domains.Box, ::Backends.Cubature) = (s.mins, s.maxs)

(s::Backends.Cubature)(f, domain, ::Missing, args; params...) =
	hcubature(point -> f(point, args...; params...), II.convert_domain(domain, s)...; s.opts...) |> first

# Cubature requires vectors of Float64, so we serialize/deserialize
function (s::Backends.Cubature)(f!, domain, result, args; params...)
    v = II.serialize_array(Float64, result)
    # Could probably use unsafe_deserialize_array here, but we prefer safety.
    fd!(point, out) = f!(II.deserialize_array(result, out), point, args...; params...)
    v .= first(hcubature(length(v), fd!, convert_domain(domain, s)...; s.opts...))
	return result
end

end # module
