module IntegrationInterfaceCubatureExt

using Cubature
using IntegrationInterface
using IntegrationInterface: serialize_array, deserialize_array

(s::IS.Cubature)(f, domain, ::Missing, args; params...) =
	hcubature(point -> f(point, args...; params...), domain...; s.opts...) |> first

# Cubature requires vectors of Float64, so we serialize/deserialize
function (s::IS.Cubature)(f!, domain, result, args; params...)
    v = serialize_array(Float64, result)
    # Could probably use unsafe_deserialize_array here, but we prefer safety.
    fd!(point, out) = f!(deserialize_array(result, out), point, args...; params...)
    v .= first(hcubature(length(v), fd!, domain...; s.opts...))
	return result
end

end # module
