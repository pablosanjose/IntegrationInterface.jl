module IntegrationInterfaceHCubatureExt

using HCubature
using IntegrationInterface

(s::IS.HCubature)(f, domain, ::Missing, args; params...) =
	hcubature(point -> f(point, args...; params...), domain...; s.opts...) |> first

(::IS.HCubature)(f!, domain, result, args; params...) =
    throw(ArgumentError("HCubature does not support in-place integration. Use StaticArrays for array-valued integrands."))

end # module
