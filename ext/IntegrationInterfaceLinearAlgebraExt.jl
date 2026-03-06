module IntegrationInterfaceLinearAlgebraExt

using LinearAlgebra
using IntegrationInterface

IntegrationInterface.unsafe_deserialize_array(::Diagonal, v::AbstractVector) =
    Diagonal(convert(Vector, v))

end # module
