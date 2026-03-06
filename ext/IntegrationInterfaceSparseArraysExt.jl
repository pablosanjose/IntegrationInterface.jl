module IntegrationInterfaceSparseArraysExt

using SparseArrays
using IntegrationInterface

IntegrationInterface.unsafe_deserialize_array(a::SparseMatrixCSC, v::AbstractVector) =
    SparseMatrixCSC(a.m, a.n, a.colptr, a.rowval, convert(Vector, v))

end # module
