
## Serialize/Deserialize ##
# Needed to convert arrays to vectors of a given type and back, in-place (non-allocating)

serialize_array(::Type{T}, a::AbstractArray) where {T} = reinterpret(T, serialize_array(a))
serialize_array(a::AbstractArray) = vec(a)

deserialize_array(v::AbstractArray, a::AbstractArray{T}) where {T} =
    deserialize_array(reinterpret(T, v), a)
deserialize_array(v::AbstractArray{T}, a::AbstractArray{T}) where {T} =
    (check_deserializer(v, a); unsafe_deserialize_array(v, a))

# assumes equal eltype T and compatible size
unsafe_deserialize_array(v::AbstractArray{T,N}, ::AbstractArray{T,N}) where {T,N} = v
unsafe_deserialize_array(v::AbstractArray, a::AbstractArray) = reshape(v, size(a))

check_deserializer(v, a) =
    size(serialize_array(a)) == size(v) ||
        throw(ArgumentError("Wrong size of serialized array, expected $(size(serialize_array(a))), got $(size(v))"))

check_deserializer(v, a::AbstractVector) =
    length(serialize_array(a)) == length(v) ||
        throw(ArgumentError("Wrong length of serialized array, expected $(length(serialize_array(a))), got $(length(v))"))
