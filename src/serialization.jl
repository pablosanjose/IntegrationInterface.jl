
## Serialize/Deserialize ##
# Needed to convert arrays to vectors of a given type and back, in-place (non-allocating)

serialize_array(::Type{T}, a::AbstractArray) where {T} = reinterpret(T, serialize_array(a))
serialize_array(a::AbstractArray) = vec(a)

deserialize_array(a::AbstractArray{T}, v::AbstractArray) where {T} =
    deserialize_array(a, reinterpret(T, v))
deserialize_array(a::AbstractArray{T}, v::AbstractArray{T}) where {T} =
    (check_deserializer(a, v); unsafe_deserialize_array(a, v))

# assumes equal eltype T and compatible size
unsafe_deserialize_array(::AbstractArray{T,N}, v::AbstractArray{T,N}) where {T,N} = v
unsafe_deserialize_array(a::AbstractArray, v::AbstractArray) = reshape(v, size(a))

check_deserializer(a, v) =
    size(serialize(a)) == size(v) ||
        argerror("Wrong size of serialized array, expected $(size(serialize(a))), got $(size(v))")

check_deserializer(a, v::AbstractVector) =
    length(serialize(a)) == length(v) ||
        argerror("Wrong length of serialized array, expected $(length(serialize(a))), got $(length(v))")
