## Fallbacks

transform_domain(d::AbstractDomain) = d

change_of_variables(t, ::AbstractDomain) = t, 1

## 1D Transforms

transform_01(t, x0, x1, t0) = (x0 + (x1 - x0) * t/t0, (x1 - x0)/t0)
transform_0Inf(t, x0, x1, t0) =
    (x0 + (x1 - x0) * (t/t0) / (1 - t/t0),
    (x1 - x0) / (t0 * (1 - t/t0)^2))
transform_InfInf(t::T, x0, x1, t0) where {T} =
    ((x0 + x1)/2 + (x1 - x0) * T(0.75) * (t/t0) / (1 - (t/t0)^2),
    T(0.75) * (x1 - x0)/t0 * (1 + (t/t0)^2) / (1 - (t/t0)^2)^2)

transform1D(t, x0::Number, x1::Infinity, t0) = transform_0Inf(t, x0, point(x1), t0)
transform1D(t, x0::Infinity, x1::Number, t0) = transform_0Inf(t, x1, point(x0), t0) .* (1, -1)
transform1D(t, x0::Infinity, x1::Infinity, t0) = transform_InfInf(t, point(x0), point(x1), t0)
# Anything number that is not real (e.g. complex, unitful) is converted for compatibility
transform1D(t, x0::Number, x1::Number, t0) = transform_01(t, x0, x1, t0)

## Domain.Box{1}

const Box1D{T1,T2} = Domain.Box{1,<:Any,Tuple{T1},Tuple{T2}}

# We don't use floats to preserve Float32 compatibility
transform_domain(::Box1D{<:Number,<:Infinity}) = Domain.Box(0, 1)
transform_domain(::Box1D{<:Infinity,<:Number}) = Domain.Box(0, 1)
transform_domain(::Box1D{<:Infinity,<:Infinity}) = Domain.Box(-1, 1)
transform_domain(::Box1D{<:Number,<:Number}) = Domain.Box(0, 1)
transform_domain(d::Box1D{<:Real,<:Real}) = d                  # no transformation

# Need only in case t is a Tuple{Number} or SVector{1}
change_of_variables(t, d::Box1D) = transform1D(only(t), d, 1)
change_of_variables(t, ::Box1D{<:Real,<:Real}) = (only(t), 1)

transform1D(t, d::Box1D, t0) = transform1D(t, first(d), last(d), t0)

## Domain.Box{N}

transform_domain(d::Domain.Box) = Domain.to_box(transform_domain.(Domain.to_1D_boxes(d)))

function change_of_variables(t, d::Domain.Box)
    xdx = change_of_variables.(t, Domain.to_1D_boxes(d))
    x, dx = first.(xdx), prod(last, xdx)
    return x, dx
end

## Domain.Simplex{N}
# Optimal change of variables requires StaticArrays, a dependency of HAdaptiveIntegration,
# so we define the corresponding dispatches in its extension.

transform_domain(s::Domain.FiniteRealSimplex) = s
transform_domain(s::Domain.InfiniteSimplex) = Domain.orthogonal(s, -1, 1)
transform_domain(s::Domain.Simplex) = Domain.orthogonal(s, 0, 1)

change_of_variables(t, ::Domain.FiniteRealSimplex) = (t, 1)

function change_of_variables(ts, d::Domain.Simplex{N}) where {N}
    origin, basis, boxes, signature = Domain.basisdata(d)
    T = eltype(ts)
    dx = Ref(T(signature))
    t0 = Ref(one(T))
    x = Tuple(simplex_t_to_x.(ts, boxes; dx, t0))
    x´ = origin .+ tuplematvec(basis, x)
    # dx *= det(basis)  # Disabled for performance, see below
    return x´, dx[]
end
# WARNING: there is a missing det(basis) in dx here. We compute and multiply by it only
# at the end for performance.
# This should be done by the backends, preferrably using StaticArrays when available!

function simplex_t_to_x(ti, bi; dx, t0)
    xi, dxi = transform1D(ti, bi, t0[])
    if last(bi) isa Infinity
        t0[] -= ti
    end
    dx[] *= dxi
    return xi
end

# Make do without StaticArrays! mat is tuple of column tuples
tuplematvec(mat::NTuple{N,NTuple{N,Any}}, vec::NTuple{N,Any}) where {N} =
    ntuple(row -> sum(col -> mat[col][row] * vec[col], 1:N), Val(N))
