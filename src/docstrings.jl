"""
`Domain` is a submodule of `IntegrationInterface` that bundles possible integration domains.
Currently it provides
    - `Domain.Box{N}`: N-dimensional hypercube domain that may be bounded or unbounded.
    - `Domain.Simplex{N}`: N-dimensional simplex domain that may be bounded or unbounded.

It also provides the `Domain.interval` function, an alternative way to build `Domain.Box`
domains
"""
Domain

"""
`Backend` is a submodule of `IntegrationInterface` that bundles available integration
backends, each relying on an external package that must be loaded for use.
Currently it provides
    - `Backend.QuadGK(; opts...)`: 1-dimensional integrals `using QuadGK` (calls `quadgk`)
    - `Backend.HCubature(; opts...)`: N-dimensional integrals `using HCubature` (calls `hcubature`)
    - `Backend.Cubature(; opts...)`: N-dimensional integrals `using Cubature` (calls `hcubature`)
    - `Backend.HAdaptiveIntegration(; opts...)`: N-dimensional integrals `using HAdaptiveIntegration` (calls `integrate`)
    - `Backend.Quadrature((nodes, weights))`: N-dimensional integrals using user-provided nodes and weights.
        - Nodes and weights are assumed to correspond to interval `[-1,1]`.
        - These can be computed with external packages, such as FastGaussQuadrature.jl
        - The form `Backend.Quadrature(nodes, weights)` is also provided.

Options `opts` are passed to the corresponding integration routine. Check the package
documentation for details.

"""
Backend

"""
    Integral(f::Function, domain; backend = default_backend(domain), result = nothing)

Create a `J::Integral` representing the integral of `f(x₁, ..., xₙ, args...; kw...)` over
`domain` in `n`-dimensional space. The backend is an object from the `Backend` submodule,
which may also depend on `args` and `kw`. The default backend is `Domain.QuadGK` for 1D domains, and
`Domain.HCubature` for higher dimensions. To evaluate the integral for a given set of `args`
and `kw`, do `J(args...; kw...)`. See also `integral` for direct evaluation.

An in-place form for `f`, namely `f!(out, x₁, ..., xₙ, args...; kw...)` that overwrites
`out` in-place with the function value of type `A`, can also be used. This mode requires
that a `result` of type `A` is provided, which will be overwritten with the result of the
integral upon evaluation.

    f |> Integral(domain; kw...)

Currying form of `integral`, equivalent to `Integral(f, domain; kw...)`

    Integral(J::Integral, domain; kw...)

Create a nested integral. Each level can have different `kw` and `domain` (see example).

# Example
```julia
julia> using QuadGK, HCubature

julia> J = Integral(x -> exp(-x), Domain.Box{1}(0, Infinity(1)))  # uses QuadGK by default
Integral
  Mutating   : false
  Domain     : Box{1,Float64}(0.0, Infinity(1.0))
  Backend    : QuadGK
  Integrand  : #61

julia> J()
1.0

julia> f(x,y) = exp(-x^2-y^2);

julia> J = f |> Integral(Domain.Box{1}(0,1)) |> Integral(Domain.Box{1}(0, Infinity(1))) # Nested!
Integral
  Mutating   : false
  Domain     : Box{1,Float64}(0.0, Infinity(1.0))
  Backend    : QuadGK
  Integrand  : Integral
    Mutating   : false
    Domain     : Box{1,Float64}(0.0, 1.0)
    Backend    : QuadGK
    Integrand  : f

julia> J´ = Integral(f, Domain.Box((0,0), (1, Infinity(1))))
Integral
  Mutating   : false
  Domain     : Box{2, Float64}((0.0, 0.0), (1.0, Infinity(1.0))))
  Backend    : HCubature
  Integrand  : f

julia> J(), J´()
(0.6618556550762793, 0.6618556550778879)
```
"""
Integral

"""
    integral(f, domain, args...; backend = default_backend(domain), result = nothing, kw...)

Create and immediately evaluate a `J::Integral` object using `args` and `kw`. This is
formally equivalent to `Integral(f, domain; backend, result)(args...; kw...)`.

# Example
```julia
julia> using QuadGK

julia> integral(x -> exp(-x), Domain.Box{1}(0, Infinity(1)))
1.0
```

# See also:
    `Integral`
"""
integral

"""
    Domain.Box((x₁ᵐⁱⁿ, x₂ᵐⁱⁿ, ..., xₙᵐⁱⁿ), (x₁ᵐᵃˣ, x₂ᵐᵃˣ, ..., xₙᵐᵃˣ))

Create an integration domain `Domain.Box{N}` for a function `f(x₁, x₂, ..., xₙ)` over an `N`
dimensional hypercube defined by the intervals `(xᵢᵐⁱⁿ, xᵢᵐᵃˣ)`. The default backend for
`Domain.Box{1}` is QuadGK, and for higher `N` it is HCubature.

    Domain.Box{N}(f::Function)

Create a `D::Domain.Functional{Box{N}}` domain that depends on external parameters.
Evaluating it as `D(args...; kw...)` produces `Box(f(args...; kw...)...)`, which should be
a `Box{N}`. Note that `Domain.Box(::Function)` is not supported, `N` must be specified.

# Examples
```julia
julia> using QuadGK, HCubature

julia> J = Integral((x; kw...) -> exp(-x), Domain.Box{1}((; x0 = 1) -> (x0, Infinity(2x0))))
Integral
  Mutating   : false
  Domain     : Functional{Box{1}}
  Backend    : QuadGK
  Integrand  : #53

julia> J(x0 = 3)
0.04978706836786264

julia> integral((x,y,z) -> cos(x+y+z), Domain.Box((0, 0, 0), (π/2, π/2, π/2)))
-1.9999999998615692

```

# See also:
    `Domain.interval`, `Domain.Simplex`, `Infinity`, `integral`

"""
Domain.Box

"""
    Domain.interval(a::Union{Number,Infinity}, b::Union{Number,Infinity})

Equivalent to `Domain.Box{1}(a, b)`. Construct a 1D domain from `a` to `b`

    Domain.interval((a₁, b₁), (a₂, b₂), ...)

Equivalent to `Domain.Box((a₁, a₂, ...), (b₁, b₂, ...))`. Construct an N-D hypercube domain
spanning interval `(aᵢ, bᵢ)` along the `i` dimension.

# Example
```julia
julia> using QuadGK

julia> integral(cos, Domain.interval(0, π/2))
1.0
```
"""
Domain.interval

"""
    Domain.Simplex((x₁, y₁, ...), (x₂, y₂, ...), ...)

Create an integration domain `Domain.Simplex{N}` of a Simplex (triangle (N=2), tetrahedron
(N=3) and higher-dimensional generalizations), defined in terms of `N+1` vertices of `N`
coordinates `vᵢ = (xᵢ, yᵢ, ...)`. One or more (but not all) of the `vᵢ` can be wrapped in
`Infinity(vᵢ...)` to construct an unbounded simplex, whose vertex `vᵢ` is shifted to
infinity along a ray passing through `vᵢ` in a direction perpendicular to its opposite
facet. The default backend for `Domain.Simplex` is `HAdaptiveIntegration`.

    Domain.Simplex{N}(f::Function)

Create a `D::Domain.Functional{Box{N}}` domain that depends on external parameters.
Evaluating it as `D(args...; kw...)` produces `Simplex(f(args...; kw...)...)`, which should
be a `Simplex{N}`. Note that `Domain.Simplex(::Function)` is not supported, `N` must be
specified.

# Examples
```julia
julia> using HAdaptiveIntegration

julia> J = Integral((x,y,z) -> exp(-z), Domain.Simplex((0,0,0), (1,0,0), (0,1,0), Infinity(0,0,1)))
Integral
  Mutating   : false
  Domain     : Simplex{3, Float64}((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), Infinity(0.0, 0.0, 1.0)))
  Backend    : HAdaptiveIntegration
  Integrand  : #17

julia> J()
0.500000000031466
```
# See also:
    `Domain.Box`, `Infinity`, `integral`

"""
Domain.Simplex

"""
    Infinity(x::Number)

Creates an `Infinity` object representing an unbounded 1D ray passing through `x`. It can be
used in place of `Inf` to create unbounded domains, e.g. `Domain.Box{1}(0, Infinity(1+im))`,
which represents a straight semi-infinite line in the complex plane starting at 0 and
extending to infinity through point `1+im`.

    Infinity(x₁, x₂, ..., xₙ)

Creates an `Infinity` object passing through a vertex in `n`-dimensional space. Useful to
build unbounded Simplex domains.

# Example
```julia
julia> using HCubature

julia> J = Integral((x, y; σ) -> 1/(2π*σ^2) * exp(-(x^2+y^2)/(2*σ^2)), Domain.Box{2}((; σ) -> ((0,0), (Infinity(σ), Infinity(σ)))))
Integral
  Mutating   : false
  Domain     : Functional{Box{2}}
  Backend    : HCubature
  Integrand  : #47

julia> J(σ = 1e5)
0.2499999998823345
```
"""
Infinity
