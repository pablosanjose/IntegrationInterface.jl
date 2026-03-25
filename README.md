# IntegrationInterface

[![Build Status](https://github.com/pablosanjose/IntegrationInterface.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/pablosanjose/IntegrationInterface.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package aims to be a lightweight, faster-loading alternative to the excellent [Integrals.jl](https://github.com/SciML/Integrals.jl) from the SciML ecosystem. IntegrationInterface.jl offers an interface to perform n-dimensional numerical integrals of scalars, arrays or more general objects, over a domain $D$.

$$J(\text{args}...; \text{kw}...) = \int_{D(\text{args}...; \text{kw}...)} d^n x f(\boldsymbol x..., \text{args}...; \text{kw}...)$$

The general interface reads
```julia
julia> J = Integral(f, domain; backend::AbstractBackend = Backend.default(domain), result = nothing)
```
This produces an `J::Integral` object representing the integral of `f` over a given `domain`.

Mutating functions `f!(out, x..., args...; kw...)` that modify an `out::A` in-place can also be used. This is useful for heap-allocated integrands of type e.g. `A::AbstractArray`. In this pass an array of type `A` with the `result` keyword.

To compute the integral for a set of `args` and `kw`, use the syntax
```julia
julia> J(args...; kw...)
```
In the mutating case, this will also write the value of the integral into `result`.

If one doesn't need to evaluate `J` repeatedly for different `args` and `kw`, or use `J` inside an integral (see Composition below), one may use the following lowercase `integral` form to obtain the value of the integral directly (equivalent to building *and* calling `J(args...; kw...)`),
```julia
julia> integral(f, domain, args...; backend::AbstractBackend = Backend.default(domain), result = nothing, kw...)
```

## Domains

We currently support bounded and unbounded hypercube and simplex domains.

- Hypercubes are created with `Domain.Box((a₁, a₂, ...), (b₁, b₂, ...))` or `Domain.interval((a₁, b₁), (a₂, b₂)...)` in terms of intervals `(aᵢ, bᵢ)` along dimension `i`.

- Simplices are defined in terms of n+1 vertices in n-dimensional space, with `Domain.Simplex(v₁, v₂,...,vₙ₊₁)`.

Functional domains that depend on `args` and `kw` can be constructed with e.g. `D(args...; kw...) = Domain.Box(...)`, etc, see below for an example. Functional domains are evaluated when calling `J(args...; kw...)`.

## Simple examples
### 1D integral
```julia
julia> using QuadGK
[ Info: Precompiling IntegrationInterfaceQuadGKExt [6a486dfe-6a5b-5d49-a9f9-02f4245ab8d6]

julia> integral(cos, Domain.interval(0, π/2)) # Defaults to Backend.QuadGK()
1.0
```

### 2D integral
```julia
julia> using HCubature, HAdaptiveIntegration
[ Info: Precompiling IntegrationInterfaceQuadGKExt [6a486dfe-6a5b-5d49-a9f9-02f4245ab8d6]
[ Info: Precompiling IntegrationInterfaceHCubatureExt [b56c1907-5b79-5f70-9c8b-cc15ee23dcd1]
[ Info: Precompiling IntegrationInterfaceHAdaptiveIntegrationExt [6c915174-9fe5-54f3-904c-38242f473220]

julia> integral((x,y) -> cos(x-y), Domain.Box((-1,-1), (1,1))) # Defaults to Backend.HCubature()
2.8322936730937722

julia> integral((x,y) -> cos(x-y), Domain.Simplex((0,0), (1, 1/2), (1/2, 1))) # Defaults to Backend.HAdaptiveIntegration()
0.36725231432888195
```

## Backends

As shown above, the integration is actually performed by backend packages that may be loaded as needed. Currently supported packages (weak dependencies) and corresponding backends in the `Backends` submodule are

- QuadGK.jl: `Backend.QuadGK(; opts...)` (calls `quadgk` and `quadgk!`, default for `Domain.Box{1}` domains)
- HCubature.jl:  `Backend.HCubature(; opts...)` (calls `hcubature`, default for `Domain.Box{N}` domains with `N ≠ 1`)
- Cubature.jl: `Backend.Cubature(; opts...)` (calls `hcubature`)
- HAdaptiveIntegration.jl: `Backend.HAdaptiveIntegration(; opts...)` (calls `integrate`, default for `Domain.Simplex`)

We also provide a `Backend.Quadrature((nodes, weights))` backend that can be used e.g. with the FastGaussQuadrature.jl package that computes nodes and weights for a 1D integral in the [-1, 1] integration domain. The `Backend.Quadrature` solver then uses these values for integrals over any `Domain.Box{N}` for any `N`.
```julia
julia> using FastGaussQuadrature, QuadGK

julia> f(x) = cos(x);

julia> J1 = Integral(f, Domain.interval(3,5); backend = Backend.Quadrature(gausslegendre(10)))
Integral
  Mutating   : false
  Domain     : Box{1,Float64}(3.0, 5.0)
  Backend    : Quadrature
  Integrand  : f

julia> J2 = Integral(f, Domain.interval(3,5); backend = Backend.QuadGK())
Integral
  Mutating   : false
  Domain     : Box{1,Float64}(3.0, 5.0)
  Backend    : QuadGK
  Integrand  : f

julia> J1(), J2()
(-1.1000442827230061, -1.1000442827230057)
```

## Composition of `Integral`s
We can also express nested integrals such as

$$J(\text{args}...; \text{kw}...) = \int_{D_n} dx_n\dots\int_{D_1} dx_1 f(\boldsymbol{x}, \text{args}...; \text{kw}...)$$

As a concrete example, consider

$$J = \int_2^3 dy\int_0^1 dx (x-y)^2\cos(x+y) $$

This integral can be evaluated either as one adaptive `HCubature` or two nested `QuadGK` (the default backend):
```julia
julia> f(x,y) = (x-y)^2 * cos(x+y)
f (generic function with 2 methods)

julia> J1 = f |> Integral(Domain.interval(0,1)) |> Integral(Domain.interval(2,3))
Integral
  Mutating   : false
  Domain     : Box{1,Float64}(2.0, 3.0)
  Backend    : QuadGK
  Integrand  : Integral
    Mutating   : false
    Domain     : Box{1,Float64}(0.0, 1.0)
    Backend    : QuadGK
    Integrand  : f

julia> J2 = Integral(f, Domain.interval((0,1), (2,3)); backend = Backend.HCubature())
Integral
  Mutating   : false
  Domain     : Box{2, Float64}((0.0, 2.0), (1.0, 3.0)))
  Backend    : HCubature
  Integrand  : f

julia> (J1(), J2())
(-3.800374064781164, -3.800374064812097)
```
Note the currying syntax used above for `J1`. It is equivalent to `J1 = Integral(Integral(f, Domain.Box(0,1)), Domain.Box(2,3))`.

We can combine functional domains and integral composition to integrate over non-Box domains. To do so, we make the inner domain depend on outer integration variables. For example,

$$J = \int_{-1}^1 dy\int_{-\sqrt{1-y^2}}^{\sqrt{1-y^2}} dx (x-y)^2\cos(x+y) $$

can be expressed as
```julia
julia> f(x,y) = (x-y)^2 * cos(x+y);

julia> J = f |> Integral(y -> Domain.Box(-sqrt(1-y^2), sqrt(1-y^2))) |> Integral(Domain.Box(-1,1))
Integral
  Mutating   : false
  Domain     : Box{1,Float64}(-1.0, 1.0)
  Backend    : Default
  Integrand  : Integral
    Mutating   : false
    Domain     : Functional
    Backend    : Default
    Integrand  : f

julia> J()
1.324825188363749

```

## Infinity

Only a few backends support using `Inf` to express unbounded domains. For those that don't support it, we provide `Infinity(point::Number)` that can be used in domain bounds instead. `Infinity(point)` represents an unbounded ray passing though `point` (which may be Real or not). These `Infinite` bounds are dealt with using an appropriate change of variables that takes `point` into account. As an example, consider a 2D half-plane `(; σ = 1) -> Domain.Box((0, -Infinity(σ)), (Infinity(σ), Infinity(σ)))`. Note that it is a functional domain that depends on a keyword argument `σ`. We can integrate a Gaussian of width `σ` over `D(σ)`, which gives `π` for `σ = 1`

```julia
julia> f(x, y; σ = 1) = exp(-0.5*(x^2+y^2)/σ^2);

julia> D(; σ = 1) = Domain.Box((0, -Infinity(σ)), (Infinity(σ), Infinity(σ)));

julia> J = Integral(f, D)
Integral
  Mutating   : false
  Domain     : Functional
  Backend    : Default
  Integrand  : f

julia> J(; σ = 1)
3.141592652846476
```

Simplices can be made unbounded by wrapping one or more (but not all) of its vertices `vᵢ` in `Infinity(vᵢ...)`. This corresponds to moving the vertex to infinity along a ray passing through it in the direction perpedicular to the opposite facet. For example, this integrates `exp(-z)` in a semi-infinite vertical prism with a triangular base
```julia
julia> integral((x,y,z) -> exp(-z), Domain.Simplex((1,0,0), (-1,-1,0), (-1,1,0), Infinity(0,0,1)))
-2.000000000125864
```
Note that the integral preserves the sign of the simplex volume, which is negative in this case due to the order chosen for the vertices.

## Error estimation

Many backends support computing both the value of the integral and its estimated error. To access the latter we may use the `witherror` command. Instead of evaluating a `J::Integral` object directly with `J(args...; kw...)`, we do `witherror(J, args...; kw...)` or `J |> witherror(args...; kw...)`, which returns `(value, error)`

Example:
```julia
julia> using QuadGK

julia> J = Integral((x; λ = 1) -> exp(-x/λ), (; λ = 1) -> Domain.Box(0, Infinity(λ)))
Integral
  Mutating   : false
  Domain     : Functional
  Backend    : Default
  Integrand  : #26

julia> value, error = witherror(J; λ = 2)
(2.0, 9.014765349073633e-11)

julia> value, error = J |> witherror(; λ = 2)   # alternative `currying` form
(2.0, 9.014765349073633e-11)
```
