# IntegrationInterface

[![Build Status](https://github.com/pablosanjose/IntegrationInterface.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/pablosanjose/IntegrationInterface.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package aims to be a lightweight, faster-loading alternative to the excellent [Integrals.jl](https://github.com/SciML/Integrals.jl) from the SciML ecosystem. IntegrationInterface.jl offers an interface to perform n-dimensional numerical integrals of scalars, arrays or more general objects, over a domain $D$.

$$J(\text{args}...; \text{kw}...) = \int_{D(\text{args}...; \text{kw}...)} d^n x f(\boldsymbol x..., \text{args}...; \text{kw}...)$$

The general interface reads
```julia
julia> J = Integral(f, domain; backend::AbstractBackend = default_backend(domain), result = nothing)
```
This produces an `J::Integral` object. Currently, only bounded and unbounded hypercube domains are possible. They are created with `Domain.Box((a₁, a₂, ...), (b₁, b₂, ...))` or `Domain.interval((a₁, b₁), (a₂, b₂)...)` in terms of intervals `(aᵢ, bᵢ)` along dimension `i`. Functions of `args` can be passed to the constructor of a domain to make it depend on arguments and keywords passed to `J`, see `Domain.Box`.

Mutating functions `f!(out, x..., args...; kw...)` that modify an `out::A` in-place can also be used. This is useful for heap-allocated integrands of type e.g. `A::AbstractArray`. In this pass an array of type `A` with the `result` keyword.

To compute the integral for a set of `args` and `kw`, use the syntax
```julia
julia> J(args...; kw...)
```
In the mutating case, this will also write the value of the integral into `result`.

If one doesn't need to evaluate `J` repeatedly for different `args` and `kw`, or use `J` inside an integral (see Composition below), one may use the following lowercase `integral` form to obtain the value of the integral directly (equivalent to building *and* calling `J(args...; kw...)`),
```julia
julia> integral(f, domain, args...; backend::AbstractBackend = default_backend(domain), result = nothing, kw...)
```

## Simple examples
### 1D integral
```julia
julia> using QuadGK
[ Info: Precompiling IntegrationInterfaceQuadGKExt [6a486dfe-6a5b-5d49-a9f9-02f4245ab8d6]

julia> integral(cos, Domain.interval(0, π/2))                      # Defaults to Backend.QuadGK()
1.0
```

### 2D integral
```julia
julia> using HCubature
[ Info: Precompiling IntegrationInterfaceQuadGKExt [6a486dfe-6a5b-5d49-a9f9-02f4245ab8d6]
[ Info: Precompiling IntegrationInterfaceHCubatureExt [b56c1907-5b79-5f70-9c8b-cc15ee23dcd1]

julia> integral((x,y) -> cos(x-y), Domain.Box((-1,-1), (1,1)))  # Defaults to Backend.HCubature()
2.8322936730937722
```

## Backends

As shown above, the integration is actually performed by backend packages that may be loaded as needed. Currently supported packages (weak dependencies) and corresponding backends in `IntegralBackends` are

- QuadGK.jl: `Backend.QuadGK(; opts...)` (calls `quadgk` and `quadgk!`, default for `Domain.Box{1}` domains)
- HCubature.jl:  `Backend.HCubature(; opts...)` (calls `hcubature`, default for `Domain.Box{N}` domains with `N ≠ 1`)
- Cubature.jl: `Backend.Cubature(; opts...)` (calls `hcubature`)

We also provide a `Backend.Quadrature((nodes, weights))` backend that can be used with the FastGaussQuadrature.jl package that computes nodes and weights for a 1D integral in the [-1, 1] integration domain. Nodes and weights are then scaled appropriately to the domain provided
```julia
julia> using FastGaussQuadrature, QuadGK

julia> f(x) = cos(x);

julia> J1 = Integral(f, Domain.interval(3,5); backend = Backend.Quadrature(gausslegendre(10)))
Integral
  Mutating   : false
  Domain     : Box{1}(3, 5)
  Backend    : Quadrature
  Integrand  : f

julia> J2 = Integral(f, Domain.interval(3,5); backend = Backend.QuadGK())
Integral
  Mutating   : false
  Domain     : Box{1}(3, 5)
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
  Domain     : Box{1}(2, 3)
  Backend    : QuadGK
  Integrand  : Integral
    Mutating   : false
    Domain     : Box{1}(0, 1)
    Backend    : QuadGK
    Integrand  : f

julia> J2 = Integral(f, Domain.interval((0,1), (2,3)); backend = Backend.HCubature())
Integral
  Mutating   : false
  Domain     : Box((0, 2), (1, 3))
  Backend    : HCubature
  Integrand  : f

julia> (J1(), J2())
(-3.800374064781164, -3.800374064812097)
```
Note the currying syntax used above for `J1`. It is equivalent to `J1 = Integral(Integral(f, Domain.Box{1}(0,1)), Domain.Box{1}(2,3))`.

## Functional domains
We can also make domains depend on `args` and `kw`s. In the above, this allows the inner domain to depend on outer integration variables. This provides one way to integrate over non-Box domains. For example,

$$J = \int_{-1}^1 dy\int_{-\sqrt{1-y^2}}^{\sqrt{1-y^2}} dx (x-y)^2\cos(x+y) $$

can be expressed as
```julia
julia> J = f |> Integral(Domain.Box{1}(y -> (-sqrt(1-y^2), sqrt(1-y^2)))) |> Integral(Domain.Box{1}(-1,1))
Integral
  Mutating   : false
  Domain     : Box{1}(-1, 1)
  Backend    : QuadGK
  Integrand  : Integral
    Mutating   : false
    Domain     : Functional{Box{1}}
    Backend    : QuadGK
    Integrand  : f

julia> J()
1.324825188363749

```
where we have passed a function to the `Domain.Box{1}` instead of the box bounds.

## Infinity

Some backends support using `Inf` to express unbounded domains. If that is not supported, we provide `Infinity(point::Number)` that can be used in box bounds instead. It represents an unbounded ray passing though `point` (which may be Real or not). These `Infinite` bounds are dealt with using an appropriate change of variables that takes `point` into account. As an example, consider a 2D half-plane `D = Domain.Box((; σ = 1) -> ((0, -Infinity(σ)), (Infinity(σ), Infinity(σ))))`. Note that it is a `Domain.Functional` object that depends on a keyword argument `σ`. We can integrate a Gaussian over `D`, which gives `π` for `σ = 1`

```julia
julia> f(x, y; σ = 1) = exp(-0.5*(x^2+y^2)/σ^2);

julia> D = Domain.Box((; σ = 1) -> ((0, -Infinity(σ)), (Infinity(σ), Infinity(σ))));

julia> J = Integral(f, D)
Integral
  Mutating   : false
  Domain     : Functional{Box}
  Backend    : HCubature
  Integrand  : f

julia> J(; σ = 1)
3.141592652846476
```
