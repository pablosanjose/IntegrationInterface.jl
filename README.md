# IntegrationInterface

[![Build Status](https://github.com/pablosanjose/IntegrationInterface.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/pablosanjose/IntegrationInterface.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package aims to be a lightweight, faster-loading alternative to the excellent [Integrals.jl](https://github.com/SciML/Integrals.jl) from the SciML ecosystem. IntegrationInterface.jl offers an interface to perform n-dimensional numerical integrals of scalars, arrays or more general objects, over a domain $D$.

$$J(\text{args}...; \text{params}...) = \int_{D(\text{args}...; \text{params}...)} d^n x f(\boldsymbol x..., \text{args}...; \text{params}...)$$

The general interface reads
```julia
julia> J = integral(f, domain; solver::AbstractBackend = default_solver(domain), result = missing)
```
This produces an `J::Integral` object. Possible domains are produced with `Domain.Line` or`Domain.Box`. Here `Backend` and `Domain` are exported submodules of `IntegrationInterface`. Functions of `args` can be passed to the constructor of a domain to make it depend on arguments passed to `J`.

Mutating functions `f!(out, x..., args...; params...)` that modify an `out::A` in-place can also be used. This is useful for heap-allocated integrands of type e.g. `A::AbstractArray`. In this pass an array of type `A` with the `result` keyword.

To compute the integral for a set of `args` and `params`, use the syntax
```julia
julia> J(args...; params...)
```
In the mutating case, this will also write the value of the integral into `result`.

Example:
```julia
julia> using HCubature
[ Info: Precompiling IntegrationInterfaceQuadGKExt [6a486dfe-6a5b-5d49-a9f9-02f4245ab8d6]
[ Info: Precompiling IntegrationInterfaceHCubatureExt [b56c1907-5b79-5f70-9c8b-cc15ee23dcd1]

julia> f(x,y) = cos(x-y);

julia> J = integral(f, Domain.Box((-1,-1), (1,1)); solver = Backend.HCubature())
Integral
  Mutating   : false
  Domain     : Box((-1, -1), (1, 1))
  Solver     : HCubature
  Integrand  : f

julia> J()
2.8322936730937722
```

As shown above, the integration is actually performed by backend packages that may be loaded as needed. Currently supported packages (weak dependencies) and corresponding solvers in `IntegralSolvers` are

- QuadGK.jl: `Backend.QuadGK(; opts...)` (calls `quadgk` and `quadgk!`, default for `Domain.Line` domains)
- HCubature.jl:  `Backend.HCubature(; opts...)` (calls `hcubature`, default for `Domain.Box` domains)
- Cubature.jl: `Backend.Cubature(; opts...)` (calls `hcubature`)

We also provide a `Backend.Quadrature((nodes, weights))` solver that can be used with the FastGaussQuadrature.jl package that computes nodes and weights for a 1D integral in the [-1, 1] integration domain. Nodes and weights are then scaled appropriately to the domain provided
```julia
julia> using FastGaussQuadrature, QuadGK

julia> f(x) = cos(x);

julia> J1 = integral(f, Domain.Line(3,5); solver = Backend.Quadrature(gausslegendre(10)))
Integral
  Mutating   : false
  Domain     : Line(3, 5)
  Solver     : Quadrature
  Integrand  : f

julia> J2 = integral(f, Domain.Line(3,5); solver = Backend.QuadGK())
Integral
  Mutating   : false
  Domain     : Line(3, 5)
  Solver     : QuadGK
  Integrand  : f

julia> J1(), J2()
(-1.1000442827230061, -1.1000442827230057)
```

We can also express nested integrals such as

$$J(\text{args}...; \text{params}...) = \int_{D_n} dx_n\dots\int_{D_1} dx_1 f(\boldsymbol{x}, \text{args}...; \text{params}...)$$

As a concrete example, consider

$$J = \int_2^3 dy\int_0^1 dx (x-y)^2\cos(x+y) $$

This integral can be evaluated either as one adaptive `HCubature` or two nested `QuadGK` (the default solver):
```julia
julia> f(x,y) = (x-y)^2 * cos(x+y)
f (generic function with 2 methods)

julia> J1 = f |> integral(Domain.Line(0,1)) |> integral(Domain.Line(2,3))
Integral
  Mutating   : false
  Domain     : Line(2, 3)
  Solver     : QuadGK
  Integrand  : Integral
    Mutating   : false
    Domain     : Line(0, 1)
    Solver     : QuadGK
    Integrand  : f

julia> J2 = integral(f, Domain.Box((0,2), (1,3)); solver = Backend.HCubature())
Integral
  Mutating   : false
  Domain     : Box((0, 2), (1, 3))
  Solver     : HCubature
  Integrand  : f

julia> (J1(), J2())
(-3.800374064781164, -3.800374064812097)
```
Note the currying syntax used above for `J1`. It is equivalent to `J1 = integral(integral(f, Domain.Line(0,1)), Domain.Line(2,3))`.

We can also make inner domains depend on outer integration variables. For example,

$$J = \int_{-1}^1 dy\int_{-\sqrt{1-y^2}}^{\sqrt{1-y^2}} dx (x-y)^2\cos(x+y) $$

can be expressed as
```julia
julia> J = f |> integral(Domain.Line(y -> (-sqrt(1-y^2), sqrt(1-y^2)))) |> integral(Domain.Line(-1,1))
Integral
  Mutating   : false
  Domain     : Line(-1, 1)
  Solver     : QuadGK
  Integrand  : Integral
    Mutating   : false
    Domain     : Functional{Line}
    Solver     : QuadGK
    Integrand  : f

julia> J()
1.324825188363749

```
where we have passed a function to the `Domain.Line` instead of the line nodes.

Some backends support using `Inf` to express unbounded domains. If that is not supported, we provide `Infinity(point::Number)` that can be used instead. It represents an unbounded ray passing though `point` (which may be Real or not). These `Infinite` bounds are dealt with using an appropriate change of variables that takes `point` into account. As an example, consider a 2D half-plane `D = Domain.Box((; σ = 1) -> ((0, -Infinity(σ)), (Infinity(σ), Infinity(σ))))`. Note that it is a `Domain.Functional` object that depends on `σ`. We can integrate a Gaussian over `D`, which gives `π` for `σ = 1`

```julia
julia> f(x, y; σ = 1) = exp(-0.5*(x^2+y^2)/σ^2);

julia> J = integral(f, D)
Integral
  Mutating   : false
  Domain     : Functional{Box}
  Solver     : HCubature
  Integrand  : f

julia> J(; σ = 1)
3.141592652846476
```
