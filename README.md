# IntegrationInterface

[![Build Status](https://github.com/pablosanjose/IntegrationInterface.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/pablosanjose/IntegrationInterface.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package aims to be a lightweight, faster-loading alternative to the excellent [Integrals.jl](https://github.com/SciML/Integrals.jl) from the SciML ecosystem. IntegrationInterface.jl offers an interface to perform n-dimensional numerical integrals of scalars, arrays or more general objects, over a domain $D$.

$$J(\text{args}...; \text{params}...) = \int_D d^n x f(\boldsymbol x, \text{args}...; \text{params}...)$$

The general interface reads
```julia
julia> J = integral(f, domain; solver::AbstractBackend = Backends.QuadGK(), result = missing)
```
This produces an `J::Integral` object. Possible domains are produced with `Domains.Segment`, `Domains.Box` or tuples thereof (for nested integrals, see below). Here `Backends` and `Domains` are exported submodules of `IntegrationInterface`

Mutating functions `f!(out, x, args...; params...)` that modify an `out::A` in-place can also be used. This is useful for heap-allocated integrands of type e.g. `A::AbstractArray`. In this pass an array of type `A` with the `result` keyword.

To compute the integral for a set of `args` and `params`, use the syntax
```julia
julia> J(args...; params...)
```
In the mutating case, this will also write the value of the integral into `result`.

The integration is actually performed by backend packages that may be loaded as needed. Currently supported packages (weak dependencies) and corresponding solvers in `IntegralSolvers` are

- QuadGK.jl: `Backends.QuadGK(; opts...)` (calls `quadgk` and `quadgk!`, for 1D integrals over a `Domains.Segment`)
- Cubature.jl: `Backends.Cubature(; opts...)` (calls `hcubature`, for n-D integrals over a `Domains.Box`)
- HCubature.jl:  `Backends.HCubature(; opts...)` (calls `hcubature`, for n-D integrals in a `Domains.Box`)

We also provide a `Backends.Quadrature((nodes, weights))` solver that can be used with the FastGaussQuadrature.jl package that computes nodes and weights for a 1D integral in the [-1, 1] integration domain. Nodes and weights are then scaled appropriately to the domain provided
```julia
julia> using FastGaussQuadrature

julia> f(x) = cos(x);

julia> J1 = integral(f, Domains.Segment(3,5); solver = Backends.Quadrature(gausslegendre(10)))
Integral: callable object representing a numerical integral over a domain
  Mutating   : false
  Domain     : Segment(3, 5)
  Solver     : Quadrature

julia> J2 = integral(f, Domains.Segment(3,5); solver = Backends.QuadGK())
Integral: callable object representing a numerical integral over a domain
  Mutating   : false
  Domain     : Segment(3, 5)
  Solver     : QuadGK

julia> J1(), J2()
(-1.1000442827230061, -1.1000442827230057)
```

We can also express nested integrals using different solvers and domains for each of them,

$$J(\text{args}...; \text{params}...) = \int_{D_n} dx_n\dots\int_{D_1} dx_1 f(\boldsymbol{x}, \text{args}...; \text{params}...)$$

As a concrete example, consider

$$J = \int_2^3 dy\int_0^1 dx (x-y)^2\cos(x+y) $$

This integral can be evaluated as one adaptive `HCubature` or two nested `QuadGK`:
```julia
julia> f((x,y)) = (x-y)^2 * cos(x+y)
f (generic function with 1 method)

julia> J1 = integral(f, Domains.Segment(0,1), Domains.Segment(2,3); solver = (Backends.QuadGK(), Backends.QuadGK()))
Integral: callable object representing a numerical integral over a domain
  Mutating   : false
  Domain     : ([0, 1], [2, 3])
  Solver     : Multi(QuadGK, QuadGK)

julia> J2 = integral(f, Domains.Box((0,2), (1,3)); solver = Backends.HCubature())
Integral: callable object representing a numerical integral over a domain
  Mutating   : false
  Domain     : ([0, 2], [1, 3])
  Solver     : HCubature

julia> (J1(), J2())
(-3.800374064781164, -3.800374064812097)

```
Note that the syntax of the domain is different for two nested `Backends.QuadGK` (intervals for each integral, from inner to outermost) versus `Backends.HCubature` (min/max corners of the rectangle domain). When using `Multi`, `f(r, ...)` and `f!(out, r, ...)` must accept points `r` of type `Tuple`. Note that when using nested solvers through `Multi`, we can make the domain for solver `i` depend on the coordinates `(r[i+1],..., r[n])` of the outer integrals. We achieve this by passing a function as the domain. As an example, the integral of the above `f` on a unit circle,

$$J = \int_{-1}^1 dy\int_{-\sqrt{1-y^2}}^{\sqrt{1-y^2}} dx (x-y)^2\cos(x+y) $$

can be expressed as
```julia
julia> J = integral(f, y -> Domains.Segment(-sqrt(1-y^2), sqrt(1-y^2)), Domains.Segment(-1,1), solver = (Backends.QuadGK(), Backends.QuadGK()));

julia> J()
1.324825188363749

```
