"""
`Domain` is a submodule of `IntegrationInterface` that bundles possible integration domains.
Currently it provides
    - `Domain.Box{N}`: N-dimensional hypercube domain, created with `Domain.box`, that may be bounded or unbounded.
"""
Domain

"""
    Domain.Box1D(a::Union{Number,Infinity}, b::Union{Number,Infinity})

Create a 1D integration domain `Domain.Box{1}` from `a` to `b`.

    Domain.Box1D(xs::Union{Number,Infinity}...)

Create a sum of `Domain.Box{1}` subdomains corresponding to intervals `(xs[1], xs[2])`,
`(xs[2], xs[3])`, `(xs[3], xs[4])`,... . For real-valued variables, or for complex variables
of holomorphic functions, this is equivalent to an integral from `first(xs)` to `last(xs)`.

# Examples
```julia
julia> using QuadGK

julia> J = integral(cos, Domain.Box1D(0, π/2))
Integral
  Mutating   : false
  Domain     : Box{1}(0, 1.5707963267948966)
  Backend    : QuadGK
  Integrand  : cos

julia> J()
1.0
```
# See also:
    `Box`, `Infinity`, `integral`

"""
Domain.Box1D

"""
    Domain.Box((x₁ᵐⁱⁿ, x₂ᵐⁱⁿ, ..., xₙᵐⁱⁿ), (x₁ᵐᵃˣ, x₂ᵐᵃˣ, ..., xₙᵐᵃˣ))

Create an integration domain `Domain.Box{N}` for a function `f(x₁, x₂, ..., xₙ)` over an `N`
dimensional hypercube defined by the intervals `(xᵢᵐⁱⁿ, xᵢᵐᵃˣ)`.

# Examples
```julia
julia> using HCubature

julia> J = integral((x,y,z) -> cos(x+y+z), Domain.Box((0, 0, 0), (π/2, π/2, π/2)))
Integral
  Mutating   : false
  Domain     : Box{3}((0, 0, 0), (1.5707963267948966, 1.5707963267948966, 1.5707963267948966)))
  Backend    : HCubature
  Integrand  : #1

julia> J()
-1.9999999998615692
```

# See also:
    `Box1D`, `Infinity`, `integral`

"""
Domain.Box

"""
    Infinity(x::Number)

Creates an `Infinity` object representing an unbounded 1D ray passing through `x`. It can be
used in place of `Inf` to create unbounded domains, e.g. `Domain.Box1D(0, Infinity(1+im))`,
which represents a straight semi-infinite line in the complex plane starting at 0 and
extending to infinity through point `1+im`.

# Example
```julia
julia> using HCubature

julia> J = integral((x, y; σ) -> 1/(2π*σ^2) * exp(-(x^2+y^2)/(2*σ^2)), Domain.Box((; σ) -> ((0,0), (Infinity(σ), Infinity(σ)))))
Integral
  Mutating   : false
  Domain     : Functional{Box}
  Backend    : HCubature
  Integrand  : #39

julia> J(σ = 1e5)
0.2499999998823345
```
"""
Infinity
