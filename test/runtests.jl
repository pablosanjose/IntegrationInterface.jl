using IntegrationInterface
using Test
using QuadGK
using HCubature
using Cubature
using HAdaptiveIntegration
using FastGaussQuadrature
using Unitful

const II = IntegrationInterface

@testset begin "integral"
    f(x) = cos(x)
    f(x, p; σ = 2, λ = 0) = (σ + cos(p*x))*exp(-abs(x)*λ)
    @test integral(f, Domain.Box{1}(0,π/2)) ≈ 1
    @test integral(f, Domain.Box{1}(0,π/2), 2; σ = 1) ≈ π/2
    @test f |> integral(Domain.Box{1}(0,π/2), 2; σ = 1) ≈ π/2
    @test integral(x -> exp(-x), Domain.Box{1}(0, Infinity(1))) ≈ 1
    @test integral((x,y) -> exp(-x-y), Domain.Box((0, 0), (Infinity(1), Infinity(1)))) ≈ 1
    @test integral((x; σ = 1) -> exp(-x/σ), Domain.Box{1}((; σ=1) -> (0, Infinity(σ))); σ = 4) ≈ 4
    @test integral((x, a; σ = 1) -> a*exp(-x/σ), Domain.Box{1}((a; σ=1) -> (0, Infinity(σ))), 4; σ = 4) ≈ 16
end

@testset begin "Unitful"
    # 1D
    f(x) = √(1 - x^2) * u"A"
    J1 = Integral(f, Domain.Box{1}(-1,1))
    J2 = Integral(f, Domain.Box{1}(-1,1); backend = Backend.Quadrature(gausslegendre(50)))
    @test J1() ≈ π/2 * u"A"
    @test J2() ≈ π/2 * u"A" atol = 1e-4*u"A"
    f(x) = √(1u"A"^2 - x^2)
    J1 = Integral(f, Domain.Box{1}(-1u"A",1u"A"))
    J2 = Integral(f, Domain.Box{1}(-1u"A",1u"A"); backend = Backend.Quadrature(gausslegendre(50)))
    @test J1() ≈ π/2 * u"A"^2
    @test J2() ≈ π/2 * u"A"^2 atol = 1e-4*u"A"^2
    # 2D
    for backend in (Backend.Quadrature(gausslegendre(50)), Backend.HCubature(), Backend.Cubature(), Backend.HAdaptiveIntegration())
        J = Integral((x,y)->cos(ustrip(x+y)), Domain.Box((0u"A", 0u"A"), (π/2 * u"A", π * u"A")); backend)
        backend isa Backend.Quadrature ? (@test (J() ≈ -2u"A"^2)) : (@test_throws ArgumentError J())
        J = Integral((x,y)->x*y*u"A", Domain.Box((0,0),(1,1)); backend)
        backend isa Backend.Quadrature ? (@test (J() ≈ 0.25u"A")) : (@test_throws ArgumentError J())
    end
end

@testset "Quadrature" begin
    backend = Backend.Quadrature(gausslegendre(50))
    # basics 1D
    f(x) = cos(x)
    f(x, p; σ = 2, λ = 0) = (σ + cos(p*x))*exp(-abs(x)*λ)
    J = Integral(f, Domain.Box{1}(0,π/2); backend)
    @test II.backend(J) isa Backend.Quadrature     #default
    @test J() ≈ 1
    @test J(2; σ = 1) ≈ π/2
    @test J(1; σ = 0, λ = 1) ≈ 0.5*(1+exp(-π/2))

    # basics 2D
    f(x, y) = cos(x+y)
    f(x, y, p; σ = 2, λ = 0) = (cos(x+p*y)+λ)*exp(-(x^2+(p*y)^2)/(2*σ^2))
    J = Integral(f, Domain.Box((0,0),(π/2,π)); backend)
    @test II.backend(J) isa Backend.Quadrature   # default
    @test J() ≈ -2
    @test J(2; λ = 2, σ = 1) ≈ 1.4803924093 atol = 1e-6
    f(x,y) = exp(-x^2-y^2)
    J = Integral(f, Domain.Box((-Infinity(1),-Infinity(1)),(Infinity(1), Infinity(1))); backend)
    @test J() ≈ π

    # Infs
    J = Integral(f, Domain.Box{1}(0, Inf); backend)
    @test_throws ArgumentError J(1; σ = 0, λ = 1) ≈ 0.5
    J = Integral(f, Domain.Box{1}(-Inf, Inf); backend)
    @test_throws ArgumentError J(1; σ = 0, λ = 1) ≈ 1.0
    J = Integral(f, Domain.Box{1}(0,Infinity(1)); backend)
    @test J(1; σ = 0, λ = 1) ≈ 0.5 atol = 1e-3
    J = Integral(f, Domain.Box{1}(-Infinity(1), Infinity(1)); backend)
    @test J(1; σ = 0, λ = 1) ≈ 1.0 atol = 1e-3
     # mixing Inf with Infinity is ambiguous
    @test_throws ArgumentError Integral(f, Domain.Box{1}(-Infinity(1), Inf))

    # bounds as args
    g(x, _...) = exp(-0.5*x^2)
    J = Integral(g, Domain.Box{1}((a,b) -> (a,b)); backend)
    @test J(-Infinity(1), Infinity(1)) ≈ √(2π)

    # complex version
    h(x, _...) = exp(x)
    J = Integral(h, Domain.Box{1}((a,b) -> (a,b)); backend)
    @test J(1+2im, 3.0+3im) ≈  h(3.0+3im, 0, 0) - h(1+2im, 0, 0)

    # in-place version
    result = zeros(Float64, 10)
    g!(out, x, _...) = (out .= exp.(x .* eachindex(out)))
    J = Integral(g!, Domain.Box{1}((a,b) -> (a,b)); result, backend)
    @test J(0, 1) === result ≈ (exp.(eachindex(result)) .- 1) ./ eachindex(result)
end

@testset "QuadGK" begin
    # basics
    f(x) = cos(x)
    f(x, p; σ = 2, λ = 0) = (σ + cos(p*x))*exp(-abs(x)*λ)
    J = Integral(f, Domain.Box{1}(0,π/2))
    @test II.backend(J) isa Backend.QuadGK     #default
    @test J() ≈ 1
    @test J(2; σ = 1) ≈ π/2
    @test J(1; σ = 0, λ = 1) ≈ 0.5*(1+exp(-π/2))

    # Infs
    J = Integral(f, Domain.Box{1}(0, Inf))
    @test J(1; σ = 0, λ = 1) ≈ 0.5
    J = Integral(f, Domain.Box{1}(-Inf, Inf))
    @test J(1; σ = 0, λ = 1) ≈ 1.0
    J = Integral(f, Domain.Box{1}(0,Infinity(1)))
    @test J(1; σ = 0, λ = 1) ≈ 0.5
    J = Integral(f, Domain.Box{1}(-Infinity(1), Infinity(1)))
    @test J(1; σ = 0, λ = 1) ≈ 1.0
     # mixing Inf with Infinity is ambiguous
    @test_throws ArgumentError Integral(f, Domain.Box{1}(-Infinity(1), Inf))

    # bounds as args
    g(x, _...) = exp(-0.5*x^2)
    J = Integral(g, Domain.Box{1}((a,b) -> (a,b)))
    @test J(0, Inf) ≈ -J(0, -Inf) ≈ - J(Inf, 0) ≈  J(-Inf, 0) ≈ √(π/2)
    @test J(-Inf, Inf) ≈ -J(Inf, -Inf) ≈ J(-Infinity(1), Infinity(1)) ≈ √(2π)

    # complex version
    h(x, _...) = exp(x)
    J = Integral(h, Domain.Box{1}((a,b) -> (a,b)))
    @test J(1+2im, 3.0+3im) ≈  h(3.0+3im, 0, 0) - h(1+2im, 0, 0)

    # in-place version
    result = zeros(Float64, 10)
    g!(out, x, _...) = (out .= exp.(x .* eachindex(out)))
    J = Integral(g!, Domain.Box{1}((a,b) -> (a,b)); result)
    @test J(0, 1) === result ≈ (exp.(eachindex(result)) .- 1) ./ eachindex(result)
end

@testset "HCubature" begin
    f(x, y) = cos(x+y)
    f(x, y, p; σ = 2, λ = 0) = (cos(x+p*y)+λ)*exp(-(x^2+(p*y)^2)/(2*σ^2))
    J = Integral(f, Domain.Box((0,0),(π/2,π)))
    @test II.backend(J) isa Backend.HCubature   # default
    @test J() ≈ -2
    @test J(2; λ = 2, σ = 1) ≈ 1.4803924093

    # bounds as args
    g(x, y, _...) = cos(x+y)
    J = Integral(g, Domain.Box{2}((a,b) -> ((a,a), (b,b))))
    @test J(-1, 2) ≈ 4*cos(1)*sin(1.5)^2

    # complex domains (auto-converted to real)
    @test J(-1-im, 1+im) ≈ 4*sin(1 + im)^2

    # Infs
    g(x, y, _...) = exp(-0.5*(x^2+y^2))
    J = Integral(g, Domain.Box{2}((a,b) -> ((a,a), (b,b))))
    @test_throws ArgumentError J(-Inf, Inf)     # HCubature doesn't like unbounded domains
    @test J(-Infinity(1), Infinity(1)) ≈ 2π     # HCubature doesn't like unbounded domains

     # in-place version
    result = zeros(Float64, 10)
    g!(out, x, y, _...) = (out .= exp.((x+y) .* eachindex(out)))
    J = Integral(g!, Domain.Box{2}((a,b) -> ((a,a), (b,b))); result)
    @test_throws ArgumentError J(0, 1)          # HCubature doesn't do in-place
end

@testset "Cubature" begin
    f(x, y) = cos(x+y)
    f(x, y, p; σ = 2, λ = 0) = (cos(x+p*y)+λ)*exp(-(x^2+(p*y)^2)/(2*σ^2))
    J = Integral(f, Domain.Box((0,0),(π/2,π)); backend = Backend.Cubature())
    @test II.backend(J) isa Backend.Cubature   # default
    @test J() ≈ -2
    @test J(2; λ = 2, σ = 1) ≈ 1.4803924093

    # bounds as args
    g(x, y, _...) = cos(x+y)
    J = Integral(g, Domain.Box{2}((a,b) -> ((a,a), (b,b))); backend = Backend.Cubature())
    @test J(-1, 2) ≈ 4*cos(1)*sin(1.5)^2

    # complex domains (auto-converted to real)
    @test_throws ArgumentError J(-1-im, 1+im)   # Cubature doesn't support complex functions

    g!(out, x, y, _...) = (out .= cos(x+y))
    result = [0.0+0im]
    J = Integral(g!, Domain.Box{2}((a,b) -> ((a,a), (b,b))); result, backend = Backend.Cubature())
    @test J(-1-im, 1+im) === result ≈ [4*sin(1 + im)^2]

    # Infs
    g(x, y, _...) = exp(-0.5*(x^2+y^2))
    J = Integral(g, Domain.Box{2}((a,b) -> ((a,a), (b,b))); backend = Backend.Cubature())
    @test_throws ArgumentError J(-Inf, Inf)     # HCubature doesn't like unbounded domains
    @test J(-Infinity(1), Infinity(1)) ≈ 2π     # HCubature doesn't like unbounded domains

    # in-place version
    result = zeros(Float64, 10)
    g!(out, x, y, _...) = (out .= exp.(-(x+y) .* eachindex(out)))
    J = Integral(g!, Domain.Box{2}((a,b) -> ((a,a), (b,b))); result, backend = Backend.Cubature())
    @test J(0, 1) === result ≈ ((exp.( .- eachindex(result)) .- 1) ./ eachindex(result)) .^ 2
end

@testset "HAdaptiveIntegration" begin
    f(x, y; σ = 1, _...) = cis(σ * x * y) * exp(-x^2-y)
    J = Integral(f, Domain.Simplex{2}((; a = 1, b = 1, _...) -> ((0, 0), (a, 0), (0, b))))
    @test II.backend(J) isa Backend.HAdaptiveIntegration   # default
    @test J() ≈ 0.309567883488632 + 0.0230970240588801im
    @test J(; σ = 2) ≈ 0.305201476947410 + 0.0457369920083456im
    @test J(; a = 3, b = 4) ≈ 0.685318616096836 + 0.291629495280167im

    J = Integral(f, Domain.Simplex((0, 0), (1, 0), Infinity(0, 2)))
    @test J() ≈ 0.618821963308144 + 0.231710996331552im

    J = Integral(f, Domain.Simplex((0,0), (1, 1/2), (1/2, 1)))
    @test J() ≈ 0.174223048821802 + 0.038391266900971im

    J = Integral(f, Domain.Box{2}((; a = 1, b = 1, _...) -> ((0, 0), (a, b))); backend = Backend.HAdaptiveIntegration())
    @test II.backend(J) isa Backend.HAdaptiveIntegration
    @test J() ≈ 0.309567883488632 + 0.0230970240588801im
    @test J(; σ = 2) ≈ 0.305201476947410 + 0.0457369920083456im
    @test J(; a = 3, b = 4) ≈ 0.685318616096836 + 0.291629495280167im
    @test J(; σ = im, b = Infinity(2)) ≈ 0.544256200588910
    @test_throws ArgumentError J(; b = Inf)
end

@testset "Nested" begin
    # 2D
    f(x, y) = cos(x^2+y)*exp(-0.5*(x^2+2y^2))
    J1 = f |> Integral(Domain.Box{1}(-Inf, Inf)) |> Integral(Domain.Box{1}(-Infinity(1),Infinity(1)))
    J2 = f |> Integral(Domain.Box((-Infinity(1), -Infinity(1)),(Infinity(1), Infinity(1))))
    @test J1() ≈ J2()
    J1 = f |> Integral(Domain.Box{1}(0, 20)) |> Integral(Domain.Box{1}(-Infinity(1),Infinity(5)))
    J2 = f |> Integral(Domain.Box((0, -Infinity(1)),(20, Infinity(1))))
    @test J1() ≈ J2()
    # 3D
    f(x, y, z) = cos(x+2y+3z)
    J1 = f |> Integral(Domain.Box{1}(-1, 1)) |> Integral(Domain.Box((0, 0),(1, 2+im)))
    J2 = f |> Integral(Domain.Box((-1, 0),(1, 1))) |> Integral(Domain.Box{1}(0, 2+im))
    J3 = f |> Integral(Domain.Box{1}(-1, 1)) |> Integral(Domain.Box{1}(0, 1)) |> Integral(Domain.Box{1}(0, 2+im))
    J4 = f |> Integral(Domain.Box((-1, 0, 0),(1, 1, 2+im)))
    @test J1() ≈ J2() ≈ J3() ≈ J4()
    f(x, y, z) = exp(-abs(x+2y+3z))
    J1 = f |> Integral(Domain.Box{1}(-1, 1)) |> Integral(Domain.Box((0, -Infinity(2+3im)),(1, 2+im)))
    J2 = f |> Integral(Domain.Box((-1, 0),(1, 1))) |> Integral(Domain.Box{1}(-Infinity(2+3im), 2+im))
    J3 = f |> Integral(Domain.Box{1}(-1, 1)) |> Integral(Domain.Box{1}(0, 1)) |> Integral(Domain.Box{1}(-Infinity(2+3im), 2+im))
    J4 = f |> Integral(Domain.Box((-1, 0, -Infinity(2+3im)),(1, 1, 2+im)))
    @test J1() ≈ J2() ≈ J3() ≈ J4()
end

@testset "Domain sums" begin
    f(x) = 1/(x^2+1)
    J = Integral(f, Domain.Box{1}(0, 1, 1+2im, -1+2im, -1, 0))
    @test J() ≈ π
    f(x) = 1/(x-im)
    @test J() ≈ 2π*im
    J = Integral(f, Domain.Box{1}([0, 1, 1+2im, -1+2im, -1, 0]); backend = Backend.QuadGK())
    @test J() ≈ 2π*im
end

@testset "type promotion" begin
    @test Domain.Box(Int16(1),Infinity(Int32(3))) isa Domain.Box{1,Float64}
    @test Domain.Box((1, 2.0f0), (3, Infinity(4))) isa Domain.Box{2,Float32}
    @test Domain.Box((1, 2.0f0), (3, Infinity(4))) isa Domain.Box{2,Float32}
    @test_throws ArgumentError Domain.Box((1u"A", 2.0f0), (3, Infinity(4)))
    J = Integral(cos, Domain.Box{1}(-Int32(1),Int16(1)); backend = Backend.Quadrature(gausslegendre(50)))
    @test J() isa Float64
    J = Integral(cos, Domain.Box{1}(-Float32(π/2),Float32(π/2)); backend = Backend.Quadrature(gausslegendre(50)))
    @test J() isa Float32
    J = Integral(cos, Domain.Box{1}(-Float32(π/2),Float32(π/2)); backend = Backend.QuadGK())
    @test J() isa Float32
    J = Integral((x,y) -> cos(x-y), Domain.Box((0f0,0f0),(-Float32(π/2),Float32(π))); backend = Backend.HCubature())
    @test J() isa Float32
    J = Integral((x,y) -> cos(x-y), Domain.Box((0f0,0f0),(-Float32(π/2),Float32(π))); backend = Backend.Cubature())
    @test_broken J() isa Float32  # Cubature does not respect types
end
