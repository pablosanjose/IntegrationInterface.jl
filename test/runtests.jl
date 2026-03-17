using IntegrationInterface
using Test
using QuadGK
using HCubature
using Cubature

const II = IntegrationInterface

@testset "QuadGK" begin
    # basics
    f(x) = cos(x)
    f(x, p; σ = 2, λ = 0) = (σ + cos(p*x))*exp(-abs(x)*λ)
    J = integral(f, Domain.Box1D(0,π/2))
    @test II.backend(J) isa Backend.QuadGK     #default
    @test J() ≈ 1
    @test J(2; σ = 1) ≈ π/2
    @test J(1; σ = 0, λ = 1) ≈ 0.5*(1+exp(-π/2))

    # Infs
    J = integral(f, Domain.Box1D(0, Inf))
    @test J(1; σ = 0, λ = 1) ≈ 0.5
    J = integral(f, Domain.Box1D(-Inf, Inf))
    @test J(1; σ = 0, λ = 1) ≈ 1.0
    J = integral(f, Domain.Box1D(0,Infinity(1)))
    @test J(1; σ = 0, λ = 1) ≈ 0.5
    J = integral(f, Domain.Box1D(-Infinity(1), Infinity(1)))
    @test J(1; σ = 0, λ = 1) ≈ 1.0
     # mixing Inf with Infinity is ambiguous
    @test_throws ArgumentError integral(f, Domain.Box1D(-Infinity(1), Inf))

    # bounds as args
    g(x, _...) = exp(-0.5*x^2)
    J = integral(g, Domain.Box1D((a,b) -> (a,b)))
    @test J(0, Inf) ≈ -J(0, -Inf) ≈ - J(Inf, 0) ≈  J(-Inf, 0) ≈ √(π/2)
    @test J(-Inf, Inf) ≈ -J(Inf, -Inf) ≈ J(-Infinity(1), Infinity(1)) ≈ √(2π)

    # complex version
    h(x, _...) = exp(x)
    J = integral(h, Domain.Box1D((a,b) -> (a,b)))
    @test J(1+2im, 3.0+3im) ≈  h(3.0+3im, 0, 0) - h(1+2im, 0, 0)

    # in-place version
    result = zeros(Float64, 10)
    g!(out, x, _...) = (out .= exp.(x .* eachindex(out)))
    J = integral(g!, Domain.Box1D((a,b) -> (a,b)); result)
    @test J(0, 1) === result ≈ (exp.(eachindex(result)) .- 1) ./ eachindex(result)
end

@testset "HCubature" begin
    f(x, y) = cos(x+y)
    f(x, y, p; σ = 2, λ = 0) = (cos(x+p*y)+λ)*exp(-(x^2+(p*y)^2)/(2*σ^2))
    J = integral(f, Domain.Box((0,0),(π/2,π)))
    @test II.backend(J) isa Backend.HCubature   # default
    @test J() ≈ -2
    @test J(2; λ = 2, σ = 1) ≈ 1.4803924093

    # bounds as args
    g(x, y, _...) = cos(x+y)
    J = integral(g, Domain.Box((a,b) -> ((a,a), (b,b))))
    @test J(-1, 2) ≈ 4*cos(1)*sin(1.5)^2

    # complex domains (auto-converted to real)
    @test J(-1-im, 1+im) ≈ 4*sin(1 + im)^2

    # Infs
    g(x, y, _...) = exp(-0.5*(x^2+y^2))
    J = integral(g, Domain.Box((a,b) -> ((a,a), (b,b))))
    @test_throws ArgumentError J(-Inf, Inf)     # HCubature doesn't like unbounded domains
    @test J(-Infinity(1), Infinity(1)) ≈ 2π     # HCubature doesn't like unbounded domains

     # in-place version
    result = zeros(Float64, 10)
    g!(out, x, y, _...) = (out .= exp.((x+y) .* eachindex(out)))
    J = integral(g!, Domain.Box((a,b) -> ((a,a), (b,b))); result)
    @test_throws ArgumentError J(0, 1)          # HCubature doesn't do in-place
end

@testset "Cubature" begin
    f(x, y) = cos(x+y)
    f(x, y, p; σ = 2, λ = 0) = (cos(x+p*y)+λ)*exp(-(x^2+(p*y)^2)/(2*σ^2))
    J = integral(f, Domain.Box((0,0),(π/2,π)); backend = Backend.Cubature())
    @test II.backend(J) isa Backend.Cubature   # default
    @test J() ≈ -2
    @test J(2; λ = 2, σ = 1) ≈ 1.4803924093

    # bounds as args
    g(x, y, _...) = cos(x+y)
    J = integral(g, Domain.Box((a,b) -> ((a,a), (b,b))); backend = Backend.Cubature())
    @test J(-1, 2) ≈ 4*cos(1)*sin(1.5)^2

    # complex domains (auto-converted to real)
    @test_throws ArgumentError J(-1-im, 1+im)   # Cubature doesn't support complex functions

    g!(out, x, y, _...) = (out .= cos(x+y))
    result = [0.0+0im]
    J = integral(g!, Domain.Box((a,b) -> ((a,a), (b,b))); result, backend = Backend.Cubature())
    @test J(-1-im, 1+im) === result ≈ [4*sin(1 + im)^2]

    # Infs
    g(x, y, _...) = exp(-0.5*(x^2+y^2))
    J = integral(g, Domain.Box((a,b) -> ((a,a), (b,b))); backend = Backend.Cubature())
    @test_throws ArgumentError J(-Inf, Inf)     # HCubature doesn't like unbounded domains
    @test J(-Infinity(1), Infinity(1)) ≈ 2π     # HCubature doesn't like unbounded domains

    # in-place version
    result = zeros(Float64, 10)
    g!(out, x, y, _...) = (out .= exp.(-(x+y) .* eachindex(out)))
    J = integral(g!, Domain.Box((a,b) -> ((a,a), (b,b))); result, backend = Backend.Cubature())
    @test J(0, 1) === result ≈ ((exp.( .- eachindex(result)) .- 1) ./ eachindex(result)) .^ 2
end

@testset "Nested" begin
    # 2D
    f(x, y) = cos(x^2+y)*exp(-0.5*(x^2+2y^2))
    J1 = f |> integral(Domain.Box1D(-Inf, Inf)) |> integral(Domain.Box1D(-Infinity(1),Infinity(1)))
    J2 = f |> integral(Domain.Box((-Infinity(1), -Infinity(1)),(Infinity(1), Infinity(1))))
    @test J1() ≈ J2()
    J1 = f |> integral(Domain.Box1D(0, 20)) |> integral(Domain.Box1D(-Infinity(1),Infinity(5)))
    J2 = f |> integral(Domain.Box((0, -Infinity(1)),(20, Infinity(1))))
    @test J1() ≈ J2()
    # 3D
    f(x, y, z) = cos(x+2y+3z)
    J1 = f |> integral(Domain.Box1D(-1, 1)) |> integral(Domain.Box((0, 0),(1, 2+im)))
    J2 = f |> integral(Domain.Box((-1, 0),(1, 1))) |> integral(Domain.Box1D(0, 2+im))
    J3 = f |> integral(Domain.Box1D(-1, 1)) |> integral(Domain.Box1D(0, 1)) |> integral(Domain.Box1D(0, 2+im))
    J4 = f |> integral(Domain.Box((-1, 0, 0),(1, 1, 2+im)))
    @test J1() ≈ J2() ≈ J3() ≈ J4()
    f(x, y, z) = exp(-abs(x+2y+3z))
    J1 = f |> integral(Domain.Box1D(-1, 1)) |> integral(Domain.Box((0, -Infinity(2+3im)),(1, 2+im)))
    J2 = f |> integral(Domain.Box((-1, 0),(1, 1))) |> integral(Domain.Box1D(-Infinity(2+3im), 2+im))
    J3 = f |> integral(Domain.Box1D(-1, 1)) |> integral(Domain.Box1D(0, 1)) |> integral(Domain.Box1D(-Infinity(2+3im), 2+im))
    J4 = f |> integral(Domain.Box((-1, 0, -Infinity(2+3im)),(1, 1, 2+im)))
    @test J1() ≈ J2() ≈ J3() ≈ J4()
end


@testset "Domain sums" begin
    f(x) = 1/(x^2+1)
    J = integral(f, Domain.Box1D(0, 1, 1+2im, -1+2im, -1, 0))
    @test J() ≈ π
    f(x) = 1/(x-im)
    @test J() ≈ 2π*im
    J = integral(f, Domain.Box1D([0, 1, 1+2im, -1+2im, -1, 0]); backend = Backend.QuadGK())
    @test J() ≈ 2π*im
end
