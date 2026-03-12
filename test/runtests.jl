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
    J = integral(f, Domain.Segment(0,π/2))
    @test J() ≈ 1
    @test J(2; σ = 1) ≈ π/2
    @test J(1; σ = 0, λ = 1) ≈ 0.5*(1+exp(-π/2))
    #default
    @test II.solver(J) isa Backend.QuadGK
    # Infs
    J = integral(f, Domain.Segment(0,Inf))
    @test J(1; σ = 0, λ = 1) ≈ 0.5
    J = integral(f, Domain.Segment(-Inf,Inf))
    @test J(1; σ = 0, λ = 1) ≈ 1.0
    # bounds as args
    g(x, a, b) = exp(-0.5*x^2)
    J = integral(g, Domain.Segment((a,b) -> (a,b)))
    @test J(0, Inf) ≈ -J(0, -Inf) ≈ - J(Inf, 0) ≈  J(-Inf, 0) ≈ √(π/2)
    @test J(-Inf, Inf) ≈ -J(Inf, -Inf) ≈ √(2π)
    # complex version
    g(x, a, b) = exp(x)
    @test J(1+2im, 3.0+3im) ≈  g(3.0+3im, 0, 0) - g(1+2im, 0, 0)
    # in-place version
    result = zeros(Float64, 10)
    inds = eachindex(result)
    g!(out, x, _...) = (out .= exp.(x .* inds))
    J = integral(g!, Domain.Segment((a,b) -> (a,b)); result)
    @test J(0, 1) ≈ (exp.(inds) .- 1) ./ inds
end

@testset "HCubature" begin
    f(x, y) = cos(x+y)
    f(x, y, p; σ = 2, λ = 0) = (cos(x+p*y)+λ)*exp(-(x^2+(p*y)^2)/(2*σ^2))
    J = integral(f, Domain.Box((0,0),(π/2,π)))
    @test J() ≈ -2
    @test J(2; λ = 2, σ = 1) ≈ 1.4803924093
    #default
    @test II.solver(J) isa Backend.HCubature
end

@testset "Cubature" begin
    f(x, y) = cos(x+y)
    f(x, y, p; σ = 2, λ = 0) = (cos(x+p*y)+λ)*exp(-(x^2+(p*y)^2)/(2*σ^2))
    J = integral(f, Domain.Box((0,0), (π/2,π)); solver = Backend.Cubature())
    @test J() ≈ -2
    @test J(2; λ = 2, σ = 1) ≈ 1.4803924093
    @test II.solver(J) isa Backend.Cubature
end
