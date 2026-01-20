using MonteCarloGlauber
using Random
using StaticArrays
using Test

@testset "MonteCarloGlauber.jl" begin
    rng = MersenneTwister(1234)

    @testset "NucleiWoodSaxon3D basics" begin
        aa = NucleiWoodSaxon3D(rng, 4, 0.6, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10, 0.0)
        @test size(aa) == (4, 2)
        @test eltype(aa) == Float64

        sample = rand(rng, aa)
        @test size(sample) == (4, 2)
        @test eltype(sample) == Float64
    end

    @testset "IntegratedWoodSaxon basics" begin
        bb = IntegratedWoodSaxon(4, 0.6, 5.0, 1.0, 0.0, true, 1e-3, rng, 10)
        @test size(bb) == (4, 2)
        @test eltype(bb) == Float64

        sample = rand(rng, bb)
        @test size(sample) == (4, 2)
        @test eltype(sample) == Float64
    end

    @testset "IntegratedWoodSaxonInterp basics" begin
        cc = IntegratedWoodSaxonInterp(4, 0.6, 5.0, 1.0, 0.0, 10, 1e-3, 10, rng, true)
        @test size(cc) == (4, 2)
        @test eltype(cc) == Float64

        sample = rand(rng, cc)
        @test size(sample) == (4, 2)
        @test eltype(sample) == Float64
    end

    @testset "Threaded wrapper" begin
        base = IntegratedWoodSaxon(4, 0.6, 5.0, 1.0, 0.0, true, 1e-3, rng, 10)
        dd = Threaded(base, 2)
        @test size(dd) == (4, 2)
        @test eltype(dd) == Float64

        sample = rand(rng, dd, 6)
        @test size(sample) == (6,)
        @test eltype(sample) == Matrix{Float64}
    end

    @testset "Participant helpers" begin
        part1 = [SVector(0.0, 0.0), SVector(0.5, -0.25)]
        part2 = [SVector(0.25, 0.25)]
        shape1 = [1.0, 1.0]
        shape2 = [1.0]
        participant = Participant(part1, part2, shape1, shape2, 1, 0.5, 1.0, 0.0, 5.0, 5.0, 0.0, 1.0)

        mult, x_cm, y_cm = center_of_mass(participant, 8, 8)
        @test mult > 0
        @test isfinite(x_cm)
        @test isfinite(y_cm)
    end

    @testset "InverseFunction" begin
        invf = InverseFunction(x -> 2 * x)
        @test isapprox(invf(4.0), 2.0; atol=1e-6)
    end
end
