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
        bb = IntegratedWoodSaxon(4, 0.6, 5.0, 1.0, 0.0, true, 1.0e-3, rng, 10)
        @test size(bb) == (4, 2)
        @test eltype(bb) == Float64

        sample = rand(rng, bb)
        @test size(sample) == (4, 2)
        @test eltype(sample) == Float64
    end

    @testset "IntegratedWoodSaxonInterp basics" begin
        cc = IntegratedWoodSaxonInterp(4, 0.6, 5.0, 1.0, 0.0, 10, 1.0e-3, 10, rng, true)
        @test size(cc) == (4, 2)
        @test eltype(cc) == Float64

        sample = rand(rng, cc)
        @test size(sample) == (4, 2)
        @test eltype(sample) == Float64
    end

    @testset "Threaded wrapper" begin
        base = IntegratedWoodSaxon(4, 0.6, 5.0, 1.0, 0.0, true, 1.0e-3, rng, 10)
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

    @testset "Participants sampling" begin
        n1 = NucleiWoodSaxon3D(rng, 4, 0.6, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10, 0.0)
        n2 = NucleiWoodSaxon3D(rng, 4, 0.6, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10, 0.0)

        participants = Participants(n1, n2, 0.5, 200.0, 1.0, 0.0; Nr = 8, Nth = 8)
        event = rand(rng, participants)

        @test event isa Participant
        @test event.n_coll > 0
        @test impactParameter(event) == event.b
        @test multiplicity(event) == event.multiplicity
        @test length(event.part1) <= n1.N_nucleon
        @test length(event.part2) <= n2.N_nucleon

        bmax = 3 * (n1.R + n2.R) + 6 * 0.5
        @test 0.0 <= event.b <= bmax
    end

    @testset "Centralities selection" begin
        function fake_participant(mult)
            part1 = [SVector(0.0, 0.0)]
            part2 = [SVector(0.0, 0.0)]
            shape1 = [1.0]
            shape2 = [1.0]
            return Participant(part1, part2, shape1, shape2, 1, 0.5, 1.0, 0.0, 1.0, 1.0, 0.0, mult)
        end

        events = [fake_participant(mult) for mult in 100.0:-1.0:1.0]
        mult = multiplicity.(events)
        event_perm = sortperm(mult, rev = true)
        mult_sorted = mult[event_perm]
        n_event_per_bin = length(events) ÷ 100
        borders = mult_sorted[1:n_event_per_bin:n_event_per_bin * 100]

        @test length(borders) == 100
        @test borders[1] == 100.0
        @test borders[end] == 1.0
        @test issorted(borders, rev = true)
    end

    @testset "InverseFunction" begin
        invf = InverseFunction(x -> 2 * x)
        @test isapprox(invf(4.0), 2.0; atol = 1.0e-6)
    end
end
