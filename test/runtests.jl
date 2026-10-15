using MonteCarloGlauber
using Random
using StaticArrays
using Statistics
using QuadGK
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

    @testset "Woods-Saxon radial density is sampled exactly" begin
        # Lead: R = 6.65, α = 0.54. Compare the sampled radial histogram with the analytic
        # Woods-Saxon profile shell by shell, and the mean radius. The previous Metropolis
        # sampler failed this test by ~10σ in the inner shells.
        R, α = 6.65, 0.54
        ρ(r) = 1 / (1 + exp((r - R) / α))
        Z = quadgk(r -> r^2 * ρ(r), 0, 30)[1]
        n = Lead()
        rng = MersenneTwister(42)
        N = 100_000
        rs = [hypot(rand(rng, n.sampler)...) for _ in 1:N]
        for lo in 0.0:1.0:9.0
            expected = N * quadgk(r -> r^2 * ρ(r), lo, lo + 1)[1] / Z
            observed = count(r -> lo <= r < lo + 1, rs)
            @test abs(observed - expected) < 4 * sqrt(expected)
        end
        mean_r_exact = quadgk(r -> r^3 * ρ(r), 0, 30)[1] / Z
        @test abs(mean(rs) - mean_r_exact) < 0.02
        # projected (2D) samplers: ⟨x²+y²⟩ = (2/3)⟨r²⟩
        r2_exact = quadgk(r -> r^4 * ρ(r), 0, 30)[1] / Z
        itp = IntegratedWoodSaxonInterp(208, α, R, 1.0, 0.0, 10, 1.0e-3, 100, rng, true)
        r2_2d = mean(sum(abs2, rand(rng, itp.sampler)) for _ in 1:N)
        @test abs(r2_2d / (2 / 3 * r2_exact) - 1) < 0.02
    end

    @testset "Rejection sampler envelope is a guaranteed bound" begin
        # Regression: a very diffuse profile (α = R) whose peak a coarse grid scan missed.
        n = NucleiWoodSaxon3D(10, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        @test n.sampler.bound >= 1 / (1 + exp(-1.0))              # true maximum of ρ
        ρ(r) = 1 / (1 + exp(r - 1))
        Z = quadgk(r -> r^2 * ρ(r), 0, 40)[1]
        trng = MersenneTwister(5)
        rs = [hypot(rand(trng, n.sampler)...) for _ in 1:200_000]
        @test abs(mean(rs) - quadgk(r -> r^3 * ρ(r), 0, 40)[1] / Z) < 0.012
        p_core = quadgk(r -> r^2 * ρ(r), 0, 0.5)[1] / Z
        @test abs(mean(rs .< 0.5) - p_core) < 3.5 * sqrt(p_core / length(rs))

        # bound ≥ target at random points in the box, for deformed, modulated and projected targets
        for nuc in (
                Uranium(),
                NucleiWoodSaxon3D(10, 0.5, 5.0, 2.0, 0.3, 0.2, 0.1, 0.1, 0.4),     # w > 0, triaxial
                IntegratedWoodSaxon(10, 0.6, 5.0, 1.0, 0.2, true, 1.0e-4, trng, 10),
                IntegratedWoodSaxonInterp(10, 0.6, 5.0, 1.0, 0.2, 10, 1.0e-4, 60, trng, true),
            )
            s = nuc.sampler
            for _ in 1:500
                x = ntuple(i -> (2 * rand(trng) - 1) * s.half_width[i], length(s.half_width))
                @test s.target_probability(x..., s.params...) <= s.bound
            end
        end

        # an envelope violation is reported, not silently clipped
        bad = MonteCarloGlauber.RejectionSampler((x) -> 2.0, (), (1.0,), 1.0)
        @test_throws ErrorException rand(trng, bad)
        @test_throws ArgumentError MonteCarloGlauber.RejectionSampler((x) -> 1.0, (), (1.0,), 0.0)
    end

    @testset "Deformed density: trig-free evaluation matches the spherical-harmonic formula" begin
        # reference implementation with explicit angles, as written before the optimisation
        function reference(x, y, z, α, R, ρ₀, w, beta2, beta3, beta4, gamma)
            r = hypot(x, y, z)
            theta = acos(z / r)
            phi = sign(y) * acos(x / hypot(x, y))
            ctheta = cos(theta)
            Y20 = sqrt(5) / 4 * 1 / sqrt(pi) * (3 * ctheta^2 - 1)
            Y22 = sqrt(15 / 32) * 1 / sqrt(pi) * (1 - ctheta^2) * cos(2 * phi)
            Y30 = 1 / 4 * sqrt(7 / pi) * (5 * ctheta^3 - 3 * ctheta)
            Y40 = 3 / 16 * 1 / sqrt(pi) * (35 * ctheta^4 - 30 * ctheta^2 + 3)
            Reff = R * (1 + beta2 * (cos(gamma) * Y20 + sqrt(2) * sin(gamma) * Y22) + beta3 * Y30 + beta4 * Y40)
            return ρ₀ * (1 + w * (r / Reff)^2) / (1 + exp((r - Reff) / α))
        end
        pars = (0.6, 6.81, 1.0, 0.1, 0.28, 0.093, 0.05, 0.3)   # all deformations and triaxiality non-zero
        trng = MersenneTwister(7)
        for _ in 1:200
            x, y, z = 20 .* (rand(trng, 3) .- 0.5)
            @test MonteCarloGlauber.density_WS_deformed(x, y, z, pars...) ≈ reference(x, y, z, pars...) rtol = 1.0e-12
        end
        u = Uranium()
        @test u.sampler.target_probability === MonteCarloGlauber.density_WS_deformed_coeff
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

    @testset "Pair-loop cutoff and default impact-parameter range" begin
        w = 0.5
        n = Lead()
        participants = Participants(n, n, w, 2760.0, 1.0, 0.0; Nr = 8, Nth = 8)
        # collision probability at the cutoff distance is at the 1e-7 threshold
        b2 = participants.b2_cutoff
        @test MonteCarloGlauber.binary_impact_parameter_probability(b2, participants) ≈ 1.0e-7 rtol = 1.0e-3
        @test MonteCarloGlauber.binary_impact_parameter_probability(0.0, participants) > 0.9
        # default b_max covers both nuclei plus the collision range (Pb-Pb: ≈ 26 fm)
        bmax = maximum(participants.inpact_parameter_magitude)
        @test 2 * n.R + sqrt(b2) < bmax < 3 * (2 * n.R)
        @test MonteCarloGlauber.nuclear_extent(Oxygen()) > 0
        @test MonteCarloGlauber.nuclear_extent(n) ≈ n.R + 8 * n.α
        # sampled events stay inside the range and have at least one collision
        evts = [rand(rng, participants) for _ in 1:20]
        @test all(0.0 <= e.b <= bmax for e in evts)
        @test all(e.n_coll > 0 for e in evts)
    end

    @testset "Light ion tabulated events" begin
        oxygen = Oxygen()
        neon = Neon()

        @test oxygen isa TabulatedEvent
        @test neon isa TabulatedEvent
        @test isfile(oxygen.name)
        @test isfile(neon.name)
        @test oxygen.N_nucleon > 0
        @test neon.N_nucleon > 0
        @test oxygen.N_configs > 0
        @test neon.N_configs > 0
        @test size(oxygen) == (oxygen.N_nucleon, 2)
        @test size(neon) == (neon.N_nucleon, 2)

        sample_o = rand(rng, oxygen)
        sample_n = rand(rng, neon)
        @test size(sample_o) == (oxygen.N_nucleon, 2)
        @test size(sample_n) == (neon.N_nucleon, 2)
        @test eltype(sample_o) == Float64
        @test eltype(sample_n) == Float64
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
        borders = mult_sorted[1:n_event_per_bin:(n_event_per_bin * 100)]

        @test length(borders) == 100
        @test borders[1] == 100.0
        @test borders[end] == 1.0
        @test issorted(borders, rev = true)
    end

    @testset "InverseFunction" begin
        invf = InverseFunction(x -> 2 * x)
        @test isapprox(invf(4.0), 2.0; atol = 1.0e-6)
    end

    @testset "Generate Background and Two Point Function" begin
        bins = [10, 20]
        fmGeV = 5.0
        entropy(T) = 47.5 * 4 * T^3 * pi^2 / 90 * fmGeV^3
        entropyToTemp = InverseFunction(entropy)
        dSdT(T) = 1 / (3 * 47.5 * 4 * T^2 * pi^2 / 90 * fmGeV^2)

        bg, twpt = MonteCarloGlauber.generate_bg_twpt_fct(
            entropyToTemp,
            dSdT,
            1.0,
            Lead(),
            Lead(), 0.5,
            1.0,
            0.0,
            2760, bins, [2, 3]; minBiasEvents = 1000, r_grid = 0:1.0:5
        )
        @test all(isreal.(bg))
        @test all(isreal.(twpt))

    end

    @testset "Generate Background and ncoll" begin

        bins = [10, 20]
        fmGeV = 5.0
        entropy(T) = 47.5 * 4 * T^3 * pi^2 / 90 * fmGeV^3
        entropyToTemp = InverseFunction(entropy)

        bg, ncoll = MonteCarloGlauber.generate_bg_ncoll(
            entropyToTemp,
            1.0,
            Lead(),
            Lead(),
            0.5,
            1.0,
            0.0,
            2760,
            bins;
            minBiasEvents = 1000,
            r_grid = 0:1.0:2
        )

        @test all(isreal.(bg))
        @test all(isreal.(ncoll))

        @test size(bg) == (length(bins), length(0:1.0:2))
        @test size(ncoll) == (length(bins), length(0:1.0:2))
    end
end
