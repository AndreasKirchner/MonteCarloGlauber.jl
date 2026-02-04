"""
    Participant

Container for a single sampled collision event.

`Participant` stores participant positions for both nuclei, per-participant weights,
binary-collision count, sampling parameters, and derived quantities such as the
impact parameter and total multiplicity. The object is callable: `evt(x, y)`
returns the fluctuating thickness at `(x, y)`.
"""
struct Participant{T, S, V, M, C, D, F} <: ValueSupport
    part1::Vector{SVector{2}{T}}
    part2::Vector{SVector{2}{S}}
    shape1::Vector{V}
    shape2::Vector{M}
    n_coll::Int64
    sub_nucleon_width::C
    shape_parameter::D
    p::F
    R1::Float64
    R2::Float64
    b::Float64
    multiplicity::Float64
end


"""
    impactParameter(evt::Participant)

Return the impact parameter `b` for a sampled event.
"""
impactParameter(x::Participant{T, S, V, M, C, D, F}) where {T, S, V, M, C, D, F} = x.b

"""
    multiplicity(evt::Participant)

Return the total multiplicity (entropy-like weight) for a sampled event.
"""
multiplicity(x::Participant{T, S, V, M, C, D, F}) where {T, S, V, M, C, D, F} = x.multiplicity

"""
    n_coll(evt::Participant)

Return the number of binary nucleon-nucleon collisions in the event.
"""
n_coll(x::Participant{T, S, V, M, C, D, F}) where {T, S, V, M, C, D, F} = x.n_coll

"""
    n_part(evt::Participant)

Return the number of participating nucleons in the event.
"""
n_part(x::Participant{T, S, V, M, C, D, F}) where {T, S, V, M, C, D, F} = length(x.part1) + length(x.part2)


Base.eltype(::Participant{T, S, V, M, C, D, F}) where {T, S, V, M, C, D, F} = promote_type(T, S)

@inline @fastmath function Tp(x, y, w)
    invw2 = 1 / (2 * w^2)
    1 / (pi) * invw2 * exp(-(x^2 + y^2) * invw2)
end

@inline @fastmath function pmeanpos(a, b, p)
    return (0.5 * (a^p + b^p))^(1 / p)
end

@inline @fastmath function pmeanzer(a, b)
    return sqrt(a * b)
end


"""
    pmean(a, b, p)

Combine two non‑negative thickness values `a` and `b` using a generalized power mean:
- If `p == 0` this returns the geometric mean `sqrt(a*b)`.
- Otherwise it returns ((a^p + b^p)/2)^(1/p) (the symmetric power mean used for the participant thickness combination).
"""
@inline function pmean(a, b, p)
    if p == 0
        return sqrt(a * b)
    end

    return pmeanpos(a, b, p)
end


function cross_section_from_energy(sqrtS) #returns the cross section for sqrtS in GeV
    return (0.073491826 * log(sqrtS) - 0.19313457) * log(sqrtS) + 3.123737545
end

@inline function (f::Participant{T, S, V, M, C, D, F})(x, y) where {T, S, V, M, C, D, F}

    return fluctuating_thickness(x, y, f)

end

"""
    fluctuating_thickness(x, y, f::Participant)

Compute the fluctuating transverse thickness at coordinates `(x,y)` for a given `Participant` event.

The thickness contributions from each participating nucleon are modelled as normalized Gaussians with width `f.sub_nucleon_width` and per‑participant weights (`f.shape1` and `f.shape2`). After summing the contributions from each nucleus the two thicknesses `ta` and `tb` are combined using `pmean(ta,tb,f.p)`.

Returns
- A scalar giving the local combined thickness value at `(x,y)`.
"""
function fluctuating_thickness(x::Num1, y::Num2, f::Participant{T, S, V, M, C, D, F}) where {Num1 <: Real, Num2 <: Real, T, S, V, M, C, D, F}

    part1 = f.part1
    part2 = f.part2
    shape1 = f.shape1
    shape2 = f.shape2
    w = f.sub_nucleon_width
    p = f.p

    ta = zero(eltype(f))
    tb = zero(eltype(f))

    @inbounds @fastmath for i in eachindex(part1)
        pa_x, pa_y = part1[i]
        ga = shape1[i]
        ta += ga * Tp(x - pa_x, y - pa_y, w) #TODO check this shift
    end

    @inbounds @fastmath for i in eachindex(part2)
        pa_x, pa_y = part2[i]
        ga = shape2[i]
        tb += ga * Tp(x - pa_x, y - pa_y, w)
    end

    return pmean(ta, tb, p) #norm((ta,tb),p)
end

"""
    fluctuating_thickness(x, y, participants::Vector{Participant})

Compute the fluctuating thickness at `(x,y)` for a vector of `Participant` events, returning an array of thickness values (one per event).

This is a convenience vectorized wrapper around the single‑event `fluctuating_thickness`.
"""
function fluctuating_thickness(x::Num1, y::Num2, f::Vector{Participant{T, S, V, M, C, D, F}}) where {Num1 <: Real, Num2 <: Real, T, S, V, M, C, D, F}
    return map(f) do f_i
        fluctuating_thickness(x, y, f_i)
    end
end


"""
    Participants{A,B,C,D,E,F,G,H,L}

Container describing a two‑nucleus collision setup and providing a `Sampleable` that generates
`Participant` events.

Fields
- `nucl1`, `nucl2` — nucleus samplers for the two colliding nuclei.
- `sub_nucleon_width` — Gaussian width (σ) used for the transverse matter distribution of sub‑nucleonic degrees of freedom.
- `inpact_parameter_magitude` — distribution for the impact parameter magnitude (radial distribution).
- `inpact_parameter_angle` — distribution for the impact parameter angle (azimuthal orientation).
- `sigma_gg` — per‑gluon (partonic) cross section used in computing binary collision probabilities.
- `shape_parameter` — Gamma distribution shape parameter used to sample per‑participant fluctuation weights.
- `total_cross_section` — inelastic nucleon‑nucleon cross section at the chosen energy (σ_NN).
- `p` — exponent used when combining thickness functions (p‑mean: `p = 0` implies geometric mean).
- `accumulation_preparation` — precomputed data used to accelerate center‑of‑mass and accumulation calculations.

Notes
- Implements `Sampleable{Univariate, Participant}`, and `rand(rng, participants)` returns a sampled `Participant` event describing one collision configuration.
"""
struct Participants{A, B, C, D, E, F, G, H, L} <: Sampleable{Univariate, Participant}
    nucl1::A
    nucl2::B
    sub_nucleon_width::C
    inpact_parameter_magitude::D
    inpact_parameter_angle::G
    sigma_gg::E
    shape_parameter::F
    total_cross_section::Float64
    p::L
    accumulation_preparation::H
end

"""
    Participants(n1,n2,w,s_NN,k,p[,b])

Constructor for the Participants distribution. It takes the following arguments:
    n1: Nucleus 1
    n2: Nucleus 2
    w: Subnucleon width
    s_NN: Energy in GeV
    k: Shape parameter
    p: Norm parameter
    b: Impact parameter distribution (optional) 
"""
function Participants

    end

function Participants(n1, n2, w, s_NN, k, p, b::Tuple{T1, T2}; Nr = 64, Nth = 32) where {T1 <: Real, T2 <: Real}

    sigma_NN = cross_section_from_energy(s_NN)
    f(sigmagg, p) = totalcross_section(w, sigmagg, sigma_NN)
    u0 = one(sigma_NN) 
    prob = NonlinearProblem{false}(f, u0)
    sol = solve(prob, SimpleNewtonRaphson())
    sigg = sol.u

    R1 = n1.R
    R2 = n2.R
    prep = prepare_accumulation(R1, R2, Nr, Nth)

    return Participants(n1, n2, w, truncated(TriangularDist(0, b[2], b[2]), b[1], b[2]), Uniform(0, 2pi), sigg, k, sigma_NN, p, prep)

end


dimension(s::Participants{A, B, C, D, E, F, G, L, H}) where {A, B, C, D, E, F, G, L, H} = 1


function Participants(n1, n2, w, s_NN, k, p; Nr = 32, Nth = 32)

    R1 = n1.R
    R2 = n2.R

    b = (0, 3(R1 + R2) + 6 * w)
    return Participants(n1, n2, w, s_NN, k, p, b; Nr = Nr, Nth = Nth)
end

function Base.eltype(::Participants{A, B, C, D, E, F, G, H, L}) where {A, B, C, D, E, F, G, H, L}
    T = promote_type(C, E, F, L)
    return Participant{T, T, T, T, C, F, L}
end


function Base.copy(s::Participants{A, B, C, D, E, F, G, H, L}) where {A, B, C, D, E, F, G, H, L}


    b = unique(Distributions.params(s.inpact_parameter_magitude))


    return Participants(
        copy(s.nucl1), copy(s.nucl2), s.sub_nucleon_width,
        truncated(TriangularDist(0, b[2], b[2]), b[1], b[2]), Uniform(0, 2pi), s.sigma_gg, s.shape_parameter, s.total_cross_section, s.p, deepcopy(s.accumulation_preparation)
    )

end


@inline function totalcross_section(w, sigmaGG, sigmaNN) #double check


    C = 6 #truncation for integral

    return C^2 / 4 + expinti(-exp(-C^2 / 4) * sigmaGG / (4 * pi * w^2)) - expinti(-sigmaGG / (4 * pi * w^2)) - sigmaNN / (4 * pi * w^2)

end


@inline @fastmath function binary_impact_parameter_probability(b2, nucleos::Participants{A, B, C, D, E, F, G, H, L}) where {A, B, C, D, E, F, G, H, L}
    w = nucleos.sub_nucleon_width
    gasussd = 1 / (4 * w^2)
    Tnn = gasussd * Base.Math.exp_fast(-b2 * gasussd) / pi
    return one(Tnn) - Base.Math.exp_fast(-nucleos.sigma_gg * Tnn)

end


#3d case

"""
    Distributions.rand(rng::AbstractRNG, nucleos::Participants)

Sample a single collision `Participant` event from the configured `Participants` object.

Arguments
- `rng` — random number generator used for draws and acceptance tests.
- `nucleos` — a `Participants` instance describing the two nuclei, sub‑nucleon width, cross sections, and other sampling parameters.

Behavior
- Draws nucleus configurations `n1`, `n2` and an impact-parameter orientation `θ_b`.
- Samples an impact-parameter magnitude `b` from `nucleos.inpact_parameter_magitude` and shifts the two nuclei by ±b/2.
- For each nucleon pair computes the binary collision probability using `binary_impact_parameter_probability` and accepts collisions via Bernoulli draws.
- If at least one binary collision is discovered, builds a `Participant` containing participating nucleons, samples per-nucleon Gamma weights, recenters to the event center-of-mass (via `center_of_mass!`), and returns the populated `Participant`.

Notes
- The function uses rejection sampling to ensure the generated event contains at least one binary collision.

Example
```julia
julia> nucleos = Participants(nucleus1, nucleus2, 0.4, 200.0, 1.0, 0.0)
julia> evt = Distributions.rand(Random.default_rng(), nucleos)
julia> evt.b  # impact-parameter magnitude
0.7
```
"""
function Distributions.rand(rng::AbstractRNG, nucleos::Participants{NUCL1, NUCL2, C, D, E, F, G, H, L}) where {NUCL1, NUCL2, C, D, E, F, G, H, L}


    R1 = nucleos.nucl1.R
    R2 = nucleos.nucl2.R
    n1 = rand(rng, nucleos.nucl1)
    n2 = rand(rng, nucleos.nucl2)
    θ_b = rand(rng, nucleos.inpact_parameter_angle)
    s_th, c_th = sincos(θ_b)
    re1 = Vector{SVector{2, eltype(nucleos.nucl1)}}()
    re2 = Vector{SVector{2, eltype(nucleos.nucl2)}}()

    while true
        b = rand(rng, nucleos.inpact_parameter_magitude)
        b_vec = SVector{2}(b * c_th / 2, b * s_th / 2)
        ncoll = 0
        @inbounds for nucl1 in axes(n1, 1)
            pos1 = SVector{2}(n1[nucl1, 1], n1[nucl1, 2])
            pos1rotshift = pos1 - b_vec
            @inbounds for nucl2 in axes(n2, 1)
                pos2 = SVector{2}(n2[nucl2, 1], n2[nucl2, 2])
                pos2rotshift = pos2 + b_vec


                v = pos1rotshift - pos2rotshift
                impact_par = dot(v, v) 
                probability = binary_impact_parameter_probability(impact_par, nucleos)

                #accepted
                if rand(rng, Bernoulli(probability))
                    push!(re1, pos1rotshift)
                    push!(re2, pos2rotshift)
                    ncoll += 1
                end
            end
        end

        if ncoll > 0
            r1 = unique(re1)
            r2 = unique(re2)
            k = nucleos.shape_parameter
            distribution = Gamma(k, 1 / k)
            shape_1 = rand(rng, distribution, length(r1))
            shape_2 = rand(rng, distribution, length(r2))
            part = Participant(r1, r2, shape_1, shape_2, ncoll, nucleos.sub_nucleon_width, nucleos.shape_parameter, nucleos.p, R1, R2, b, 0.0)
            #CoM=center_of_mass(Participant)
            mult, x_cm, y_cm = center_of_mass!(part, nucleos.accumulation_preparation)
            r_cm = SVector{2}(x_cm / mult, y_cm / mult)
            @inbounds for i in eachindex(r1)
                r1[i] = r1[i] - r_cm
            end
            @inbounds for i in eachindex(r2)
                r2[i] = r2[i] - r_cm
            end

            return Participant(r1, r2, shape_1, shape_2, ncoll, nucleos.sub_nucleon_width, nucleos.shape_parameter, nucleos.p, R1, R2, b, mult)
        end
    end
    return
end

"""
    Distributions.rand(nucleos::Participants)

Convenience wrapper that samples a `Participant` using the default RNG (`Random.default_rng()`).
This is equivalent to `Distributions.rand(Random.default_rng(), nucleos)`.
"""
function Distributions.rand(nucleos::Participants{A, B, C, D, E, F, G, H, L}) where {A, B, C, D, E, F, G, H, L}
    return Distributions.rand(Random.default_rng(), nucleos)
end

function Distributions.rand!(rng::AbstractRNG, nucleos::Participants{NUCL1, NUCL2, C, D, E, F, H, L}, A::AbstractVector{T}) where {NUCL1, NUCL2, C, D, E, F, L, H, T}

    for i in eachindex(A)
        A[i] = rand(rng, nucleos)
    end

    return A

end
