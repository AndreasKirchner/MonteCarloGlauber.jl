"""
    NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,B,samp}

Representation of a deformed 3D Woods–Saxon nucleus suitable for sampling.

Fields
- `N_nucleon::Int64` — number of nucleons per sample.
- `α` — surface diffuseness parameter.
- `R` — nuclear radius parameter.
- `ρ₀` — central density normalization.
- `w` — quadratic modulation coefficient.
- `beta2, beta3, beta4` — deformation parameters (spherical-harmonic multipoles).
- `gamma` — triaxiality angle used with `beta2`.
- `burning_time` — number of burn-in transitions for the internal MCMC sampler.
- `sampler` — internal `Metropolis_Hastings` sampler that proposes 3D locations according to the deformed density.
- `d_min` — minimum inter-nucleon separation (when zero no separation enforcement is performed).

Behavior
- Samples are drawn in 3D using the internal sampler and then a random 3×3 rotation is applied; the object returns transverse (x,y) coordinates (an `(N_nucleon, 2)` array) by projecting the rotated 3D positions.

Notes
- The type implements `Sampleable{ArrayLikeVariate{2}, Continuous}` and is intended for use with `rand`/`Distributions`.
"""
struct NucleiWoodSaxon3D{T, R, C, W, D1, D2, D3, D4, B, samp} <: Sampleable{ArrayLikeVariate{2}, Continuous}
    N_nucleon::Int64
    α::T
    R::R
    ρ₀::C
    w::W
    beta2::D1
    beta3::D2
    beta4::D3
    gamma::D4
    burning_time::B
    sampler::samp
    d_min::Float64
end

"""
    density_WS_deformed(x, y, z, α, R, ρ₀, w, beta2, beta3, beta4, gamma)

Evaluate the deformed 3D Woods–Saxon density at Cartesian coordinates `(x,y,z)`.

The local effective radius is modified by spherical harmonics up to l=4:

Reff(θ,φ) = R * (1 + β₂*(cosγ·Y₂₀ + √2 sinγ·Y₂₂) + β₃·Y₃₀ + β₄·Y₄₀).

Arguments
- `x,y,z` — Cartesian coordinates where the density is evaluated.
- `α` — surface diffuseness parameter (skin depth).
- `R` — base nuclear radius.
- `ρ₀` — central density normalization.
- `w` — quadratic radial modulation coefficient.
- `beta2, beta3, beta4` — deformation coefficients (multipoles).
- `gamma` — triaxiality angle used to mix Y20/Y22 components for β₂.

Returns
- Scalar density value ρ(x,y,z) = ρ₀ (1 + w (r/Reff)^2) / (1 + exp((r - Reff)/α)).
"""
function density_WS_deformed(x, y, z, α, R, ρ₀, w, beta2, beta3, beta4, gamma)
    #spherical coordinates
    r = hypot(x, y, z)
    theta = acos(z / r)
    phi = sign(y) * acos(x / hypot(x, y))

    #define often used functions
    ctheta = cos(theta)
    #define spherical harmonics
    Y20 = sqrt(5) / 4 * 1 / sqrt(pi) * (3 * ctheta^2 - 1)
    Y22 = sqrt(15 / 32) * 1 / sqrt(pi) * (1 - ctheta^2) * cos(2 * phi)
    Y30 = 1 / 4 * sqrt(7 / pi) * (5 * ctheta^3 - 3 * ctheta)
    Y40 = 3 / 16 * 1 / sqrt(pi) * (35 * ctheta^4 - 30 * ctheta^2 + 3)

    Reff = R * (1 + beta2 * (cos(gamma) * Y20 + sqrt(2) * sin(gamma) * Y22) + beta3 * Y30 + beta4 * Y40)

    return ρ₀ * (1 + w * (r / Reff)^2) / (1 + exp((r - Reff) / α))
end


dimension(::NucleiWoodSaxon3D{T, R, C, W, D1, D2, D3, D4, B, samp}) where {T, R, C, W, D1, D2, D3, D4, B, samp} = 2
nucleon(s::NucleiWoodSaxon3D{T, R, C, W, D1, D2, D3, D4, B, samp}) where {T, R, C, W, D1, D2, D3, D4, B, samp} = s.N_nucleon

Base.size(s::NucleiWoodSaxon3D{T, R, C, W, D1, D2, D3, D4, B, samp}) where {T, R, C, W, D1, D2, D3, D4, B, samp} = (s.N_nucleon, 2)
function Base.eltype(s::NucleiWoodSaxon3D{T, R, C, W, D1, D2, D3, D4, B, samp}) where {T, R, C, W, D1, D2, D3, D4, B, samp}
    return promote_type(T, R, C, W, D1, D2, D3, D4)
end


"""
    NucleiWoodSaxon3D(rng::AbstractRNG, N_nucleon, α, R, ρ₀, w, beta2, beta3, beta4, gamma, burning_time, dmin)

Construct a `NucleiWoodSaxon3D` with explicit RNG, burn-in and minimum-separation parameters.

Arguments
- `rng` — RNG used to seed proposals and acceptance tests in the internal MCMC sampler.
- `N_nucleon` — number of nucleons to sample per draw.
- `α, R, ρ₀, w, beta2, beta3, beta4, gamma` — Woods–Saxon and deformation parameters (see `density_WS_deformed`).
- `burning_time` — number of burn-in iterations for the internal `Metropolis_Hastings` sampler.
- `dmin` — minimum inter-nucleon separation (stored on the object; not applied here by the constructor).

Behavior
- Constructs a 3-component Normal proposal (σ=R) and a `Metropolis_Hastings` sampler targeting `density_WS_deformed`, performs the requested burn-in, and returns a fully configured `NucleiWoodSaxon3D` instance ready for sampling via `rand`.
"""
function NucleiWoodSaxon3D(rng::AbstractRNG, N_nucleon, α, R, ρ₀, w, beta2, beta3, beta4, gamma, burning_time, dmin)

    #draw=(Uniform(-2R,2R),Uniform(-2R,2R),Uniform(-2R,2R))
    draw = (Normal(0, R), Normal(0, R), Normal(0, R))

    params = (α, R, ρ₀, w, beta2, beta3, beta4, gamma)

    sampler = Metropolis_Hastings(rng, density_WS_deformed, params, draw, burning_time)

    return NucleiWoodSaxon3D(N_nucleon, α, R, ρ₀, w, beta2, beta3, beta4, gamma, burning_time, sampler, dmin)

end

function NucleiWoodSaxon3D(N_nucleon, α, R, ρ₀, w, beta2, beta3, beta4, gamma, dmin)
    rng = Random.default_rng()
    burning_time = 1000
    return NucleiWoodSaxon3D(rng, N_nucleon, α, R, ρ₀, w, beta2, beta3, beta4, gamma, burning_time, dmin)
end

"""
    NucleiWoodSaxon3D(N_nucleon,α,R,ρ₀,w,beta2,beta3,beta4,gamma)

Create a 3D distribution of nucleons in a nucleus with a Wood-Saxon density profile. The density profile is deformed by the parameters `beta2`, `beta3`, `beta4`, and `gamma`. The number of nucleons is `N_nucleon`, the radius is `R`, the diffuseness is `α`, the central density is `ρ₀`, and the skin depth is `w`.
"""
function NucleiWoodSaxon3D(N_nucleon, α, R, ρ₀, w, beta2, beta3, beta4, gamma)
    rng = Random.default_rng()
    burning_time = 1000
    dmin = 0.0
    return NucleiWoodSaxon3D(rng, N_nucleon, α, R, ρ₀, w, beta2, beta3, beta4, gamma, burning_time, dmin)
end


function Base.copy(s::NucleiWoodSaxon3D{T, R, C, W, D1, D2, D3, D4, B, samp}) where {T, R, C, W, D1, D2, D3, D4, B, samp}
    return NucleiWoodSaxon3D(s.N_nucleon, s.α, s.R, s.ρ₀, s.w, s.beta2, s.beta3, s.beta4, s.gamma, s.d_min)
end


"""
    Distributions._rand!(rng::AbstractRNG, s::NucleiWoodSaxon3D, x::DenseMatrix{<:Real})

Fill `x` in-place with a single sampled nuclear configuration from `s`.

Arguments
- `rng` — random number generator used for proposals and acceptance tests.
- `s` — the `NucleiWoodSaxon3D` sampler.
- `x` — preallocated matrix to be filled; expected shape `(s.N_nucleon, 2)` representing transverse `(x,y)` coordinates.

Behavior
- For each nucleon the function draws a 3‑vector from the internal `Metropolis_Hastings` sampler via `rand(rng, s.sampler)`.
- A single random `RotMatrix{3}` is sampled and applied to all nucleon positions; the rotated 3‑vector is projected to `(x,y)` and written into `x`.
- If `s.d_min > 0` a minimum-inter-nucleon-separation policy may be applied (implementation detail); when zero no separation policy is enforced.

Returns
- The mutated input matrix `x`.

Notes
- This is a low‑level, allocation‑efficient sampling hook used by the package's `rand` APIs.
"""
function Distributions._rand!(rng::AbstractRNG, s::NucleiWoodSaxon3D{T, R, C, W, D1, D2, D3, D4, B, samp}, x::DenseMatrix{<:Real}) where {T, R, C, W, D1, D2, D3, D4, B, samp}

    if s.d_min == 0.0
        rotmatrix = rand(rng, RotMatrix{3})
        for i in Base.axes(x, 1)
            xpos = rand(rng, s.sampler)
            rot = rotmatrix * SVector{3}(xpos)
            for j in Base.axes(x, 2)
                x[i, j] = rot[j]
            end
        end

        return x
    end

    if s.d_min > 0

        rotmatrix = rand(rng, RotMatrix{3})
        for i in Base.axes(x, 1)
            xpos = rand(rng, s.sampler)
            rot = rotmatrix * SVector{3}(xpos)
            for j in Base.axes(x, 2)
                x[i, j] = rot[j]
            end
        end

        return x

    end

end
