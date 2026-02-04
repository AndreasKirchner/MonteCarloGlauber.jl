"""
    NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,B,samp}

A 3‑dimensional deformed Woods–Saxon nucleus sampler that produces `N_nucleon` 3‑vectors per draw.

Fields
- `N_nucleon::Int64` — number of nucleons returned per draw (each draw yields an `(N_nucleon, 3)` matrix of coordinates).
- `α::T` — surface diffuseness parameter controlling the skin thickness.
- `R::R` — nuclear radius parameter.
- `ρ₀::C` — central density normalization.
- `w::W` — quadratic radial modulation coefficient.
- `beta2, beta3, beta4` — deformation parameters (spherical-harmonic multipoles) affecting the nuclear surface.
- `gamma::D4` — orientation (triaxiality) angle used with `beta2` to mix Y20 and Y22 components.
- `burning_time::B` — number of burn-in transitions executed when the internal sampler is constructed.
- `sampler::samp` — internal `Metropolis_Hastings` sampler used to draw 3D candidate positions from the deformed projected density.
- `d_min::Float64` — optional minimum inter-nucleon separation parameter (units of length); when zero no explicit separation policy is applied.

Notes
- This sampler produces 3D Cartesian coordinates (x,y,z). Sampling rotates each configuration by a random `RotMatrix{3}` before returning.
"""
struct NucleiWoodSaxon3D{T, R, C, W, D1, D2, D3, D4, B, samp} <: Sampleable{ArrayLikeVariate{3}, Continuous}
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


dimension(::NucleiWoodSaxon3D{T, R, C, W, D1, D2, D3, D4, B, samp}) where {T, R, C, W, D1, D2, D3, D4, B, samp} = 3
nucleon(s::NucleiWoodSaxon3D{T, R, C, W, D1, D2, D3, D4, B, samp}) where {T, R, C, W, D1, D2, D3, D4, B, samp} = s.N_nucleon

Base.size(s::NucleiWoodSaxon3D{T, R, C, W, D1, D2, D3, D4, B, samp}) where {T, R, C, W, D1, D2, D3, D4, B, samp} = (s.N_nucleon, 3)
function Base.eltype(s::NucleiWoodSaxon3D{T, R, C, W, D1, D2, D3, D4, B, samp}) where {T, R, C, W, D1, D2, D3, D4, B, samp}
    return promote_type(T, R, C, W, D1, D2, D3, D4)
end


"""
    NucleiWoodSaxon3D(rng::AbstractRNG, N_nucleon, α, R, ρ₀, w, beta2, beta3, beta4, gamma, burning_time, dmin)

Construct a `NucleiWoodSaxon3D` sampler.

Arguments
- `rng` — random number generator used to seed the internal MCMC sampler.
- `N_nucleon` — number of nucleons per sample.
- `α, R, ρ₀, w` — Woods–Saxon parameters.
- `beta2, beta3, beta4, gamma` — deformation parameters and triaxiality angle.
- `burning_time` — number of burn-in transitions for the internal `Metropolis_Hastings` sampler.
- `dmin` — minimum inter-nucleon separation (units of length); may be zero to disable.

Behavior
- Builds a 3‑component Normal proposal (σ = R for each Cartesian component) and a `Metropolis_Hastings` sampler targeted on `density_WS_deformed`.
"""
function NucleiWoodSaxon3D(rng::AbstractRNG, N_nucleon, α, R, ρ₀, w, beta2, beta3, beta4, gamma, burning_time, dmin)

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
    NucleiWoodSaxon3D(N_nucleon, α, R, ρ₀, w, beta2, beta3, beta4, gamma)

Convenience constructor using sensible defaults:
- `rng = Random.default_rng()`
- `burning_time = 1000`
- `dmin = 0.0` (no minimum separation)

Returns a ready-to-use `NucleiWoodSaxon3D` instance that produces `(N_nucleon, 3)` arrays of Cartesian coordinates when sampled with `rand`.

Example
```julia
julia> s = NucleiWoodSaxon3D(100, 0.5, 5.4, 0.16, 0.0, 0.0, 0.0, 0.0)
julia> sample = rand(s)  # single configuration; shape (100, 3)
```
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

Fill `x` in-place with a single sampled nuclear configuration. The expected shape of `x` is `(s.N_nucleon, 3)`.

Behavior
- For each nucleon (row) `i` the function draws a 3‑vector `xpos = rand(rng, s.sampler)` from the internal `Metropolis_Hastings` sampler.
- A single random `RotMatrix{3}` is sampled and applied to every drawn configuration so nuclei are returned with random orientation.
- If `s.d_min` is greater than zero a minimum-inter-nucleon-separation policy may be applied (implementation-dependent); when zero no separation policy is enforced.

Notes
- This is the low-level, in-place sampling hook used by `Distributions.rand` and the package's high-level APIs.
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
