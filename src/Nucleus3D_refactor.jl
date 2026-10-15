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
- `burning_time` — unused, kept for backward compatibility.
- `sampler` — internal `RejectionSampler` that draws exact 3D locations from the deformed density.
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
    return density_WS_deformed_coeff(x, y, z, α, R, ρ₀, w, deformation_coefficients(beta2, beta3, beta4, gamma)...)
end

"""
    deformation_coefficients(beta2, beta3, beta4, gamma)

Fold the deformation parameters and the spherical-harmonic normalisations into four
coefficients `(c20, c22, c30, c40)` such that

Reff(θ,φ) = R * (1 + c20 (3cos²θ − 1) + c22 (1 − cos²θ) cos2φ + c30 (5cos³θ − 3cosθ) + c40 (35cos⁴θ − 30cos²θ + 3)).

Computed once per nucleus so the density evaluation needs no trigonometric calls.
"""
function deformation_coefficients(beta2, beta3, beta4, gamma)
    c20 = beta2 * cos(gamma) * sqrt(5 / pi) / 4
    c22 = beta2 * sqrt(2) * sin(gamma) * sqrt(15 / (32 * pi))
    c30 = beta3 * sqrt(7 / pi) / 4
    c40 = beta4 * 3 / (16 * sqrt(pi))
    return (c20, c22, c30, c40)
end

"""
    density_WS_deformed_coeff(x, y, z, α, R, ρ₀, w, c20, c22, c30, c40)

Trig-free evaluation of the deformed Woods–Saxon density using the coefficients from
`deformation_coefficients`. Uses `cosθ = z/r` and `cos2φ = (x² − y²)/(x² + y²)` directly.
This is the target function of the internal rejection sampler.
"""
@inline function density_WS_deformed_coeff(x, y, z, α, R, ρ₀, w, c20, c22, c30, c40)
    r = hypot(x, y, z)

    # spherical nucleus: Reff == R, skip the angular part
    if iszero(c20) && iszero(c22) && iszero(c30) && iszero(c40)
        return ρ₀ * (1 + w * (r / R)^2) / (1 + exp((r - R) / α))
    end

    ρ2 = x^2 + y^2
    c = r > 0 ? z / r : one(r)            # cosθ
    c2 = c^2
    cos2φ = ρ2 > 0 ? (x^2 - y^2) / ρ2 : one(ρ2)

    Reff = R * (1 + c20 * (3 * c2 - 1) + c22 * (1 - c2) * cos2φ + c30 * (5 * c2 * c - 3 * c) + c40 * (35 * c2^2 - 30 * c2 + 3))

    return ρ₀ * (1 + w * (r / Reff)^2) / (1 + exp((r - Reff) / α))
end


"""
    density_WS_deformed_bound(L, α, R, ρ₀, w, beta2, beta3, beta4, gamma)

Guaranteed upper bound of `density_WS_deformed` on the cube `[-L, L]³`.

The Fermi factor is at most 1, and `Reff ≥ R_min = R (1 − |β₂|(|cosγ| max Y₂₀ + √2 |sinγ| max|Y₂₂|) − |β₃| max Y₃₀ − |β₄| max Y₄₀)`,
so `ρ ≤ ρ₀ (1 + max(w, 0) · 3L² / R_min²)`. For `w = 0` (all built-in nuclei) this is simply `ρ₀`.
"""
function density_WS_deformed_bound(L, α, R, ρ₀, w, beta2, beta3, beta4, gamma)
    maxY20 = sqrt(5 / (4 * pi))          # at θ = 0
    maxY22 = sqrt(15 / (32 * pi))        # at θ = π/2, φ = 0
    maxY30 = sqrt(7 / (4 * pi))          # at θ = 0
    maxY40 = 3 / (2 * sqrt(pi))          # at θ = 0
    Rmin = R * (1 - abs(beta2) * (abs(cos(gamma)) * maxY20 + sqrt(2) * abs(sin(gamma)) * maxY22) - abs(beta3) * maxY30 - abs(beta4) * maxY40)
    Rmin > 0 || throw(ArgumentError("deformation too large: effective radius can reach zero (R_min = $Rmin)"))
    return ρ₀ * (1 + max(w, zero(w)) * 3 * L^2 / Rmin^2)
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
- `burning_time` — ignored; kept for backward compatibility (the sampler is exact and needs no burn-in).
- `dmin` — minimum inter-nucleon separation (stored on the object; not applied here by the constructor).

Behavior
- Constructs a `RejectionSampler` targeting `density_WS_deformed_coeff` on a box that contains the nucleus, and returns a fully configured `NucleiWoodSaxon3D` instance ready for sampling via `rand`.
"""
function NucleiWoodSaxon3D(rng::AbstractRNG, N_nucleon, α, R, ρ₀, w, beta2, beta3, beta4, gamma, burning_time, dmin)

    # box large enough to contain the deformed surface plus a 10α tail (ρ/ρ₀ ≈ e^-10 at the edge)
    L = float(R * (1 + abs(beta2) + abs(beta3) + abs(beta4)) + 10 * α)

    params = (α, R, ρ₀, w, deformation_coefficients(beta2, beta3, beta4, gamma)...)

    bound = density_WS_deformed_bound(L, α, R, ρ₀, w, beta2, beta3, beta4, gamma)
    sampler = RejectionSampler(density_WS_deformed_coeff, params, (L, L, L), bound)

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
- For each nucleon the function draws an independent 3‑vector from the internal `RejectionSampler` via `rand(rng, s.sampler)`.
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
