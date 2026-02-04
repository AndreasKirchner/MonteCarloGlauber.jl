"""
    IntegratedWoodSaxonInterp{T,R,C,W,B,Inte,samp}

A sampleable object representing a two-dimensional integrated Wood–Saxon nucleon density.

Fields
- `N_nucleon::Int64` — number of nucleons to draw per sample (returned sample has size `(N_nucleon, 2)`).
- `α::T` — diffuseness parameter (controls surface thickness).
- `R::R` — nuclear radius scale.
- `ρ₀::C` — central density normalization.
- `w::W` — surface diffuseness parameter.
- `burning_time::B` — number of burn-in (warmup) steps for the internal Metropolis–Hastings sampler.
- `retol::Float64` — relative tolerance passed to the numerical quadrature used when building the interpolant.
- `npoint::Int` — number of grid points per axis used to compute the interpolant.
- `interpolation::Inte` — scaled B-spline interpolant (2D) approximating the integrated density used as a target for sampling.
- `sampler::samp` — sampler object (typically `Metropolis_Hastings`) used to produce correlated samples.
- `rotated::Bool` — if `true`, each sampled configuration is rotated by a random 2×2 rotation matrix before being returned.

Notes
- Implements `Sampleable{ArrayLikeVariate{2}, Continuous}` and provides `size`, `dimension`, and sampling via `rand` / `Distributions._rand!`.
"""
struct IntegratedWoodSaxonInterp{T, R, C, W, B, Inte, samp} <: Sampleable{ArrayLikeVariate{2}, Continuous}
    N_nucleon::Int64
    α::T
    R::R
    ρ₀::C
    w::W
    burning_time::B
    retol::Float64
    npoint::Int
    interpolation::Inte
    sampler::samp
    rotated::Bool
end

dimension(::IntegratedWoodSaxonInterp{T, R, C, W, B, Inte, samp}) where {T, R, C, W, B, Inte, samp} = 2
nucleon(s::IntegratedWoodSaxonInterp{T, R, C, W, B, Inte, samp}) where {T, R, C, W, B, Inte, samp} = s.N_nucleon

Base.size(s::IntegratedWoodSaxonInterp{T, R, C, W, B, Inte, samp}) where {T, R, C, W, B, Inte, samp} = (s.N_nucleon, 2)
function Base.eltype(s::IntegratedWoodSaxonInterp{T, R, C, W, B, Inte, samp}) where {T, R, C, W, B, Inte, samp}
    return promote_type(T, R, C, W)
end


"""
    IntegratedWoodSaxonInterp(N_nucleon, α, R, ρ₀, w, burning_time, retol, Npoint; rng=Random.default_rng(), rotated=true)

Construct an `IntegratedWoodSaxonInterp` by numerically integrating the Wood–Saxon profile on a 2D grid and constructing a scaled B-spline interpolant used as the target density for sampling.

Arguments
- `N_nucleon::Integer` — number of nucleons per draw.
- `α` — diffuseness parameter (controls surface thickness).
- `R` — nuclear radius scale (used also to set grid extent and proposal width).
- `ρ₀` — central density normalization.
- `w` — surface diffuseness parameter.
- `burning_time::Integer` — number of burn-in (warmup) iterations for the Metropolis–Hastings sampler.
- `retol::Float64` — relative tolerance for the numerical quadrature (`integrate_charge`).
- `Npoint::Integer` — number of grid points per axis used to evaluate the integrated density.
- `rng::AbstractRNG` (optional) — RNG used to seed the internal sampler (default `Random.default_rng()`).
- `rotated::Bool` (optional) — if `true`, produced configurations are randomly rotated before being returned.

This function:
1. builds a square grid spanning ±3R in each axis,
2. evaluates the integrated density at each grid point using `integrate_charge`,
3. constructs a cubic B-spline interpolant and scales it to the physical grid,
4. creates a `Metropolis_Hastings` sampler using the interpolant and returns the fully populated `IntegratedWoodSaxonInterp`.

Returns
- `IntegratedWoodSaxonInterp` instance ready to produce (N_nucleon, 2) samples via `rand`.
"""
function IntegratedWoodSaxonInterp(N_nucleon, α, R, ρ₀, w, burning_time, retol, Npoint, rng = Random.default_rng(), rotated = true)

    range = 3 * R

    A_x1 = -range:(2 * range / Npoint):range
    A_x2 = -range:(2 * range / Npoint):range

    A = [integrate_charge(x1, x2, α, R, ρ₀, w, retol) for x1 in A_x1, x2 in A_x2]
    itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))

    sitp = Interpolations.scale(itp, A_x1, A_x2)

    draw = (truncated(Normal(0, R), -3R, 3R), truncated(Normal(0, R), -3R, 3R))

    params = ()
    sampler = Metropolis_Hastings(rng, sitp, params, draw, burning_time)

    return IntegratedWoodSaxonInterp(N_nucleon, α, R, ρ₀, w, burning_time, retol, Npoint, sitp, sampler, rotated)

end

"""
    IntegratedWoodSaxonInterp(N_nucleon, α, R, ρ₀, w; rng=Random.default_rng(), rotated=true)

Create an `IntegratedWoodSaxonInterp` using sensible defaults for integration and sampling parameters.

Defaults
- `burning_time = 1000`
- `retol = 1.0e-3`
- `Npoint = 100`

Use this lightweight constructor for typical usage; it delegates to the full constructor which builds the interpolant and the Metropolis–Hastings sampler. Sample configurations with `rand` (or `Distributions._rand!`) to obtain an `(N_nucleon, 2)` matrix of positions.
"""
function IntegratedWoodSaxonInterp(N_nucleon, α, R, ρ₀, w, rng = Random.default_rng(), rotated = true)
    burning_time = 1000
    retol = 1.0e-3
    Npoint = 100
    return IntegratedWoodSaxonInterp(N_nucleon, α, R, ρ₀, w, burning_time, retol, Npoint, rng, rotated)
end

function Base.copy(s::IntegratedWoodSaxonInterp{T, R, C, W, B, Inte, samp}) where {T, R, C, W, B, Inte, samp}
    return IntegratedWoodSaxonInterp(s.N_nucleon, s.α, s.R, s.ρ₀, s.w)
end


"""
    Distributions._rand!(rng::AbstractRNG, s::IntegratedWoodSaxonInterp, x::DenseMatrix{<:Real})

Fill the provided matrix `x` with a single sample from the integrated Wood–Saxon distribution stored in `s`.

Behavior
- Samples are produced by drawing positions from `s.sampler` (a `Metropolis_Hastings` sampler) using `rand(rng, s.sampler)`.
- If `s.rotated` is `true` a random 2×2 rotation (`RotMatrix{2}`) is sampled and applied to each nucleon position before writing into `x`.
- The function writes directly into `x` and returns it. The expected shape of `x` is `(s.N_nucleon, 2)`.

Notes
- This method implements the low-level sampling hook used by `rand`/`Distributions` and mutates the provided buffer for efficiency.
"""
function Distributions._rand!(rng::AbstractRNG, s::IntegratedWoodSaxonInterp{T, R, C, W, B, Inte, samp}, x::DenseMatrix{<:Real}) where {T, R, C, W, B, Inte, samp}
    if s.rotated
        rotmatrix = rand(rng, RotMatrix{2})
        for i in Base.axes(x, 1)
            xpos = rand(rng, s.sampler)
            rot = rotmatrix * SVector{2}(xpos)
            for j in Base.axes(x, 2)
                x[i, j] = rot[j]
            end
        end
        return x
    end

    if !s.rotated

        for i in Base.axes(x, 1)
            xpos = rand(rng, s.sampler)
            for j in Base.axes(x, 2)
                x[i, j] = xpos[j]
            end
        end
        return x
    end


end
