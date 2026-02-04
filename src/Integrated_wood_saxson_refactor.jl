"""
    IntegratedWoodSaxon{T,R,C,W,samp}

Representation of a two-dimensional, integrated Wood–Saxon nucleon distribution suitable for sampling.

Fields
- `N_nucleon::Int64` — number of nucleons returned per draw (samples have shape `(N_nucleon, 2)`).
- `α::T` — surface diffuseness (or deformation parameter feeding into the radial dependence).
- `R::R` — characteristic nuclear radius used to set the spatial scale.
- `ρ₀::C` — central density normalization.
- `w::W` — optional quadratic deformation coefficient that modulates the radial profile.
- `retol::Float64` — relative tolerance passed to the quadrature used when building the projected density.
- `burning_time::Int64` — number of burn-in (warmup) iterations used by the internal `Metropolis_Hastings` sampler.
- `sampler::samp` — sampler object (typically `Metropolis_Hastings`) that draws correlated 2D positions from the projected density.
- `rotated::Bool` — when `true`, each sampled configuration is post-multiplied by a random `RotMatrix{2}` so the orientation is randomized.

Notes
- Implements `Sampleable{ArrayLikeVariate{2}, Continuous}` and is intended to be used with `rand`/`Distributions` to obtain `(N_nucleon, 2)` coordinate samples.
- The sampler targets the integrated (projected) density obtained by integrating the Woods–Saxon density along the longitudinal axis.
"""
struct IntegratedWoodSaxon{T, R, C, W, samp} <: Sampleable{ArrayLikeVariate{2}, Continuous}
    N_nucleon::Int64
    α::T
    R::R
    ρ₀::C
    w::W
    retol::Float64
    burning_time::Int64
    sampler::samp
    rotated::Bool
end

dimension(::IntegratedWoodSaxon{T, R, C, W, samp}) where {T, R, C, W, samp} = 2
nucleon(s::IntegratedWoodSaxon{T, R, C, W, samp}) where {T, R, C, W, samp} = s.N_nucleon

Base.size(s::IntegratedWoodSaxon{T, R, C, W, samp}) where {T, R, C, W, samp} = (s.N_nucleon, 2)
Base.eltype(s::IntegratedWoodSaxon{T, R, C, W, samp}) where {T, R, C, W, samp} = promote_type(T, R, C, W)

"""
    density_WS(x, y, z, α, R, ρ₀, w)

Evaluate the (deformed) Woods–Saxon density at cartesian coordinates `(x,y,z)`.

Arguments
- `x, y, z` — Cartesian coordinates where the density is evaluated.
- `α` — diffuseness parameter (controls surface thickness).
- `R` — radius parameter.
- `ρ₀` — central density normalization.
- `w` — optional quadratic deformation coefficient; when nonzero applies a `(r/R)^2` modulation.

Returns
- The scalar density value ρ(x,y,z) = ρ₀ (1 + w (r/R)^2) / (1 + exp((r - R)/α)), where `r = √(x^2 + y^2 + z^2)`.
"""
function density_WS(x, y, z, α, R, ρ₀, w)

    r = hypot(x, y, z)

    return ρ₀ * (1 + w * (r / R)^2) / (1 + exp((r - R) / α))

end


"""
    integrate_charge(x, y, α, R, ρ₀, w, retol, segbuf=...)

Compute the longitudinally integrated (projected) Woods–Saxon density at transverse position `(x,y)` by integrating
`density_WS(x,y,z,...)` along `z` from `-Inf` to `Inf` using `quadgk`.

Arguments
- `x, y` — transverse coordinates where the projected (2D) density is evaluated.
- `α, R, ρ₀, w` — same parameters as for `density_WS` (diffuseness, radius, normalization, deformation).
- `retol::Float64` — relative tolerance forwarded to `quadgk`.
- `segbuf` — optional segmentation buffer for `quadgk` to improve performance on repeated calls.

Returns
- A scalar giving the integrated density ∫_{-∞}^{∞} ρ(x,y,z) dz.
"""
function integrate_charge(x, y, α, R, ρ₀, w, retol, segbuf = alloc_segbuf(typeof(1.0 * Inf), promote_type(typeof(x), typeof(x), typeof(α), typeof(R), typeof(ρ₀), typeof(w)), promote_type(typeof(x), typeof(x), typeof(α), typeof(R), typeof(ρ₀), typeof(w))))

    integral, error = quadgk(z -> density_WS(x, y, z, α, R, ρ₀, w), -1.0 * Inf, 1.0 * Inf; rtol = retol, segbuf = segbuf)
    return integral
end


"""
    IntegratedWoodSaxon(N_nucleon::Int, α, R, ρ₀, w; rotated=true, retol=1e-3, rng=Random.default_rng(), burning_time=1000)

Construct an `IntegratedWoodSaxon` sampler which draws `N_nucleon` transverse positions per sample by Monte Carlo from the projected
(2D) Woods–Saxon density. The projected density is obtained by integrating the 3D Woods–Saxon profile along the longitudinal axis
using `quadgk`, and the resulting function is used as the target for a `Metropolis_Hastings` sampler.

Arguments
- `N_nucleon::Int` — number of nucleons per draw (sample shape `(N_nucleon,2)`).
- `α, R, ρ₀, w` — density parameters passed to `density_WS` and `integrate_charge`.
- `rotated::Bool` (default `true`) — if `true` apply a random 2×2 rotation to each configuration to randomize orientation.
- `retol::Float64` (default `1e-3`) — relative tolerance passed to `quadgk` during projection.
- `rng::AbstractRNG` — RNG used to seed the internal `Metropolis_Hastings` sampler.
- `burning_time::Int` — number of burn-in iterations for the internal sampler.

Notes
- This constructor builds a `Metropolis_Hastings` sampler with `integrate_charge` as the target and returns a fully configured
  `IntegratedWoodSaxon{typeof(α), typeof(R), typeof(ρ₀), typeof(w), typeof(sampler)}` ready for sampling with `rand`.
"""
function IntegratedWoodSaxon(
        N_nucleon::Int, α, R, ρ₀, w, rotated = true, retol = 1.0e-3
        , rng = Random.default_rng()
        , burning_time = 1000
    )

    segbuf = alloc_segbuf(typeof(1.0 * Inf), promote_type(typeof(α), typeof(R), typeof(ρ₀), typeof(w)), promote_type(typeof(α), typeof(R), typeof(ρ₀), typeof(w)))

    draw = (Normal(0, R), Normal(0, R))
    params = (α, R, ρ₀, w, retol, segbuf)

    sampler = Metropolis_Hastings(rng, integrate_charge, params, draw, burning_time)

    return IntegratedWoodSaxon{typeof(α), typeof(R), typeof(ρ₀), typeof(w), typeof(sampler)}(N_nucleon, α, R, ρ₀, w, retol, burning_time, sampler, rotated)
end

"""
    copy(s::IntegratedWoodSaxon)

Return a lightweight copy of `s` with a fresh internal sampler state.

The returned value is a new `IntegratedWoodSaxon` constructed with the same density parameters and defaults; the internal
`Metropolis_Hastings` sampler is re-created so the sampler state is independent from the original.
"""
function Base.copy(s::IntegratedWoodSaxon{T, R, C, W, samp}) where {T, R, C, W, samp}
    return IntegratedWoodSaxon(s.N_nucleon, s.α, s.R, s.ρ₀, s.w)
end


"""
    Distributions._rand!(rng::AbstractRNG, s::IntegratedWoodSaxon, x::DenseMatrix{<:Real})

Fill the provided matrix `x` with a single sample from the projected Woods–Saxon density stored in `s`.

Behavior
- For each nucleon (row) `i` in `x` the function draws a transverse position `xpos = rand(rng, s.sampler)` (a 2-vector).
- If `s.rotated` is `true` a random `RotMatrix{2}` is generated and applied to every nucleon position; otherwise positions are written
  directly.
- Mutates `x` in-place and returns it. The expected shape is `(s.N_nucleon, 2)`.

Notes
- This low-level `_rand!` hook is used by the `Distributions` APIs to implement efficient allocation-free sampling.
"""
function Distributions._rand!(rng::AbstractRNG, s::IntegratedWoodSaxon{T, R, C, W, samp}, x::DenseMatrix{<:Real}) where {T, R, C, W, samp}
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


"""
    Threaded(elem[, nbuffers])

Wrapper that stores thread-local copies of a sampleable object to enable
multithreaded sampling with `rand`.
"""
struct Threaded{T, N} <: Sampleable{ArrayLikeVariate{N}, Continuous}
    copy_buff::T
    dim::Int64
end

Base.size(s::Threaded{T, N}) where {T, N} = (s.copy_buff[1].N_nucleon, s.dim)
Base.eltype(s::Threaded{T, N}) where {T, N} = eltype(first(s.copy_buff))


function Threaded(elem::T, nbuffers) where {T}
    copy_buff = map(i -> copy(elem), 1:nbuffers)
    return Threaded{typeof(copy_buff), dimension(elem)}(copy_buff, dimension(elem))
end
function Threaded(elem)
    return Threaded(elem, 2 * nthreads())
end

"""
    threaded(elem)

Return a `Threaded` wrapper around `elem` to enable multithreaded sampling with `rand`.

Arguments
- `elem` — a `Sampleable` object (for example, an `IntegratedWoodSaxon`) which will be replicated into thread-local buffers.

Behavior
- Creates `2 * nthreads()` thread-local copies by default via `Threaded(elem, 2 * nthreads())`.
- Each copy has an independent internal sampler state so concurrent calls to `rand` are thread-safe and do not share sampler state.
- When sampling with `rand(rng, threaded(elem), n)` the total workload is split among the internal buffers and executed in parallel, improving throughput on multi-core machines.

Examples
```julia
julia> s = IntegratedWoodSaxon(100, 0.5, 5.4, 0.16, 0.0)   # create a sampler for 100 nucleons
julia> ts = threaded(s)                                   # wrap for multithreaded sampling
julia> samples = rand(ts, 4)                              # produce 4 samples; returned matrix concatenates them
```

Notes
- This is a lightweight convenience function; it does not change the semantics or shape of the samples returned by `elem`.
"""
threaded(elem) = Threaded(elem)

"""
    Distributions.rand(rng::AbstractRNG, s::Threaded{T,N}, n::Int64)

Generate `n` samples from `s` using multiple thread-local copies stored in `s.copy_buff`.

Behavior
- The total work is split roughly evenly across the number of buffers (`ntask = length(s.copy_buff)`).
- Each buffer produces `div(n, ntask)` samples concurrently via a threaded map/reduce (`tmapreduce`).
- Any remainder samples are produced by the first buffer and appended.
- Returns a vertically concatenated matrix of shape `(n*s.copy_buff[1].N_nucleon, s.dim)` (or shape `(n, ...)` consistent with the underlying sample shape).
"""
function Distributions.rand(rng::AbstractRNG, s::Threaded{T, N}, n::Int64) where {T, N}
    ntask = length(s.copy_buff)
    chuncks = div(n, ntask)
    reminder = n - ntask * chuncks


    result = tmapreduce(vcat, s.copy_buff) do elm
        rand(rng, elm, chuncks)
    end

    if reminder > 0
        return vcat(result, rand(rng, s.copy_buff[1], reminder))
    end
    return result
end

Distributions.rand(s::Threaded{T, N}, n::Int64) where {T, N} = Distributions.rand(Random.default_rng(), s, n)

Distributions.rand(s::Threaded{T, N}) where {T, N} = Distributions.rand(Random.default_rng(), s, 1)
