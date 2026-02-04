struct WarpedTuple{N, T} <: ValueSupport

end

"""
    Metropolis_Hastings{DIM, T, S, D, M, uType}

A Metropolis–Hastings MCMC sampler that produces `DIM`-tuples drawn from a target probability function.

Fields
- `target_probability::T` — a callable `f(x1,...,xDIM, params...)` returning a (proportional) probability density for the proposed point.
- `params::S` — tuple of extra parameters forwarded to `target_probability`.
- `draw_distrubution::D` — a tuple of proposal distributions (one per coordinate) used to propose new candidate components.
- `burning_time::Int64` — number of initial burn-in transitions executed when the sampler is constructed.
- `seed_probability::M` — cached probability value for the current state; updated after each accepted move.

Behavior
- Implements `Sampleable{Univariate, WarpedTuple{DIM,uType}}` and supports `rand(rng, s)` to produce single draws and `_rand!` for vectorized draws.
"""
mutable struct Metropolis_Hastings{DIM, T, S, D, M, uType} <: Sampleable{Univariate, WarpedTuple{DIM, uType}}
    target_probability::T
    params::S
    draw_distrubution::D
    burning_time::Int64
    seed_probability::M
end

Base.eltype(::Metropolis_Hastings{DIM, T, S, D, M, uType}) where {DIM, T, S, D, M, uType} = NTuple{DIM, uType}


"""
    Metropolis_Hastings(rng::AbstractRNG, target_probability, params, draw_distrubution, burning_time::Int)

Create a `Metropolis_Hastings` sampler and perform `burning_time` warm-up transitions.

Arguments
- `rng` — RNG used for proposals and acceptance tests.
- `target_probability` — function `f(x1,...,xDIM, params...)` returning an unnormalized probability density for candidate points.
- `params` — additional parameters forwarded to `target_probability` on each evaluation.
- `draw_distrubution` — tuple of proposal distributions used to sample candidate components (one per dimension).
- `burning_time` — number of burn-in transitions to execute before the sampler is returned.

Behavior
- An initial seed state is drawn from the proposal distributions, its probability is cached, and `burning_time` calls to `rand` are executed to equilibrate the chain.

Returns
- A configured `Metropolis_Hastings` object ready for sampling via `rand(rng, mcmc)`.
"""
function Metropolis_Hastings(rng::AbstractRNG, target_probability::T, params::S, draw_distrubution::D, burning_time::Int64) where {T, S, D}

    xtype = ntuple(i -> rand(rng, draw_distrubution[i]), length(draw_distrubution))

    seed_probability = target_probability(xtype..., params...)

    mcmc = Metropolis_Hastings{length(xtype), T, S, D, typeof(seed_probability), eltype(xtype)}(target_probability, params, draw_distrubution, burning_time, seed_probability)

    for i in 1:burning_time
        rand(rng, mcmc)
    end

    return mcmc
end

"""
    Metropolis_Hastings(target_probability, params, draw_distrubution, burning_time::Int)

Convenience constructor that uses `Random.default_rng()` as the RNG. See the full constructor for argument semantics.
"""
function Metropolis_Hastings(target_probability::T, params::S, draw_distrubution::D, burning_time::Int64) where {T, S, D}
    return Metropolis_Hastings(Random.default_rng(), target_probability, params, draw_distrubution, burning_time)
end


function Base.copy(s::Metropolis_Hastings{DIM, T, S, D, M, uType}) where {DIM, T, S, D, M, uType}

    return Metropolis_Hastings(s.target_probability, s.params, s.draw_distrubution, s.burning_time)
end

"""
    mcmcaccept(rng, prob, s::Metropolis_Hastings)

Attempt to propose and accept a new candidate state for the MCMC chain using the Metropolis–Hastings criterion.

Arguments
- `rng` — RNG used to draw proposals and acceptance random variates.
- `prob` — current cached probability of the chain's state.
- `s` — the `Metropolis_Hastings` sampler instance providing proposal distributions and the target probability function.

Returns
- `(x_propose, newprob)` — the accepted candidate `x_propose` (an NTuple of length `DIM`) and its probability `newprob`.

Behavior
- Proposals are drawn component-wise from `s.draw_distrubution`.
- Uses the standard Metropolis acceptance ratio `newprob / prob`. If `newprob > prob` the proposal is accepted; otherwise it's accepted with probability `newprob/prob`.
- This function loops until it accepts a proposal and returns the accepted state.
"""
function mcmcaccept(rng::AbstractRNG, prob, s::Metropolis_Hastings{DIM, T, S, D, M, uType}) where {DIM, T, S, D, M, uType}

    coll = false
    while !coll

        x_propose = map(s.draw_distrubution) do x
            rand(rng, x)
        end

        newprob = s.target_probability(x_propose..., s.params...)

        deltaS = newprob / prob

        if deltaS > 1
            coll = true
            return (x_propose, newprob)

        end

        if (rand(rng) < deltaS)
            coll = true
            return (x_propose, newprob)
        end


    end

    return
end


"""
    Distributions.rand(rng::AbstractRNG, s::Metropolis_Hastings)

Produce a single MCMC draw from the Metropolis–Hastings sampler `s`.

Behavior
- Uses the cached `s.seed_probability` as the current state's probability and calls `mcmcaccept` to obtain an accepted proposal and its probability.
- Updates `s.seed_probability` with the newly accepted state's probability and returns the sample (an `NTuple{DIM, uType}`).

Notes
- This function is the hook used by `rand` and `Distributions` APIs and is intentionally allocation-light to support efficient repeated sampling.
"""
function Distributions.rand(rng::AbstractRNG, s::Metropolis_Hastings{DIM, T, S, D, M, uType}) where {DIM, T, S, D, M, uType}
    prob = s.seed_probability
    x, newprob = mcmcaccept(rng, prob, s)
    s.seed_probability = newprob
    return x
end

function Distributions._rand!(rng::AbstractRNG, s::Metropolis_Hastings{DIM, T, S, D, M, uType}, out::Vector{NTuple{DIM, uType}}) where {DIM, T, S, D, M, uType}

    @inbounds for i in eachindex(out)
        out[i] = rand(rng, s)
    end
    return out
end
