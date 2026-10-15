"""
    RejectionSampler{DIM, T, S, V}

Exact sampler for an unnormalised density `target_probability(x₁, …, x_DIM, params...)`
supported on the box `[-half_width[i], half_width[i]]` in each coordinate.

Candidates are drawn uniformly in the box and accepted with probability
`target / bound`, where `bound ≥ max target` is estimated on a grid at
construction time (with a safety margin). Accepted points are independent
draws from the normalised target, so no burn-in and no chain state is needed.
"""
struct RejectionSampler{DIM, T, S, V}
    target_probability::T
    params::S
    half_width::NTuple{DIM, V}
    bound::V
end

"""
    RejectionSampler(target_probability, params, half_width::NTuple; nscan = 24, safety = 1.1)

Build a `RejectionSampler`. The upper bound of the target is found by scanning an
`nscan^DIM` grid over the box (an even `nscan` avoids evaluating exactly at the origin,
where `density_WS_deformed` is undefined) and multiplying by `safety`.
"""
function RejectionSampler(target_probability::T, params::S, half_width::NTuple{DIM, V}; nscan = 24, safety = 1.1) where {DIM, T, S, V}
    grid = ntuple(i -> range(-half_width[i], half_width[i], length = nscan), DIM)
    bound = -Inf
    for x in Iterators.product(grid...)
        p = target_probability(x..., params...)
        if !isnan(p) && p > bound
            bound = p
        end
    end
    isfinite(bound) && bound > 0 || error("RejectionSampler: could not find a positive bound for the target density")
    return RejectionSampler{DIM, T, S, V}(target_probability, params, half_width, V(safety * bound))
end

Base.eltype(::RejectionSampler{DIM, T, S, V}) where {DIM, T, S, V} = NTuple{DIM, V}

Base.copy(s::RejectionSampler) = s

function Random.rand(rng::AbstractRNG, s::RejectionSampler{DIM, T, S, V}) where {DIM, T, S, V}
    while true
        x = ntuple(i -> (2 * rand(rng, V) - 1) * s.half_width[i], Val(DIM))
        p = s.target_probability(x..., s.params...)
        if rand(rng, V) * s.bound < p
            return x
        end
    end
    return
end

Random.rand(s::RejectionSampler) = rand(Random.default_rng(), s)
