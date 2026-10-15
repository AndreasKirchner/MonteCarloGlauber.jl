"""
    RejectionSampler{DIM, T, S, V}

Exact sampler for an unnormalised density `target_probability(x₁, …, x_DIM, params...)`
supported on the box `[-half_width[i], half_width[i]]` in each coordinate.

Candidates are drawn uniformly in the box and accepted with probability
`target / bound`. `bound` must be a guaranteed upper bound of the target on the box
(an *envelope*); the constructor of each nucleus type supplies one analytically.
Accepted points are independent draws from the normalised target, so no burn-in and
no chain state is needed. If a candidate ever exceeds `bound` the envelope is wrong and
`rand` throws rather than silently sampling `min(target, bound)`.
"""
struct RejectionSampler{DIM, T, S, V}
    target_probability::T
    params::S
    half_width::NTuple{DIM, V}
    bound::V

    function RejectionSampler(target_probability::T, params::S, half_width::NTuple{DIM, V}, bound::Real) where {DIM, T, S, V}
        isfinite(bound) && bound > 0 || throw(ArgumentError("RejectionSampler: bound must be positive and finite, got $bound"))
        all(h -> isfinite(h) && h > 0, half_width) || throw(ArgumentError("RejectionSampler: half widths must be positive and finite, got $half_width"))
        return new{DIM, T, S, V}(target_probability, params, half_width, V(bound))
    end
end

Base.eltype(::RejectionSampler{DIM, T, S, V}) where {DIM, T, S, V} = NTuple{DIM, V}

Base.copy(s::RejectionSampler) = s

function Random.rand(rng::AbstractRNG, s::RejectionSampler{DIM, T, S, V}) where {DIM, T, S, V}
    while true
        x = ntuple(i -> (2 * rand(rng, V) - 1) * s.half_width[i], Val(DIM))
        p = s.target_probability(x..., s.params...)
        if p > s.bound
            error("RejectionSampler: target density $p at $x exceeds its envelope $(s.bound); the bound supplied by the nucleus constructor is not an upper bound")
        end
        if rand(rng, V) * s.bound < p
            return x
        end
    end
    return
end

Random.rand(s::RejectionSampler) = rand(Random.default_rng(), s)
