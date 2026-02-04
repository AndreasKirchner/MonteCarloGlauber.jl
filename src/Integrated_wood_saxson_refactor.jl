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

function density_WS(x, y, z, α, R, ρ₀, w)

    r = hypot(x, y, z)

    return ρ₀ * (1 + w * (r / R)^2) / (1 + exp((r - R) / α))

end


function integrate_charge(x, y, α, R, ρ₀, w, retol, segbuf = alloc_segbuf(typeof(1.0 * Inf), promote_type(typeof(x), typeof(x), typeof(α), typeof(R), typeof(ρ₀), typeof(w)), promote_type(typeof(x), typeof(x), typeof(α), typeof(R), typeof(ρ₀), typeof(w))))

    integral, error = quadgk(z -> density_WS(x, y, z, α, R, ρ₀, w), -1.0 * Inf, 1.0 * Inf; rtol = retol, segbuf = segbuf)
    return integral
end


"""
    IntegratedWoodSaxon(N_nucleon::Int,α,R,ρ₀,w,rotated=true,retol=1e-3
    ,rng=Random.default_rng()
    ,burning_time=1000)

Create a sampleable object that represents a Wood Saxon distribution with a given number of nucleons.
The distribution is deformed by the parameters `α`, `R`, `ρ₀`, and `w`. The number of nucleons is `N_nucleon`.
"""
function IntegratedWoodSaxon(
        N_nucleon::Int, α, R, ρ₀, w, rotated = true, retol = 1.0e-3
        , rng = Random.default_rng()
        , burning_time = 1000
    )

    segbuf = alloc_segbuf(typeof(1.0 * Inf), promote_type(typeof(α), typeof(R), typeof(ρ₀), typeof(w)), promote_type(typeof(α), typeof(R), typeof(ρ₀), typeof(w)))

    #draw=(Uniform(-3R,3R),Uniform(-3R,3R))
    draw = (Normal(0, R), Normal(0, R))
    params = (α, R, ρ₀, w, retol, segbuf)

    sampler = Metropolis_Hastings(rng, integrate_charge, params, draw, burning_time)

    return IntegratedWoodSaxon{typeof(α), typeof(R), typeof(ρ₀), typeof(w), typeof(sampler)}(N_nucleon, α, R, ρ₀, w, retol, burning_time, sampler, rotated)
end

function Base.copy(s::IntegratedWoodSaxon{T, R, C, W, samp}) where {T, R, C, W, samp}
    return IntegratedWoodSaxon(s.N_nucleon, s.α, s.R, s.ρ₀, s.w)
end


#here we define the rand of this object

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
    threaded(elm)

Create a sampleable object that will mutithread the call of rand.
"""
threaded(elm) = Threaded(elm)

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
