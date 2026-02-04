"""
    TabulatedEvent(filename)

Load a set of precomputed light-ion configurations from an HDF5 file and
construct a `Sampleable` over transverse nucleon positions.

File format
- The HDF5 file should contain a dataset named `configs` with shape
  `(3, N_nucleon, N_configs)`, where the first axis indexes Cartesian coords
  `(x,y,z)`, the second axis indexes nucleons, and the third axis indexes
  different configuration realizations.

Behavior
- Each sample selects a config index uniformly at random, applies a random
  3D rotation (if `rotated==true`), projects to the transverse plane by taking
  the first two coordinates, and returns an `(N_nucleon, 2)` matrix of positions.

Example
```julia
julia> s = TabulatedEvent("configs.h5")
julia> pos = rand(s)  # pos has shape (N_nucleon, 2)
```
"""
struct TabulatedEvent{B} <: Sampleable{ArrayLikeVariate{2}, Continuous}
    N_nucleon::Int64
    N_configs::Int64
    R::Float64
    name::String
    array::B
    rotated::Bool

end

function TabulatedEvent(filename)
    fid = h5open(filename, "r")
    configs = fid["configs"]
    radius = Float64(
        mean(
            map(
                i -> mean(
                    map(eachcol(configs[:, :, i])) do ll
                        hypot(ll...)
                    end
                ), 1:size(configs)[3]
            )
        )
    )

    N_nucleon = size(configs)[2]
    N_configs = size(configs)[3]
    array = Float64.(read(configs))
    close(fid)
    return TabulatedEvent(N_nucleon, N_configs, radius, filename, array, true)
end

function drawIndeces(nevents, s::TabulatedEvent{B}) where {B}
    return rand(1:s.N_configs, nevents)
end


Base.size(s::TabulatedEvent{B}) where {B} = (s.N_nucleon, 2)


"""
    Distributions._rand!(rng::AbstractRNG, s::TabulatedEvent, x::DenseMatrix{<:Real})

Fill the provided matrix `x` in-place with a single sampled configuration from
`s`.

Behavior
- Selects a random configuration index, samples a random 3×3 rotation matrix,
  applies the rotation to the 3×N config array, and writes the first two
  components of each rotated nucleon into the rows of `x`.
- The expected shape of `x` is `(s.N_nucleon, 2)` and the function returns
  the mutated `x`.

Notes
- This low-level `_rand!` method provides an allocation-friendly sampling
  primitive and is used internally by `Distributions.rand`/`rand`.
"""
function Distributions._rand!(rng::AbstractRNG, s::TabulatedEvent{B}, x::DenseMatrix{<:Real}) where {B}
    index = drawIndeces(1, s)[1]
    rotmatrix = rand(rng, RotMatrix{3})
    config = s.array[:, :, index]
    for i in Base.axes(x, 1)
        rot = rotmatrix * config
        for j in Base.axes(x, 2)
            x[i, j] = rot[j, i]
        end

    end
    return x
end


"""
    Distributions.rand(s::TabulatedEvent)

Convenience wrapper that samples a `TabulatedEvent` configuration using the
default RNG (`Random.default_rng()`). Equivalent to
`Distributions.rand(Random.default_rng(), s)`.
"""
function Distributions.rand(s::TabulatedEvent{B}) where {B}
    return Distributions.rand(Random.default_rng(), s)
end

function Base.copy(s::TabulatedEvent{B}) where {B}
    return TabulatedEvent(s.N_nucleon, s.N_configs, s.R, s.name, copy(s.array), s.rotated)
end
