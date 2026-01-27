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


function Base.copy(s::TabulatedEvent{B}) where {B}
    return TabulatedEvent(s.N_nucleon, s.N_configs, s.R, s.name, copy(s.array), s.rotated)
end
