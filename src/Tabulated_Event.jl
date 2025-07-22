using HDF5

struct TabulatedEvent{B} <:Sampleable{ArrayLikeVariate{2},Continuous}
    N_nucleon::Int64
    N_configs::Int64
    R::Float64
    name::String
    array::B
    rotated::Bool

end 

function TabulatedEvent(filename)
    fid = h5open(filename, "r")
    configs=fid["configs"]
    radius=Float64(mean(map(i->mean(map(eachcol(configs[:,:,i])) do ll
        hypot(ll...)
        end),1:size(configs)[3])))
    #TabulatedEvent(size(configs)[2],size(configs)[3],radius,filename,Float64.(read(configs)),true)
    N_nucleon=size(configs)[2]
    N_configs=size(configs)[3]
    array=Float64.(read(configs))
    close(fid)
    TabulatedEvent(N_nucleon,N_configs,radius,filename,array,true)
end

function drawIndeces(nevents,s::TabulatedEvent{B}) where {B}
    return rand(1:s.N_configs,nevents)
 end


 Base.size(s::TabulatedEvent{B}) where {B} = (s.N_nucleon,2)


function Distributions._rand!(rng::AbstractRNG, s::TabulatedEvent{B}, x::DenseMatrix{<:Real}) where {B}
    index=drawIndeces(1,s)[1]
    rotmatrix=rand(rng,RotMatrix{3})
    config=s.array[:,:,index]
    for i in Base.axes(x,1)
        #@show i
        rot = rotmatrix*config
        #@show rot
        for j in Base.axes(x,2)
            x[i,j] = rot[j,i]
        end

    end
    return x
end 


function Base.copy(s::TabulatedEvent{B}) where {B}
    TabulatedEvent(s.N_nucleon,s.N_configs,s.R,s.name,s.array,s.rotated)
end 


struct Threaded{T,N} <: Sampleable{ArrayLikeVariate{N},Continuous}
    copy_buff::T
    dim::Int64  
end

Base.size(s::Threaded{T,N}) where {T,N} = (s.copy_buff[1].N_nucleon,s.dim)
Base.eltype(s::Threaded{T,N}) where {T,N} = eltype(first(s.copy_buff))


function Threaded(elem::T,nbuffers) where {T}
    copy_buff=map(i->copy(elem),1:nbuffers)
    Threaded{typeof(copy_buff),dimension(elem)}(copy_buff,dimension(elem))
end 
function Threaded(elem)
    Threaded(elem,2*nthreads())
end 

"""
    threaded(elm)

Create a sampleable object that will mutithread the call of rand.
"""
threaded(elm)=Threaded(elm)

function Distributions.rand(rng::AbstractRNG,s::Threaded{T,N},n::Int64) where {T,N}
    ntask=length(s.copy_buff)
    chuncks=div(n,ntask)
    reminder=n-ntask*chuncks
    

    result = tmapreduce(vcat,s.copy_buff) do elm
         rand(rng,elm,chuncks)
    end

    if reminder>0
        return vcat(result,rand(rng,s.copy_buff[1],reminder))
    end
    return result
end 

Distributions.rand(s::Threaded{T,N},n::Int64) where {T,N}= Distributions.rand(Random.default_rng(),s,n)

Distributions.rand(s::Threaded{T,N}) where {T,N} = Distributions.rand(Random.default_rng(),s,1)

