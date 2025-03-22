

struct WarpedTuple{N,T} <:ValueSupport
    
end 

mutable struct Metropolis_Hastings{DIM,T,S,D,M,uType}<:Sampleable{Univariate,WarpedTuple{DIM,uType}}
    targhet_probability::T
    params::S
    draw_distrubution::D
    burning_time::Int64
    seed_probability::M
end 

Base.eltype(::Metropolis_Hastings{DIM,T,S,D,M,uType}) where {DIM,T,S,D,M,uType} = NTuple{DIM,uType}


function Metropolis_Hastings(rng::AbstractRNG,targhet_probability::T,params::S,draw_distrubution::D,burning_time::Int64) where {T,S,D}

    xtype=ntuple(i->rand(rng,draw_distrubution[i]),length(draw_distrubution))

    seed_probability=targhet_probability(xtype...,params...)    
    
    mcmc=Metropolis_Hastings{length(xtype),T,S,D,typeof(seed_probability),eltype(xtype)}(targhet_probability,params,draw_distrubution,burning_time,seed_probability)

    for i in 1:burning_time
        rand(rng,mcmc)
    end 
    
    return mcmc
end 

function Metropolis_Hastings(targhet_probability::T,params::S,draw_distrubution::D,burning_time::Int64) where {T,S,D}
    Metropolis_Hastings(Random.default_rng(),targhet_probability,params,draw_distrubution,burning_time)
end 
 

function Base.copy(s::Metropolis_Hastings{DIM,T,S,D,M,uType}) where {DIM,T,S,D,M,uType}
    
    Metropolis_Hastings(s.targhet_probability,s.params,s.draw_distrubution,s.burning_time)
end

function mcmcaccept(rng::AbstractRNG,prob,s::Metropolis_Hastings{DIM,T,S,D,M,uType}) where {DIM,T,S,D,M,uType}

    coll=false
    while !coll
    
    x_propose= map(s.draw_distrubution) do x
        rand(rng,x)
    end

    newprob=s.targhet_probability(x_propose...,s.params...)

    deltaS =newprob/prob

    if deltaS>1
        coll=true 
    return (x_propose,newprob)

    end 

    if (rand(rng) <deltaS)
        coll=true
        return (x_propose,newprob)
    end  

    #return mcmcaccept(rng,prob,s)
      prob=newprob

    end

end 

function mcmcaccept_old(rng::AbstractRNG,prob,s::Metropolis_Hastings{DIM,T,S,D,M,uType}) where {DIM,T,S,D,M,uType}

    x_propose= map(s.draw_distrubution) do x
        rand(rng,x)
    end

    newprob=s.targhet_probability(x_propose...,s.params...)

    deltaS =newprob/prob

    if deltaS>1

    return (x_propose,newprob)

    end 

    if (rand(rng) <deltaS)
        return (x_propose,newprob)
    end  

    return mcmcaccept(rng,prob,s)

end

function Distributions.rand(rng::AbstractRNG, s::Metropolis_Hastings{DIM,T,S,D,M,uType}) where {DIM,T,S,D,M,uType}
    prob=s.seed_probability
    x,newprob = mcmcaccept(rng,prob,s)
    s.seed_probability=newprob
    return x
end

function Distributions._rand!(rng::AbstractRNG, s::Metropolis_Hastings{DIM,T,S,D,M,uType}, out::Vector{NTuple{DIM,uType}}) where {DIM,T,S,D,M,uType}
   
    @inbounds for i in eachindex(out) 
        out[i]=rand(rng,s)
    end 
    return out
end



