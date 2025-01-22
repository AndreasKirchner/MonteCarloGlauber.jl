using Distributions
using Random
using QuadGK
using Interpolations
using HDF5


function density_WS(x,y,z,α,R,ρ₀,w)

    r=hypot(x,y,z) 

    return ρ₀*(1+w*(r/R)^2)/(1+exp((r-R)/α))

end 

function density_WS_deformed(x,y,z,α,R,ρ₀,w,beta2,beta3,beta4,gamma)
    #spherical coordinates
    r=hypot(x,y,z)
    theta=acos(z/r)
    phi=sign(y)*acos(x/hypot(x,y))

    #define often used functions
    ctheta=cos(theta)
    #define spherical harmonics
    Y20=sqrt(5)/4*1/sqrt(pi)*(3*ctheta^2-1)
    Y22=sqrt(15/32)*1/sqrt(pi)*(1-ctheta^2)*cos(2*phi)
    Y30=1/4*sqrt(7/pi)*(5*ctheta^3 - 3*ctheta)
    Y40=3/16*1/sqrt(pi)*(35*ctheta^4-30*ctheta^2+3)

    Reff=R * (1 + beta2*(cos(gamma)*Y20+sqrt(2)*sin(gamma)*Y22) + beta3*Y30 + beta4*Y40)

    return ρ₀*(1+w*(r/Reff)^2)/(1+exp((r-Reff)/α))
end




function integrate_charge(x,y,α,R,ρ₀,w,retol,segbuf=alloc_segbuf(typeof(1.0*Inf), promote_type(typeof(x),typeof(x),typeof(α),typeof(R),typeof(ρ₀),typeof(w)), promote_type(typeof(x),typeof(x),typeof(α),typeof(R),typeof(ρ₀),typeof(w))))

    integral, error =quadgk(z -> density_WS(x,y,z,α,R,ρ₀,w), -1.0*Inf, 1.0*Inf;rtol=retol,segbuf=segbuf)
    return integral 
end

function integrate_charge(x,y,α,R,ρ₀,w,beta2,beta3,beta4,gamma,retol,segbuf=alloc_segbuf(typeof(1.0*Inf), promote_type(typeof(x),typeof(x),typeof(α),typeof(R),typeof(ρ₀),typeof(w),typeof(beta2),typeof(beta3),typeof(beta4),typeof(gamma)), promote_type(typeof(x),typeof(x),typeof(α),typeof(R),typeof(ρ₀),typeof(w),typeof(beta2),typeof(beta3),typeof(beta4),typeof(gamma))))

    integral, error =quadgk(z -> density_WS_deformed(x,y,z,α,R,ρ₀,w,beta2,beta3,beta4,gamma), -1.0*Inf, 1.0*Inf;rtol=retol,segbuf=segbuf)
    return integral 
end



abstract type Nuclus <:Sampleable{Multivariate,Continuous}

end 


struct IntegratedWoodSaxon{T,R,C,W,F,B,Seg} <: Nuclus
    N_nucleon::Int64
    α::T
    R::R
    ρ₀::C
    w::W
    range::F
    burning_time::B
    retol::Float64
    segbuf::Seg 
end





function IntegratedWoodSaxon(N_nucleon,α,R,ρ₀,w,range,burning_time,retol)
    segbuf=alloc_segbuf(typeof(1.0*Inf), promote_type(typeof(α),typeof(R),typeof(ρ₀),typeof(w)),promote_type(typeof(α),typeof(R),typeof(ρ₀),typeof(w)))
    IntegratedWoodSaxon(N_nucleon,α,R,ρ₀,w,range,burning_time,retol,segbuf)
end

function IntegratedWoodSaxon(N_nucleon,α,R,ρ₀,w,beta2,beta3,beta4,gamma,range,burning_time,retol)
    segbuf=alloc_segbuf(typeof(1.0*Inf), promote_type(typeof(α),typeof(R),typeof(ρ₀),typeof(w),typeof(beta2),typeof(beta3),typeof(beta4),typeof(gamma)),promote_type(typeof(α),typeof(R),typeof(ρ₀),typeof(w),typeof(beta2),typeof(beta3),typeof(beta4),typeof(gamma)))
    IntegratedWoodSaxonDeformed(N_nucleon,α,R,ρ₀,w,beta2,beta3,beta4,gamma,range,burning_time,retol,segbuf)
end

function IntegratedWoodSaxon(N_nucleon,α,R,ρ₀,w)
    range=3*R
    burning_time=1000
    retol=1e-3
    segbuf=alloc_segbuf(typeof(1.0*Inf), promote_type(typeof(α),typeof(R),typeof(ρ₀),typeof(w)),promote_type(typeof(α),typeof(R),typeof(ρ₀),typeof(w)))
    IntegratedWoodSaxon(N_nucleon,α,R,ρ₀,w,range,burning_time,retol,segbuf)
end

function IntegratedWoodSaxon(N_nucleon,α,R,ρ₀,w,beta2,beta3,beta4,gamma)
    range=3*R
    burning_time=1000
    retol=1e-3
    segbuf=alloc_segbuf(typeof(1.0*Inf), promote_type(typeof(α),typeof(R),typeof(ρ₀),typeof(w),typeof(beta2),typeof(beta3),typeof(beta4),typeof(gamma)),promote_type(typeof(α),typeof(R),typeof(ρ₀),typeof(w),typeof(beta2),typeof(beta3),typeof(beta4),typeof(gamma)))
    IntegratedWoodSaxonDeformed(N_nucleon,α,R,ρ₀,w,beta2,beta3,beta4,gamma,range,burning_time,retol,segbuf)
end

#more changes below needed

(f::IntegratedWoodSaxon{T,R,C,W,F,B,Seg})(x::AbstractVector{T1}) where {T,R,C,W,F,B,Seg,T1} = integrate_charge(x[1],x[2],f.α,f.R,f.ρ₀,f.w,f.retol,f.segbuf)

(f::IntegratedWoodSaxon{T,R,C,W,F,B,Seg})(x::T2,y::T3) where {T,R,C,W,F,B,Seg,T2,T3} = integrate_charge(x,y,f.α,f.R,f.ρ₀,f.w,f.retol,f.segbuf)

Base.length(s::IntegratedWoodSaxon{T,R,C,W,F,B,Seg}) where {T,R,C,W,F,B,Seg} = 2 

is_indipendnet(::IntegratedWoodSaxon{T,R,C,W,F,B,Seg}) where {T,R,C,W,F,B,Seg} =true 

uniform_point(rng::AbstractRNG,s::IntegratedWoodSaxon{T1,R,C,W,F,B,Seg}) where {T1,R,C,W,F,B,Seg}= (rand(rng)*2*s.range-s.range ,rand(rng)*2*s.range-s.range )

function sample_nenvet_with_rotation(rng,s::IntegratedWoodSaxon{T1,R,C,W,F,B,Seg},nevent) where {T1,R,C,W,F,B,Seg}

    totsize1=nevent*s.N_nucleon
        n1_=reshape(transpose(rand(rng,s,totsize1)),(nevent,s.N_nucleon,2))
    
    
    for i in axes(n1_,1)
        r = rand(rng,RotMatrix2)
        for n in axes(n1_,2)
            pos=SVector{2}(n1_[i,n,1],n1_[i,n,2])
            rot=r*pos
            n1_[i,n,:].=rot
        end 
    
    end
    return  n1_
end 
    
function sample_nevent_with_rotation_threaded(rng,s::IntegratedWoodSaxon{T1,R,C,W,F,B,Seg},nevent) where {T1,R,C,W,F,B,Seg}
    
        totsize=nevent*s.N_nucleon
    
        n = 2* Threads.nthreads()
    
        chunks=div(totsize,n)    
        
        nthreaed=stack(tmap(1:n) do i
            rand(rng,s,chunks)
        end) 
    
        n1_=permutedims(reshape(nthreaed,(2,nevent,s.N_nucleon)),(2,3,1))
    
        
        for i in axes(n1_,1)
            r = rand(rng,RotMatrix2)
            for n in axes(n1_,2)
                pos=SVector{2}(n1_[i,n,1],n1_[i,n,2])
                rot=r*pos
                n1_[i,n,:].=rot
            end 
        
        end
        return  n1_
end



struct IntegratedWoodSaxonInterp{T,R,C,W,F,B,Inte} <: Nuclus
    N_nucleon::Int64
    α::T
    R::R
    ρ₀::C
    w::W
    range::F
    burning_time::B
    retol::Float64
    npoint::Int
    interpolation::Inte
end

function IntegratedWoodSaxonInterp(N_nucleon,α,R,ρ₀,w,range,burning_time,retol,Npoint)
    
    
    A_x1 = -range:2*range/Npoint:range
    A_x2 = -range:2*range/Npoint:range
    
    A = [integrate_charge(x1,x2,α,R,ρ₀,w,retol) for x1 in A_x1, x2 in A_x2]
    itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
    
    sitp = Interpolations.scale(itp, A_x1, A_x2)

    IntegratedWoodSaxonInterp(N_nucleon,α,R,ρ₀,w,range,burning_time,retol,Npoint,sitp)
end

function IntegratedWoodSaxonInterp(N_nucleon,α,R,ρ₀,w)
    range=3*R
    burning_time=1000
    retol=1e-3
    Npoint=100
    IntegratedWoodSaxonInterp(N_nucleon,α,R,ρ₀,w,range,burning_time,retol,Npoint)
end 

(f::IntegratedWoodSaxonInterp{T,R,C,W,F,B,Inte})(x::AbstractVector{T1}) where {T,R,C,W,F,B,Inte,T1} = f.interpolation(x[1],x[2])

(f::IntegratedWoodSaxonInterp{T,R,C,W,F,B,Inte})(x::T2,y::T3) where {T,R,C,W,F,B,Inte,T2,T3} = f.interpolation(x,y)

Base.length(s::IntegratedWoodSaxonInterp{T,R,C,W,F,B,Inte}) where {T,R,C,W,F,B,Inte} = 2 

is_indipendnet(s::IntegratedWoodSaxonInterp{T,R,C,W,F,B,Inte}) where {T,R,C,W,F,B,Inte} = true 

uniform_point(rng::AbstractRNG,s::IntegratedWoodSaxonInterp{T1,R,C,W,F,B,Inte}) where {T1,R,C,W,F,B,Inte}= (rand(rng)*2*s.range-s.range ,rand(rng)*2*s.range-s.range )


struct NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B} <:Sampleable{ArrayLikeVariate{2},Continuous}
    N_nucleon::Int64
    α::T
    R::R
    ρ₀::C
    w::W
    beta2::D1
    beta3::D2
    beta4::D3
    gamma::D4
    range::F
    burning_time::B
    d_min::Float64
end

function NucleiWoodSaxon3D(N_nucleon,α,R,ρ₀,w,beta2,beta3,beta4,gamma,dmin;range=3*R,burning_time=1000) #check with a the range
    NucleiWoodSaxon3D(N_nucleon,α,R,ρ₀,w,beta2,beta3,beta4,gamma,range,burning_time,dmin)
end 


@inline (f::NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B})(x::T1,y::T2,z::T3) where {T,R,C,W,D1,D2,D3,D4,F,B,T1<:Number,T2<:Number,T3<:Number} = density_WS_deformed(x,y,z,f.α,f.R,f.ρ₀,f.w,f.beta2,f.beta3,f.beta4,f.gamma)


@inline (f::NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B})(x::AbstractVector{T1}) where {T,R,C,W,D1,D2,D3,D4,F,B,T1} = f(x[1],x[2],x[3]) 

Base.length(s::NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B}) where {T,R,C,W,D1,D2,D3,D4,F,B} = 3*s.N_nucleon
Base.size(s::NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B}) where {T,R,C,W,D1,D2,D3,D4,F,B} = (s.N_nucleon,3)


is_indipendnet(::NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B}) where {T,R,C,W,D1,D2,D3,D4,F,B} =false

uniform_point(rng::AbstractRNG,s::NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B}) where {T,R,C,W,D1,D2,D3,D4,F,B}= (rand(rng)*2*s.range-s.range ,rand(rng)*2*s.range-s.range, rand(rng)*2*s.range-s.range)



function sample_nenvet_with_rotation(rng,s::IntegratedWoodSaxonInterp{T1,R,C,W,F,B,Inte},nevent) where {T1,R,C,W,F,B,Inte}

    totsize1=nevent*s.N_nucleon
        n1_=reshape(transpose(rand(rng,s,totsize1)),(nevent,s.N_nucleon,2))
    
    
    for i in axes(n1_,1)
        r = rand(rng,RotMatrix2)
        for n in axes(n1_,2)
            pos=SVector{2}(n1_[i,n,1],n1_[i,n,2])
            rot=r*pos
            n1_[i,n,:].=rot
        end 
    
    end
    return  n1_
end 
    
function sample_nevent_with_rotation_threaded(rng,s::IntegratedWoodSaxonInterp{T1,R,C,W,F,B,Inte},nevent) where {T1,R,C,W,F,B,Inte}
    
        totsize=nevent*s.N_nucleon
    
        n = 2* Threads.nthreads()
    
        chunks=div(totsize,n)    
        
        nthreaed=stack(tmap(1:n) do i
            rand(rng,s,chunks)
        end) 
    
        n1_=permutedims(reshape(nthreaed,(2,nevent,s.N_nucleon)),(2,3,1))
    
        
        for i in axes(n1_,1)
            r = rand(rng,RotMatrix2)
            for n in axes(n1_,2)
                pos=SVector{2}(n1_[i,n,1],n1_[i,n,2])
                rot=r*pos
                n1_[i,n,:].=rot
            end 
        
        end
        return  n1_
end

function sample_nenvet_with_rotation(rng,s::NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B},nevent) where {T,R,C,W,D1,D2,D3,D4,F,B}

    totsize1=nevent*s.N_nucleon
    n1_=reshape(transpose(rand(rng,s,totsize1)),(nevent,s.N_nucleon,2))
    
    
    for i in axes(n1_,1)
        r = rand(rng,RotMatrix2)
        for n in axes(n1_,2)
            pos=SVector{2}(n1_[i,n,1],n1_[i,n,2])
            rot=r*pos
            n1_[i,n,:].=rot
        end 
    
    end
    return  n1_
end 
    
function sample_nevent_with_rotation_threaded(rng,s::NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B},nevent) where {T,R,C,W,D1,D2,D3,D4,F,B}
    
        totsize=nevent*s.N_nucleon

    
        n = 2* Threads.nthreads()
    
        chunk=div(nevent,n)   
        leftovers=nevent-chunk*n 

        nthreaed=tmap(1:n) do i
            rand(rng,s,chunk)
        end
    
        #n1_=permutedims(reshape(nthreaed,(2,nevent,s.N_nucleon)),(2,3,1))
        bulk=permutedims(stack(collect(Iterators.flatten(nthreaed))),(3,1,2))
        leftoverrand=rand(rng,s,leftovers)
        n1_=cat(bulk,permutedims(stack(leftoverrand),(3,1,2)),dims=1)
        
        for i in axes(n1_,1)
            r = rand(rng,RotMatrix3)
            for n in axes(n1_,2)
                pos=SVector{3}(n1_[i,n,1],n1_[i,n,2],n1_[i,n,3])
                rot=r*pos
                n1_[i,n,:].=rot
            end 
        
        end
        return  @view n1_[:,:,1:2]
end

#3D stuff 

struct WoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B} <:Nuclus
    N_nucleon::Int64
    α::T
    R::R
    ρ₀::C
    w::W
    beta2::D1
    beta3::D2
    beta4::D3
    gamma::D4
    range::F
    burning_time::B
end

struct TabulatedEvent{B} <:Nuclus
    N_nucleon::Int64
    N_configs::Int64
    R::Float64
    name::String
    array::B
end 


function WoodSaxon3D(N_nucleon,α,R,ρ₀,w,beta2,beta3,beta4,gamma;range=3*R,burning_time=1000) #check with a the range
    WoodSaxon3D(N_nucleon,α,R,ρ₀,w,beta2,beta3,beta4,gamma,range,burning_time)
end 

function TabulatedEvent(filename)
    fid = h5open(filename, "r")
    configs=fid["configs"]
    radius=Float64(mean(map(i->mean(map(eachcol(configs[:,:,i])) do ll
        hypot(ll...)
        end),1:size(configs)[3])))
    TabulatedEvent(size(configs)[2],size(configs)[3],radius,filename,Float64.(read(configs)))
end

@inline (f::WoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B})(x::T1,y::T2,z::T3) where {T,R,C,W,D1,D2,D3,D4,F,B,T1<:Number,T2<:Number,T3<:Number} = density_WS_deformed(x,y,z,f.α,f.R,f.ρ₀,f.w,f.beta2,f.beta3,f.beta4,f.gamma)


@inline (f::WoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B})(x::AbstractVector{T1}) where {T,R,C,W,D1,D2,D3,D4,F,B,T1} = f(x[1],x[2],x[3]) 

Base.length(s::WoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B}) where {T,R,C,W,D1,D2,D3,D4,F,B} = 3

is_indipendnet(::WoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B}) where {T,R,C,W,D1,D2,D3,D4,F,B} =true

uniform_point(rng::AbstractRNG,s::WoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B}) where {T,R,C,W,D1,D2,D3,D4,F,B}= (rand(rng)*2*s.range-s.range ,rand(rng)*2*s.range-s.range, rand(rng)*2*s.range-s.range)

function drawIndeces(nevents,s::TabulatedEvent{B}) where {B}
    return rand(1:s.N_configs,nevents)
 end

function sample_nevent_with_rotation(rng,s::TabulatedEvent{B},nevents) where {B}
    index_list=drawIndeces(nevents,s)
    n1_=permutedims(stack(map(i->s.array[:,:,i],index_list)),(3,2,1))
    for i in axes(n1_,1)
        r = rand(rng,RotMatrix3)
        for n in axes(n1_,2)
            pos=SVector{3}(n1_[i,n,1],n1_[i,n,2],n1_[i,n,3])
            rot=r*pos
            n1_[i,n,:].=rot
        end 
    end
    return @views n1_[:,:,1:2]
end

function sample_nevent_with_rotation_threaded(rng,s::TabulatedEvent{B},nevents) where {B}
    index_list=drawIndeces(nevents,s)
    n1_=permutedims(stack(map(i->s.array[:,:,i],index_list)),(3,2,1))
    for i in axes(n1_,1)
        r = rand(rng,RotMatrix3)
        for n in axes(n1_,2)
            pos=SVector{3}(n1_[i,n,1],n1_[i,n,2],n1_[i,n,3])
            rot=r*pos
            n1_[i,n,:].=rot
        end 
    end
    return @views n1_[:,:,1:2]
end

function sample_nenvet_with_rotation(rng,s::WoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B},nevent) where {T,R,C,W,D1,D2,D3,D4,F,B}

totsize1=nevent*s.N_nucleon
    n1_=reshape(transpose(rand(rng,s,totsize1)),(nevent,s.N_nucleon,3))
    @show size(n1_)


for i in axes(n1_,1)
    r = rand(rng,RotMatrix3)
    for n in axes(n1_,2)
        pos=SVector{3}(n1_[i,n,1],n1_[i,n,2],n1_[i,n,3])
        rot=r*pos
        n1_[i,n,:].=rot
    end 

end

return @views n1_[:,:,1:2]
end 



function sample_nevent_with_rotation_threaded(rng,s::WoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B},nevent) where {T,R,C,W,D1,D2,D3,D4,F,B}

    totsize=nevent*s.N_nucleon

    n = 2* Threads.nthreads()
    chunk=div(totsize,n)   
    leftover=totsize-chunk*n
    leftoverrand=rand(rng,s,leftover)

    nthreaed=stack(tmap(1:n) do i
        rand(rng,s,chunk)
    end)
    #append!(nthreaed,rand(rng,s,leftover))
    bulk=reshape(nthreaed,(3,size(nthreaed,2)*size(nthreaed,3)))
    n1_=permutedims(reshape(cat(bulk,leftoverrand,dims=2),(3,nevent,s.N_nucleon)),(2,3,1))

    
    for i in axes(n1_,1)
        r = rand(rng,RotMatrix3)
        for n in axes(n1_,2)
            pos=SVector{3}(n1_[i,n,1],n1_[i,n,2],n1_[i,n,3])
            rot=r*pos
            n1_[i,n,:].=rot
        end 
    
    end
    return @views n1_[:,:,1:2]
end


function mcmcaccept(rng::AbstractRNG,prob,s::T) where {T}

    x_propose= uniform_point(rng,s)

    newprob=s(x_propose...)

    deltaS =newprob/prob

    if deltaS>1

    return (x_propose,newprob)

    elseif deltaS<1&&(rand(rng) <deltaS)

        return (x_propose,newprob)
    else  

    return mcmcaccept(rng,prob,s)

    end 

end 



function Distributions._rand!(rng::AbstractRNG, s::T1, x::AbstractVector{T}) where {T1<:Nuclus,T}

    x_intial= uniform_point(rng,s)
    
    intial_probability=s(x_intial...)

    bunrning_index=1 

    x_point=x_intial
    probability= intial_probability

    for i in 1:s.burning_time 
        x_point,probability=mcmcaccept(rng,probability,s)
        bunrning_index+=1
    end 

    x.=x_point


end


function Distributions._rand!(rng::AbstractRNG, s::T1, A::AbstractMatrix{T}) where {T1<:Nuclus,T}
    
    x_intial= uniform_point(rng,s)
    
    intial_probability=s(x_intial...)

    bunrning_index=1 

    x_point=x_intial
    probability= intial_probability

    for i in 1:s.burning_time 
        x_point,probability=mcmcaccept(rng,probability,s)
        bunrning_index+=1
    end 

    for a in eachcol(A)
        x_point,probability=mcmcaccept(rng,probability,s)
        a.=x_point
    end 

    return A
end

function set_Projectile(projectileString)
    if projectileString=="Pb"
        return IntegratedWoodSaxonInterp(208,0.54,6.65,1,0)#chnage all to woodsaxon3d
    elseif projectileString=="Au"
        return IntegratedWoodSaxonInterp(197,0.535,6.38,1,0)
    elseif projectileString=="Xe"
        return IntegratedWoodSaxonInterp(129,0.59,5.36,1,0)
    elseif projectileString=="Cu"
        return WoodSaxon3D(63,0.596,4.2,1.,0.,0.,0.,0.,0.)#IntegratedWoodSaxonInterp(63,0.596,4.2,1,0)
    elseif projectileString=="U"
        return WoodSaxon3D(238,0.6,6.81,1,0,0.280, 0.093, 0.0,0.0)
    else
        return "Unknown projectile"
    end
end

function Uranium()
    WoodSaxon3D(238,0.6,6.81,1,0,0.280, 0.093, 0.0,0.0)
end 

function Copper()
    WoodSaxon3D(63,0.596,4.2,1.,0.,0.,0.,0.,0.)
end 

function Lead()
    return WoodSaxon3D(208,0.54,6.65,1.,0.,0.,0.,0.,0.)
end


function Gold()
    return WoodSaxon3D(197,0.535,6.38,1,0.,0.,0.,0.,0.)
end
function Xenon()
    return WoodSaxon3D(129,0.59,5.36,1,0,0.,0.,0.,0.)
end

function sample_nenvet_with_rotation(rng,s::IntegratedWoodSaxon{T1,R,C,W,F,B,Seg},nevent) where {T1,R,C,W,F,B,Seg}

    n2_=rand(rng,s,nevent)
    @show size(n2_)
        for i in axes(n2_,1)
            r = rand(rng,RotMatrix3)
            for n in axes(n2_,2)
                pos=SVector{3}(n2_[i,n,1],n2_[i,n,2],n2_[i,n,3])
                rot=r*pos
                n2_[i,n,:].=rot
            end 
        end
    return @views n2_[:,:,1:2]
   
end 

function sample_nenvet_with_rotation(rng,s::NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B},nevent) where {T,R,C,W,D1,D2,D3,D4,F,B}
    totsize1=nevent*s.N_nucleon
    #n1_=reshape(collect(Iterators.flatten(rand(rng,s,nevent))),(nevent,s.N_nucleon,3))
    n1_=permutedims(stack(rand(rng,s,nevent)),(3,1,2))
    

    for i in axes(n1_,1)
        r = rand(rng,RotMatrix3)
        for n in axes(n1_,2)
            pos=SVector{3}(n1_[i,n,1],n1_[i,n,2],n1_[i,n,3])
            rot=r*pos
            n1_[i,n,:].=rot
        end 

    end
return @views n1_[:,:,1:2]
   # n2_=rand(rng,s,nevent)
   # @show n2_
   #     for i in axes(n2_,1)
   #         r = rand(rng,RotMatrix3)
   #         for n in axes(n2_,2)
   #             pos=SVector{3}(n2_[i,n,1],n2_[i,n,2],n2_[i,n,3])
   #             rot=r*pos
   #             n2_[i,n,:].=rot
   #         end 
   #     end
   # return @views n2_[:,:,1:2]
   
end 
    
#function sample_nevent_with_rotation_threaded(rng,s::IntegratedWoodSaxon{T1,R,C,W,F,B,Seg},nevent) where {T1,R,C,W,F,B,Seg}
    

#end


#function sample_nevent_with_rotation_threaded(rng,s::NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B},nevent) where {T,R,C,W,D1,D2,D3,D4,F,B}
    

#end

function Distributions._rand!(rng::AbstractRNG, s::NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B}, A::AbstractMatrix{T1}) where {T,R,C,W,D1,D2,D3,D4,F,B,T1}
    
    x_intial= uniform_point(rng,s)
    
    intial_probability=s(x_intial...)

    bunrning_index=1 

    x_point=x_intial
    probability= intial_probability

    for i in 1:s.burning_time 
        x_point,probability=mcmcaccept(rng,probability,s)
        bunrning_index+=1
    end 

    #The matrix A has already the good dimension (N ,3 )
    x_point,probability=mcmcaccept(rng,probability,s)
    idx=1 
    A[1,:].=x_point
    
    N_n=size(A)[1]
    countingIndex=0
    while idx<N_n || countingIndex>1000 #limit for number of tries to find a new point
        #@show countingIndex
        x_point,probability=mcmcaccept(rng,probability,s)
        countingIndex+=1
        min_d=s.d_min
        for i in 1:idx 
            distance= sqrt(sum(x->(x[1]-x[2])^2,zip(view(A,i,:),x_point))) 
            min_d= min(distance,min_d)
        end
        if(min_d>= s.d_min)
            A[idx+1,:].=x_point
            idx=idx+1
        end
    end

    return A
end



function Distributions._rand!(rng::AbstractRNG, s::NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,F,B}, A::AbstractArray{<:AbstractArray{<:Real,2}}) where {T,R,C,W,D1,D2,D3,D4,F,B}
    
    x_intial= uniform_point(rng,s)
    
    intial_probability=s(x_intial...)

    bunrning_index=1 

    x_point=x_intial
    probability= intial_probability

    for i in 1:s.burning_time 
        x_point,probability=mcmcaccept(rng,probability,s)
        bunrning_index+=1
    end 

    for n_event in eachindex(A)
    #The matrix A has already the good dimension (N ,3 )
    x_point,probability=mcmcaccept(rng,probability,s)
    idx=1 
    A[n_event][1,:].=x_point
    
    N_n=size(A[i])[1]
    
    countingIndex=0
    while idx<N_n || countingIndex>1000 #limit for number of tries to find a new point
        #@show countingIndex
        x_point,probability=mcmcaccept(rng,probability,s)
        countingIndex+=1
        min_d=s.d_min
        for i in 1:idx 
            distance= sqrt(sum(x->(x[1]-x[2])^2,zip(view(A[n_event],i,:),x_point))) 
            min_d= min(distance,min_d)
        end
        if(min_d>= s.d_min)
            A[n_event][idx+1,:].=x_point
            idx=idx+1
        end
    end
    end 

    return A
end

