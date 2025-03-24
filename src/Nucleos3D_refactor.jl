
struct NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,B,samp} <:Sampleable{ArrayLikeVariate{2},Continuous}
    N_nucleon::Int64
    α::T
    R::R
    ρ₀::C
    w::W
    beta2::D1
    beta3::D2
    beta4::D3
    gamma::D4
    burning_time::B
    sampler::samp
    d_min::Float64
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


dimension(::NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,B,samp}) where {T,R,C,W,D1,D2,D3,D4,B,samp} = 2
nucleon(s::NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,B,samp}) where {T,R,C,W,D1,D2,D3,D4,B,samp} = s.N_nucleon

Base.size(s::NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,B,samp}) where {T,R,C,W,D1,D2,D3,D4,B,samp} = (s.N_nucleon,2)
function Base.eltype(s::NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,B,samp}) where {T,R,C,W,D1,D2,D3,D4,B,samp} 
    promote_type(T,R,C,W,D1,D2,D3,D4)
end 




function NucleiWoodSaxon3D(rng::AbstractRNG,N_nucleon,α,R,ρ₀,w,beta2,beta3,beta4,gamma,burning_time,dmin)

    #draw=(Uniform(-2R,2R),Uniform(-2R,2R),Uniform(-2R,2R))
    draw=(Normal(0,R),Normal(0,R),Normal(0,R))
    
    params=(α,R,ρ₀,w,beta2,beta3,beta4,gamma)

    sampler=Metropolis_Hastings(rng,density_WS_deformed,params,draw,burning_time)

    NucleiWoodSaxon3D(N_nucleon,α,R,ρ₀,w,beta2,beta3,beta4,gamma,burning_time,sampler,dmin)

end

function NucleiWoodSaxon3D(N_nucleon,α,R,ρ₀,w,beta2,beta3,beta4,gamma,dmin)
    rng=Random.default_rng()
    burning_time=1000
    NucleiWoodSaxon3D(rng,N_nucleon,α,R,ρ₀,w,beta2,beta3,beta4,gamma,burning_time,dmin)
end 

"""
    NucleiWoodSaxon3D(N_nucleon,α,R,ρ₀,w,beta2,beta3,beta4,gamma)

Create a 3D distribution of nucleons in a nucleus with a Wood-Saxon density profile. The density profile is deformed by the parameters `beta2`, `beta3`, `beta4`, and `gamma`. The number of nucleons is `N_nucleon`, the radius is `R`, the diffuseness is `α`, the central density is `ρ₀`, and the skin depth is `w`.
"""
function NucleiWoodSaxon3D(N_nucleon,α,R,ρ₀,w,beta2,beta3,beta4,gamma)
    rng=Random.default_rng()
    burning_time=100000
    dmin=0.
    NucleiWoodSaxon3D(rng,N_nucleon,α,R,ρ₀,w,beta2,beta3,beta4,gamma,burning_time,dmin)
end 



function Base.copy(s::NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,B,samp}) where {T,R,C,W,D1,D2,D3,D4,B,samp}
    NucleiWoodSaxon3D(s.N_nucleon,s.α,s.R,s.ρ₀,s.w,s.beta2,s.beta3,s.beta4,s.gamma,s.d_min)
end



function Distributions._rand!(rng::AbstractRNG, s::NucleiWoodSaxon3D{T,R,C,W,D1,D2,D3,D4,B,samp}, x::DenseMatrix{<:Real}) where {T,R,C,W,D1,D2,D3,D4,B,samp}
    
    if s.d_min==0.
    rotmatrix=rand(rng,RotMatrix{3})
        for i in Base.axes(x,1)
            xpos=rand(rng,s.sampler)
            rot=rotmatrix*SVector{3}(xpos)
            for j in Base.axes(x,2)
                x[i,j] =rot[j]
            end
        end 

    return x
    end 

    if s.d_min>0

        rotmatrix=rand(rng,RotMatrix{3})
        for i in Base.axes(x,1)
            xpos=rand(rng,s.sampler)
            rot=rotmatrix*SVector{3}(xpos)
            for j in Base.axes(x,2)
                x[i,j] =rot[j]
            end
        end 

    return x

    end 

end 