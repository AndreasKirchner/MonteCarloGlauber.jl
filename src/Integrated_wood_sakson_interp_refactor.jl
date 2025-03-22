
struct IntegratedWoodSaxonInterp{T,R,C,W,B,Inte,samp} <: Sampleable{ArrayLikeVariate{2},Continuous}
    N_nucleon::Int64
    α::T
    R::R
    ρ₀::C
    w::W
    burning_time::B
    retol::Float64
    npoint::Int
    interpolation::Inte
    sampler::samp
    rotated::Bool
end

dimension(::IntegratedWoodSaxonInterp{T,R,C,W,B,Inte,samp}) where {T,R,C,W,B,Inte,samp} = 2
nucleon(s::IntegratedWoodSaxonInterp{T,R,C,W,B,Inte,samp}) where {T,R,C,W,B,Inte,samp} = s.N_nucleon

Base.size(s::IntegratedWoodSaxonInterp{T,R,C,W,B,Inte,samp}) where {T,R,C,W,B,Inte,samp} = (s.N_nucleon,2)
function Base.eltype(s::IntegratedWoodSaxonInterp{T,R,C,W,B,Inte,samp}) where {T,R,C,W,B,Inte,samp}
    promote_type(T,R,C,W)
end 


function IntegratedWoodSaxonInterp(N_nucleon,α,R,ρ₀,w,burning_time,retol,Npoint,rng=Random.default_rng(),rotated=true)
    
    range=3*R
    
    A_x1 = -range:2*range/Npoint:range
    A_x2 = -range:2*range/Npoint:range
    
    A = [integrate_charge(x1,x2,α,R,ρ₀,w,retol) for x1 in A_x1, x2 in A_x2]
    itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
    
    sitp = Interpolations.scale(itp, A_x1, A_x2)

    draw=(truncated(Normal(0,R), -3R, 3R),truncated(Normal(0,R), -3R, 3R))
    
    params=()
    sampler=Metropolis_Hastings(rng, sitp,params,draw,burning_time)

    IntegratedWoodSaxonInterp(N_nucleon,α,R,ρ₀,w,burning_time,retol,Npoint,sitp,sampler,rotated)

end

"""
    IntegratedWoodSaxonInterp(N_nucleon,α,R,ρ₀,w,rng=Random.default_rng(),rotated=true)

Create a sampleable object that represents a Wood Saxon distribution with a given number of nucleons.
The distribution is deformed by the parameters `α`, `R`, `ρ₀`, and `w`. The number of nucleons is `N_nucleon`.
"""
function IntegratedWoodSaxonInterp(N_nucleon,α,R,ρ₀,w,rng=Random.default_rng(),rotated=true)
    burning_time=1000
    retol=1e-3
    Npoint=100
    IntegratedWoodSaxonInterp(N_nucleon,α,R,ρ₀,w,burning_time,retol,Npoint,rng,rotated)
end 

function Base.copy(s::IntegratedWoodSaxonInterp{T,R,C,W,B,Inte,samp}) where {T,R,C,W,B,Inte,samp}
    IntegratedWoodSaxonInterp(s.N_nucleon,s.α,s.R,s.ρ₀,s.w)
end



function Distributions._rand!(rng::AbstractRNG, s::IntegratedWoodSaxonInterp{T,R,C,W,B,Inte,samp}, x::DenseMatrix{<:Real}) where {T,R,C,W,B,Inte,samp}
    if s.rotated
        rotmatrix=rand(rng,RotMatrix{2})
        for i in Base.axes(x,1)
            xpos=rand(rng,s.sampler)
            rot=rotmatrix*SVector{2}(xpos)
            for j in Base.axes(x,2)
                x[i,j] =rot[j]
            end
        end 
    return x
    end 

    if !s.rotated
        
        for i in Base.axes(x,1)
            xpos=rand(rng,s.sampler)
            for j in Base.axes(x,2)
                x[i,j] =xpos[j]
            end
        end 
        return x
    end



end 