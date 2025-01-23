#This file has all the functions to run the event by event initial conditions and to generate the fluidum intiail conditions
using Distributions
using Random
#using QuadGK
#using Interpolations
using SpecialFunctions
using SimpleNonlinearSolve
using StatsBase
using StaticArrays
using LinearAlgebra
using Plots
#using BenchmarkTools
using OhMyThreads
using Rotations
#using Profile
#using PProf
using JLD2
using DifferentiationInterface
using IterTools
using LaTeXStrings

include("nucleos_distribution.jl")
include("participants_distribution.jl")

function onlyecc_properly(con::T,m;Nr=50, Nth=10) where {T<:Participant}
    R1=con.R1
    R2=con.R2
    Rmax=3*(R1+R2)
    δr=Rmax/Nr
    δθ=2pi/Nth
    r_range=range(δr/2,Rmax-δr/2,Nr)
    theta_range=range(0,2pi,Nth)
    mult,x0,y0=sum(Iterators.product(enumerate(r_range),enumerate(theta_range))) do ((ir,r),(iθ, θ))
        s,c=sincos(θ) 
        x=r*c
        y=r*s
        con(x,y)*r*δr*δθ .*SVector{3}(1,x,y)
    end
    num=sum(Iterators.product(enumerate(r_range),enumerate(theta_range))) do ((ir,r),(iθ, θ))
        s,c=sincos(θ) 
        x=r*c
        y=r*s
        con(x-x0/mult,y-y0/mult)*r*δr*δθ*r .^m .*exp.(Ref(im) .*m .*Ref(θ))
    end
    denom=sum(Iterators.product(enumerate(r_range),enumerate(theta_range))) do ((ir,r),(iθ, θ))
        s,c=sincos(θ) 
        x=r*c
        y=r*s
        con(x-x0/mult,y-y0/mult)*r*δr*δθ*r .^m
    end 
    return con.b,con.n_coll,mult,hypot(x0/mult,y0/mult),atan(y0,x0),num./denom ...
end

"""
    trento_cmd_line(nucleus1::String, nucleus2::String, w::Float64, k::Float64, p::Float64, sqrtS::Float64, nEvents::Int64, [mList]::Vec[Float64],[b]::Float64)::Vec[Float64]
    
Run TrenTo and return integrated quantities, not the profiles. The output is impact parameter, Ncoll, Center of mass radius, Center of mass angle, eccentricities corresponding to mList input

### Input

    - `nucleus1::String` -- String representing the nucleus 1. For example, "Pb"
    - `nucleus2::String` -- String representing the nucleus 2. For example, "Pb"
    - `w::Float64` -- Float64 nucleon width
    - `k::Float64` -- Float64 amount of fluctuations
    - `p::Float64` -- Float64 power of the generalized mean
    - `sqrtS::Float64` -- Float64 center of mass energy, in GeV, the sigma_NN value is obtained via a fit obtained from 2306.08665
    - `nEvents::Int64` -- Int64 number of events to generate
    - `mList::Vec[Float64]` -- Vector of Float64 values of the eccentricities to calculate
    - `b::Float64` -- Float64 impact parameter. For a number all events are generated with the same impact parameter
    - `b::Tuple{Float64,Float64}` -- Tuple of Float64 values of the impact parameter range. 

### Output

The output is a vector containing impact parameter, Ncoll, Center of mass radius, Center of mass angle, eccentricities corresponding to mList input

### Notes

nothjiong
"""
function trento_cmd_line(Projectile1,Projectile2,w,k,p,sqrtS,nevents;mList=[2,3,4,5],b=-1)
    sigma_NN=cross_section_from_energy(sqrtS)
    #if b<-1
        participants=Participants(Projectile1,Projectile2,w,sigma_NN,k,p)
   # else
    #    participants=Participants(Projectile1,Projectile2,w,b,sigma_NN,k,p)
    #end
    events=rand(participants,nevents)
    return stack(map(i->onlyecc_properly(events[i],mList),1:nevents))
#output multiplicity, center of mass, reaction plane angle, eccentricities, npart, ncoll, b
end



function cross_section_from_energy(sqrtS) #returns the cross section for sqrtS in GeV
    return (0.073491826*log(sqrtS)-.19313457)*log(sqrtS)+3.123737545
end

function construct_trento_names(part;extensionString="dat",mMode="2")
    proj1=projectile_dictionary(part.nucl1.N_nucleon)
    proj2=projectile_dictionary(part.nucl2.N_nucleon)
    w=convert_float_to_string(string(part.sub_nucleon_width))
    sigma_NN=convert_float_to_string(string(part.total_cross_section))
    p=convert_float_to_string(string(part.p))
    k=convert_float_to_string(string(part.shape_parameter))
    m=string(mMode)
    bgString="trento_BG_"*proj1*"_"*proj2*"_w_"*w*"_sNN_"*sigma_NN*"_p_"*p*"_k_"*k*"."*extensionString
    twoptString="trento_two_pt_fct_m_"*m*"_"*proj1*"_"*proj2*"_w_"*w*"_sNN_"*sigma_NN*"_p_"*p*"_k_"*k*"."*extensionString
    return bgString,twoptString
end

function convert_float_to_string(f)
    return replace(string(f), "." => "_")
end

function check_for_config(part;extensionString="dat",mMode="2",path="")
    bgString,twoptString=construct_trento_names(part;extensionString=extensionString,mMode=mMode)
    bgString=path*bgString
    twoptString=path*twoptString
    return isfile(bgString), isfile(twoptString)
end


function projectile_dictionary(massNumber)
    dict=Dict(208=>"Pb",197=>"Au",238=>"U",63=>"Cu",4=>"He",16=>"O",129=>"Xe",84=>"Kr",40=>"Ar",20=>"Ne",14=>"N",12=>"C",1=>"H")
    return dict[massNumber]
end


function multiplicity(con::T;Nr=100, Nth=50) where {T<:Participant}
    R1=con.R1
    R2=con.R2
    Rmax=3*(R1+R2)
    δr=Rmax/Nr
    δθ=2pi/Nth
    r_range=range(δr/2,Rmax-δr/2,Nr)
    theta_range=range(0,2pi,Nth)
    sum(Iterators.product(enumerate(r_range),enumerate(theta_range))) do ((ir,r),(iθ, θ))
        s,c=sincos(θ) 
        x=r*c
        y=r*s
        SVector{3}(1,x,y) .*con(x,y)*r*δr*δθ
    end 
end

function multiplicity_threads(con::Vector{T};Nr=500, Nth=50) where {T<:Participant}

    tmap(con) do rr
        multiplicity(rr;Nr=Nr, Nth=Nth)
    end
end 


function define_centralities(Projectile1,Projectile2,w,k,p,sqrtS,bins;minBiasEvents=1000000)
    x=cross_section_from_energy(sqrtS)
    participants=rand(Participants(Projectile1,Projectile2,w,x,k,p),minBiasEvents)
    CoM_mult=multiplicity_threads(participants)
    eventPerm=sortperm(CoM_mult,rev=true)
    CoM_mult_Sort=CoM_mult[eventPerm]
    #participants_Sort=participants[eventPerm]
    Nev_per_bin=Int(minBiasEvents/100)
    centBorders=map(x->x[1],CoM_mult_Sort[1:Nev_per_bin:minBiasEvents])
    return centBorders[bins]
end

function define_centralities(participants,bins)
    CoM_mult=multiplicity_threads(participants)
    eventPerm=sortperm(CoM_mult,rev=true)
    CoM_mult_Sort=CoM_mult[eventPerm]
    #participants_Sort=participants[eventPerm]
    Nev_per_bin=Int(minBiasEvents/100)
    centBorders=map(x->x[1],CoM_mult_Sort[1:Nev_per_bin:length(CoM_mult_Sort)])
    return centBorders[bins]
end

function batched_events(Projectile1,Projectile2,w,k,p,sqrtS,bins;minBiasEvents=1000000)
    x=cross_section_from_energy(sqrtS)
    participants=rand(Participants(Projectile1,Projectile2,w,x,k,p),minBiasEvents)
    CoM_mult=multiplicity_threads(participants)
    eventPerm=sortperm(CoM_mult,rev=true)
    CoM_mult_Sort=CoM_mult[eventPerm]
    participants_Sort=participants[eventPerm]
    Nev_per_bin=Int(minBiasEvents/100)
    realBinVals=Nev_per_bin*bins
    #centBorders=map(x->x[1],CoM_mult_Sort[1:Nev_per_bin:minBiasEvents])
    return split_vector_by_indices(participants_Sort,realBinVals), split_vector_by_indices(CoM_mult_Sort,realBinVals)
end


function split_vector_by_indices(vector, indices)
    chunks = []
    start_idx = 1
    for idx in indices
        push!(chunks, vector[start_idx:idx])
        start_idx = idx + 1
    end
    if start_idx <= length(vector)
        push!(chunks, vector[start_idx:end])
    end
    return chunks
end

function generatingfunction(conf_part::Vector{T},h,r,step,CoMList) where {T<:Participant}
    #sumval=0.0
    mean(enumerate(conf_part)) do (i,conf)#for i in 1:length(conf_part)#mean(conf_part) do conf
        exp(1/(2pi)* sum(0:step:2pi) do θ
        sum(eachindex(r)) do iᵣ
            # exp(im\theta) h(r,theta) = h cos + ih sin = (re h+i im h)cos + (i re- im h )sin  
            -conf(r[iᵣ]*cos(θ)-CoMList[i][2]/CoMList[i][1],r[iᵣ]*sin(θ)-CoMList[i][3]/CoMList[i][1])*h[iᵣ]#(h[iᵣ,1]*cos(θ*m)+h[iᵣ,2]*sin(θ*m))
        end
        end)
    end 
    #return sumval/length(conf_part)
end

function generatingfunction(conf_part::Vector{T},h,r,step,CoMList,f,norm) where {T<:Participant}
    #sumval=0.0
    mean(enumerate(conf_part)) do (i,conf)#for i in 1:length(conf_part)#mean(conf_part) do conf
        exp(1/(2pi)* sum(0:step:2pi) do θ
        sum(eachindex(r)) do iᵣ
            # exp(im\theta) h(r,theta) = h cos + ih sin = (re h+i im h)cos + (i re- im h )sin  
            -f(norm*conf(r[iᵣ]*cos(θ)-CoMList[i][2]/CoMList[i][1],r[iᵣ]*sin(θ)-CoMList[i][3]/CoMList[i][1]))*h[iᵣ]#(h[iᵣ,1]*cos(θ*m)+h[iᵣ,2]*sin(θ*m))
        end
        end)
    end 
    #return sumval/length(conf_part)
end

function generatingfunction_h(conf_part::Vector{T},h,r,step,CoMList,m) where {T<:Participant}
    #sumval=0.0
    mean(enumerate(conf_part)) do (i,conf)#for i in 1:length(conf_part)#mean(conf_part) do conf
        exp(1/(2pi)* sum(0:step:2pi) do θ
        sum(eachindex(r)) do iᵣ
            # exp(im\theta) h(r,theta) = h cos + ih sin = (re h+i im h)cos + (i re- im h )sin  
            -conf(r[iᵣ]*cos(θ)-CoMList[i][2]/CoMList[i][1],r[iᵣ]*sin(θ)-CoMList[i][3]/CoMList[i][1])*(h[iᵣ,1]*cos(θ*m)+h[iᵣ,2]*sin(θ*m))
        end
        end)
    end 
    #return sumval/length(conf_part)
end

function generatingfunction_h(conf_part::Vector{T},h,r,step,CoMList,m,f,norm) where {T<:Participant}
    #sumval=0.0
    mean(enumerate(conf_part)) do (i,conf)#for i in 1:length(conf_part)#mean(conf_part) do conf
        exp(1/(2pi)* sum(0:step:2pi) do θ
        sum(eachindex(r)) do iᵣ
            # exp(im\theta) h(r,theta) = h cos + ih sin = (re h+i im h)cos + (i re- im h )sin  
            -f(norm*conf(r[iᵣ]*cos(θ)-CoMList[i][2]/CoMList[i][1],r[iᵣ]*sin(θ)-CoMList[i][3]/CoMList[i][1]))*(h[iᵣ,1]*cos(θ*m)+h[iᵣ,2]*sin(θ*m))
        end
        end)
    end 
    #return sumval/length(conf_part)
end

function generate_background(Projectile1,Projectile2,w,k,p,sqrtS,bins;minBiasEvents=1000000,r_grid=0:1:10,step=2pi/10)
    batches, CoM=batched_events(Projectile1,Projectile2,w,k,p,sqrtS,bins;minBiasEvents=minBiasEvents)
    finalRes=tmap(1:length(batches)-1) do i
        DifferentiationInterface.gradient(x->log(generatingfunction(batches[i],x,r_grid,step,CoM[i])), AutoForwardDiff(), zeros(length(r_grid),2))       
    end
    return -map(x->x[:,1],finalRes)
end

function generate_background(f,norm,Projectile1,Projectile2,w,k,p,sqrtS,bins;minBiasEvents=1000000,r_grid=0:1:10,step=2pi/10)
    batches, CoM=batched_events(Projectile1,Projectile2,w,k,p,sqrtS,bins;minBiasEvents=minBiasEvents)
    finalRes=tmap(1:length(batches)-1) do i
        DifferentiationInterface.gradient(x->log(generatingfunction(batches[i],x,r_grid,step,CoM[i],f,norm)), AutoForwardDiff(), zeros(length(r_grid),2))       
    end
    return -map(x->x[:,1],finalRes)
end

function generate_background(f,norm,batches,CoM;r_grid=0:1:10,step=2pi/10)
    #batches, CoM=batched_events(Projectile1,Projectile2,w,k,p,sqrtS,bins;minBiasEvents=minBiasEvents)
    finalRes=tmap(1:length(batches)-1) do i
        DifferentiationInterface.gradient(x->log(generatingfunction(batches[i],x,r_grid,step,CoM[i],f,norm)), AutoForwardDiff(), zeros(length(r_grid),2))
    end
    return -map(x->x[:,1],finalRes)
end



function FT(mat)
    return mat[1,1]+mat[2,2]-im*(mat[1,2]+mat[2,1]) #check sign of complex part here
end

function generate_2ptfct(Projectile1,Projectile2,w,k,p,sqrtS,bins,mList;minBiasEvents=1000000,r_grid=0:1:10,step=2pi/50)
    batches, CoM=batched_events(Projectile1,Projectile2,w,k,p,sqrtS,bins;minBiasEvents=minBiasEvents)
    finalRes=stack(tmap(1:length(batches)-1) do i
        map(m->DifferentiationInterface.hessian(x->log(generatingfunction_h(batches[i],x,r_grid,step,CoM[i],m)), AutoForwardDiff(), zeros(length(r_grid),2)),mList)
    end)
    hessianTransform=map(m->map(cc->reshape(finalRes[m,cc],(length(r_grid),2,length(r_grid),2)),1:length(bins)),1:length(mList))
    rIndeces=Iterators.product(1:length(r_grid),1:length(r_grid))
    twoPtFct=map(m->map(cc->map(rr->FT(hessianTransform[m][cc][rr[1],:,rr[2],:]),rIndeces),1:length(bins)),1:length(mList))
    return twoPtFct
end

function generate_2ptfct(f,norm,Projectile1,Projectile2,w,k,p,sqrtS,bins,mList;minBiasEvents=1000000,r_grid=0:1:10,step=2pi/50)
    batches, CoM=batched_events(Projectile1,Projectile2,w,k,p,sqrtS,bins;minBiasEvents=minBiasEvents)
    finalRes=stack(tmap(1:length(batches)-1) do i
        map(m->DifferentiationInterface.hessian(x->log(generatingfunction_h(batches[i],x,r_grid,step,CoM[i],m,f,norm)), AutoForwardDiff(), zeros(length(r_grid),2)),mList)
    end)
    hessianTransform=map(m->map(cc->reshape(finalRes[m,cc],(length(r_grid),2,length(r_grid),2)),1:length(bins)),1:length(mList))
    rIndeces=Iterators.product(1:length(r_grid),1:length(r_grid))
    twoPtFct=map(m->map(cc->map(rr->FT(hessianTransform[m][cc][rr[1],:,rr[2],:]),rIndeces),1:length(bins)),1:length(mList))
    return twoPtFct
end

function generate_2ptfct(f,norm,batches, CoM,mList;r_grid=0:1:10,step=2pi/50)
    finalRes=stack(tmap(1:length(batches)-1) do i
        map(m->DifferentiationInterface.hessian(x->log(generatingfunction_h(batches[i],x,r_grid,step,CoM[i],m,f,norm)), AutoForwardDiff(), zeros(length(r_grid),2)),mList)
    end)
    hessianTransform=map(m->map(cc->reshape(finalRes[m,cc],(length(r_grid),2,length(r_grid),2)),1:length(batches)-1),1:length(mList))
    rIndeces=Iterators.product(1:length(r_grid),1:length(r_grid))
    twoPtFct=map(m->map(cc->map(rr->FT(hessianTransform[m][cc][rr[1],:,rr[2],:]),rIndeces),1:length(batches)-1),1:length(mList))
    return twoPtFct
end

struct InverseFunction{N,F}
    fun::F
end 

function InverseFunction(f)
    InverseFunction{1,typeof(f)}(f)
end

function (func::InverseFunction{1,F})(x;u0=0.1*one(x)) where F 
    
   # function f(u,p)
   #     func.fun(u)-x
   # end 
    f= let x=x
        (u,p)->func.fun(u)-x 
    end
    prob = NonlinearProblem{false}(f,u0) 
    
    return solve(prob, SimpleNewtonRaphson()).u
end

function generate_bg_two_pt_fct(f,norm,Projectile1,Projectile2,w,k,p,sqrtS,bins,m;minBiasEvents=1000000,r_grid=0:1:10,step=2pi/20)
    batches, CoM=batched_events(Projectile1,Projectile2,w,k,p,sqrtS,bins;minBiasEvents=minBiasEvents)
    bg=generate_background(f,norm,batches,CoM,r_grid=r_grid,step=step)
    twoPtFct=generate_2ptfct(f,norm,batches, CoM,m;r_grid=r_grid,step=step)
    return bg,twoPtFct
end



####minimal working example

#define EoS and its inverse
energy(T)=T^4
Te4=InverseFunction(energy)
# get the background and two point function
bg,twpt=generate_bg_two_pt_fct(Te4,100,Lead(),Lead(),1.2,3.0,0.5,2760,[10,20],[2,3,4];minBiasEvents=500,r_grid=0:1:20)

#plot the background
plot(bg[1])
#plot the two point function
heatmap(real.(twpt[1][1]))
heatmap(imag.(twpt[1][1]))