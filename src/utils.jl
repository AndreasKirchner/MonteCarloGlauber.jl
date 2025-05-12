"""
    center_of_mass(con; Nr=100, Nth=50) 

Calculate the center of mass for a given participant configuration.

# Arguments
- `con::Participant`: The participant configuration.
- `Nr`: Number of radial divisions (default: 100).
- `Nth`: Number of angular divisions (default: 50).

# Returns
- `SVector{3}`: The multiplicity ad x and y as static vector.
"""
function center_of_mass(con::T;Nr=100, Nth=50) where {T<:Participant}
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


"""
    multiplicity(con; Nr=100, Nth=50) 

Calculate the multiplicity for a given participant configuration.

# Arguments
- `con::Participant`: The participant configuration.
- `Nr`: Number of radial divisions (default: 100).
- `Nth`: Number of angular divisions (default: 50).

# Returns
- The multiplicity ad x and y as static vector.
"""
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
        con(x,y)*r*δr*δθ
    end 
end



"""
    centralities_selection(events; threaded=true) 

Selects centralities from a vector of events.

# Arguments
- `events::Vector{T}`: A vector of events, where each event is of type `T` which is a subtype of `Participant`.
- `threaded::Bool`: A keyword argument that specifies whether the operation should be threaded. Defaults to `true`.

# Returns
- The function returns the selected centralities from the given events.

# Example
"""
function centralities_selection_CoM(events::Vector{T},bins ;Threaded=true) where {T<:Participant}
    
    if Threaded 
        multi=tmap(events)   do x 
            center_of_mass(x)
        end
        mult=tmap(m->m[1],multi)
        #com1=tmap(m->m[2],multi)
        #com2=tmap(m->m[3],multi)
    else 
        multi=map(events)   do x 
            center_of_mass(x)
        end
        mult=map(m->m[1],multi)
        #com1=map(m->m[2],multi)
        #com2=map(m->m[3],multi)
    end 

    event_perm=sortperm(mult,rev=true)
    multiplicity_sorted=mult[event_perm]
    CoM_sorted=multi[event_perm]
    events_sorted=events[event_perm]

 
    min_bias=length(events)
    n_event_per_bin=Int(min_bias ÷ 100)
    realBinVals=n_event_per_bin*bins

    centrality_borders=map(x->x[1],multiplicity_sorted[1:n_event_per_bin:min_bias])

    return split_vector_by_indices(events_sorted,realBinVals), split_vector_by_indices(CoM_sorted,realBinVals)

    #return centrality_borders
end

function centralities_selection(events::Vector{T} ;threaded=true) where {T<:Participant}
    
    if threaded 
        mult=tmap(events)   do x 
            multiplicity(x)
        end
    else 
        mult=map(events)   do x 
            multiplicity(x)
        end
    end 

    event_perm=sortperm(mult,rev=true)
    multiplicity_sorted=mult[event_perm]
 
    min_bias=length(events)
    n_event_per_bin=min_bias ÷ 100

    centrality_borders=map(x->x[1],multiplicity_sorted[1:n_event_per_bin:min_bias])
    return centrality_borders
end

function construct_trento_names(part;extensionString="dat",mMode="2")
    proj1=projectile_dictionary(part.nucl1.N_nucleon)
    proj2=projectile_dictionary(part.nucl2.N_nucleon)
    w=convert_float_to_string(string(part.sub_nucleon_width))
    sigma_NN=convert_float_to_string(substring(string(part.total_cross_section),1,5))
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

substring(str, start, stop) = str[nextind(str, 0, start):nextind(str, 0, stop)]

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
    δr=last(r)/length(r)
    δθ=step
    mean(enumerate(conf_part)) do (i,conf)#for i in 1:length(conf_part)#mean(conf_part) do conf
        exp(1/(2pi)* sum(0:step:2pi) do θ
        sum(eachindex(r)) do iᵣ
            # exp(im\theta) h(r,theta) = h cos + ih sin = (re h+i im h)cos + (i re- im h )sin  
            #-f(norm*conf(r[iᵣ]*cos(θ)-CoMList[i][2]/CoMList[i][1],r[iᵣ]*sin(θ)-CoMList[i][3]/CoMList[i][1]))*h[iᵣ]*δθ#(h[iᵣ,1]*cos(θ*m)+h[iᵣ,2]*sin(θ*m))
            -conf(r[iᵣ]*cos(θ)-CoMList[i][2]/CoMList[i][1],r[iᵣ]*sin(θ)-CoMList[i][3]/CoMList[i][1])*h[iᵣ]*δθ
        end
        end)
    end 
    #return sumval/length(conf_part)
end

function generatingfunction_h(conf_part::Vector{T},h,r,step,CoMList,m,norm) where {T<:Participant}
    #sumval=0.0
    δr=last(r)/length(r)
    δθ=step
    mean(enumerate(conf_part)) do (i,conf)#for i in 1:length(conf_part)#mean(conf_part) do conf
        exp(1/(2pi)* sum(0:step:2pi) do θ
        sum(eachindex(r)) do iᵣ
            # exp(im\theta) h(r,theta) = h cos + ih sin = (re h+i im h)cos + (i re- im h )sin  
            -norm*conf(r[iᵣ]*cos(θ)-CoMList[i][2]/CoMList[i][1],r[iᵣ]*sin(θ)-CoMList[i][3]/CoMList[i][1])*(h[iᵣ,1]*cos(θ*m)+h[iᵣ,2]*sin(θ*m))*δθ
        end
        end)
    end 
end


function batched_events(Projectile1,Projectile2,w,k,p,sqrtS,bins;minBiasEvents=1000000,threaded=true)
    participants=Participants(Projectile1,Projectile2,w,x,k,p)

    if threaded
        event=rand(threaded(participants),minBiasEvents)
    else
        event=rand(participants,minBiasEvents)
    end

    return split_vector_by_indices(participants_Sort,realBinVals), split_vector_by_indices(CoM_mult_Sort,realBinVals)
end




function generate_background(f,norm,batches,CoM;r_grid=0:1:10,step=2pi/10)
    #batches, CoM=batched_events(Projectile1,Projectile2,w,k,p,sqrtS,bins;minBiasEvents=minBiasEvents)
    finalRes=tmap(1:length(batches)-1) do i
        #DifferentiationInterface.gradient(x->log(generatingfunction(batches[i],x,r_grid,step,CoM[i],f,norm)), AutoForwardDiff(), zeros(length(r_grid),2))
        DifferentiationInterface.gradient(x->log(generatingfunction(batches[i],x,r_grid,step,CoM[i])), AutoForwardDiff(), zeros(length(r_grid),2))
    end
    return map(x->f.(-norm .*x[:,1]),finalRes)
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


 function FT(mat)
    return mat[1,1]+mat[2,2]-im*(mat[1,2]+mat[2,1]) #check sign of complex part here
end


function generate_2ptfct(norm,batches, CoM,mList;r_grid=0:1:10,step=2pi/50)
    finalRes=stack(tmap(1:length(batches)-1) do i
        map(m->DifferentiationInterface.hessian(x->log(generatingfunction_h(batches[i],x,r_grid,step,CoM[i],m,norm)), AutoForwardDiff(), zeros(length(r_grid),2)),mList)
    end)
    hessianTransform=map(m->map(cc->reshape(finalRes[m,cc],(length(r_grid),2,length(r_grid),2)),1:length(batches)-1),1:length(mList))
    rIndeces=Iterators.product(1:length(r_grid),1:length(r_grid))
    twoPtFct=map(m->map(cc->map(rr->FT(hessianTransform[m][cc][rr[1],:,rr[2],:]),rIndeces),1:length(batches)-1),1:length(mList))
    return twoPtFct
end




function generate_bg_two_pt_fct(f,delta_factor,norm,Projectile1,Projectile2,w,k,p,sqrtS,bins,mList;minBiasEvents=1000000,r_grid=0:1:10,step=2pi/20,Threaded=true,n_ext_Grid=0,nFields=10)
    #batches, CoM=batched_events(Projectile1,Projectile2,w,k,p,sqrtS,bins;minBiasEvents=minBiasEvents)
    if (length(bins)+1)*100>minBiasEvents
       error("Not enough events for number of bins, increase minBiasEvents")
    end
    participants=Participants(Projectile1,Projectile2,w,sqrtS,k,p)
    #if threaded
        events=rand(threaded(participants),minBiasEvents)
    #else
    #    events=rand(participants,minBiasEvents)
    #end
    batches, CoM=centralities_selection_CoM(events,bins;Threaded=Threaded)
    bg=generate_background(f,norm,batches,CoM,r_grid=r_grid,step=step)
    twoPtFct_entropy=generate_2ptfct(norm,batches, CoM,mList;r_grid=r_grid,step=step)
    #@show  size(twoPtFct_entropy[1][1][1]),size(twoPtFct_entropy[1]),size(twoPtFct_entropy[1]),size(twoPtFct_entropy)
    twoPtFct=map(m->map(cc->map(r1->map(r2->twoPtFct_entropy[m][cc][r1,r2]*delta_factor(bg[cc][r1])*delta_factor(bg[cc][r2]),1:length(r_grid)),1:length(r_grid)),1:length(bg)),1:length(mList))
    #twoPtFct=map(m->map(cc->map(r1->map(r2->twoPtFct_entropy[m][cc][r1][r2],1:length(r_grid)),1:length(r_grid)),1:length(bg)),1:length(mList))
    nGrid=max(length(r_grid),n_ext_Grid)
    finalCorr=zeros(2,nFields,length(mList),nGrid,nGrid)
    for cc in 1:length(bg)
        for m in 1:length(mList)
            for r1 in 1:length(r_grid)
                for r2 in 1:length(r_grid)
                    finalCorr[1,1,m,r1,r2]=real(twoPtFct[m][cc][r1][r2])
                end
            end
        end
    end
    return bg,finalCorr
end

function generate_bg(f,norm,Projectile1,Projectile2,w,k,p,sqrtS,bins;minBiasEvents=1000000,r_grid=0:1:10,step=2pi/10,Threaded=true)
    #batches, CoM=batched_events(Projectile1,Projectile2,w,k,p,sqrtS,bins;minBiasEvents=minBiasEvents)
    participants=Participants(Projectile1,Projectile2,w,sqrtS,k,p)
    #if threaded
        events=rand(threaded(participants),minBiasEvents)
    #else
    #    events=rand(participants,minBiasEvents)
    #end
    batches, CoM=centralities_selection_CoM(events,bins;Threaded=Threaded)
    bg=generate_background(f,norm,batches,CoM,r_grid=r_grid,step=step)
    return bg
end


function save_bg_two_pt_fct(f,delta_factor,norm,Projectile1,Projectile2,w,k,p,sqrtS,bins,mList;minBiasEvents=1000000,r_grid=0:1:10,step=2pi/20,Threaded=true)
    participants=Participants(Projectile1,Projectile2,w,sqrtS,k,p)
    #if threaded
        events=rand(threaded(participants),minBiasEvents)
    #else
     #   events=rand(participants,minBiasEvents)
    #end
    batches, CoM=centralities_selection_CoM(events,bins;Threaded=Threaded)
    bg=generate_background(f,norm,batches,CoM,r_grid=r_grid,step=step)
    twoPtFct_entropy=generate_2ptfct(norm,batches, CoM,mList;r_grid=r_grid,step=step)
    twoPtFct=map(m->map(cc->map(r1->map(r2->twoPtFct_entropy[m][cc][r1,r2]*delta_factor(bg[cc][r1])*delta_factor(bg[cc][r2]),1:length(r_grid)),1:length(r_grid)),1:length(bg)),1:length(mList))
    bgString,twoptString=construct_trento_names(participants;extensionString="dat",mMode="2")
    writedlm(bgString,bg)
    for i in 1:length(mList)
        bgString,twoptString=construct_trento_names(participants;extensionString="dat",mMode=mList[i])
        writedlm(twoptString,real.(twoPtFct[i]))
    end

end