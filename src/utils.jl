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
        SVector{3}(1.,x,y) .*con(x,y)*r*δr*δθ
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
function centralities_selection_CoM(events::Vector{T},bins; Threaded=true ) where {T<:Participant}
    
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
   # multiplicity_sorted=mult[event_perm]
    CoM_sorted=multi[event_perm]
    events_sorted=events[event_perm]

 
    min_bias=length(events)
    n_event_per_bin=Int(min_bias ÷ 100)
    realBinVals=n_event_per_bin*bins

    #centrality_borders=map(x->x[1],multiplicity_sorted[1:n_event_per_bin:min_bias])

    #return split_vector_by_indices(events_sorted,realBinVals), split_vector_by_indices(CoM_sorted,realBinVals)
    batches= split_vector_by_indices(events_sorted,realBinVals)
    CoM= split_vector_by_indices(CoM_sorted,realBinVals)
    return batches, CoM

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

function construct_trento_names(part;extensionString="dat",mMode="2",cc="0-100",path="./")
    proj1=projectile_dictionary(part.nucl1.N_nucleon)
    proj2=projectile_dictionary(part.nucl2.N_nucleon)
    w=convert_float_to_string(string(part.sub_nucleon_width))
    sigma_NN=convert_float_to_string(substring(string(part.total_cross_section),1,5))
    p=convert_float_to_string(string(part.p))
    k=convert_float_to_string(string(part.shape_parameter))
    m=string(mMode)
    bin=string(cc)
    bgString=path*"trento_BG_"*proj1*"_"*proj2*"_w_"*w*"_sNN_"*sigma_NN*"_p_"*p*"_k_"*k*"_cc_"*bin*"."*extensionString
    twoptString=path*"trento_two_pt_fct_m_"*m*"_"*proj1*"_"*proj2*"_w_"*w*"_sNN_"*sigma_NN*"_p_"*p*"_k_"*k*"_cc_"*bin*"."*extensionString
    return bgString,twoptString
end


function convert_float_to_string(f)
    return replace(string(f), "." => "_")
end

substring(str, start, stop) = str[nextind(str, 0, start):nextind(str, 0, stop)]

function check_for_config(part;extensionString="dat",mMode="2",path="./",cc="0-100")
    bgString,twoptString=construct_trento_names(part;extensionString=extensionString,mMode=mMode,cc=cc)
    bgString=path*bgString
    twoptString=path*twoptString
    return isfile(bgString), isfile(twoptString)
end


function projectile_dictionary(massNumber)
    dict=Dict(208=>"Pb",197=>"Au",238=>"U",63=>"Cu",4=>"He",16=>"O",129=>"Xe",84=>"Kr",40=>"Ar",20=>"Ne",14=>"N",12=>"C",1=>"H")
    return dict[massNumber]
end

function split_vector_by_indices_old(vector, indices)
    chunks = empty(vector,eltype(vector))
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

function split_vector_by_indices(vector, indices)
    # Combine the start index (1) and all indices+1 as starting points,
    # and all indices and the end index as ending points.
    
    # 1. Prepare all the index ranges.
    # The start indices are 1 and all indices + 1
    starts = vcat(1, indices .+ 1)
    
    # The end indices are all indices and the end of the vector
    ends = vcat(indices, length(vector))

    # 2. Use a generator/comprehension to create the chunks.
    # This automatically infers the concrete type of the resulting vector.
    chunks = [vector[s:e] for (s, e) in zip(starts, ends) if s <= e]
    
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




function generate_background(f,norm,batches,CoM;r_grid=0.:1.:10.,step=2pi/10)
    #batches, CoM=batched_events(Projectile1,Projectile2,w,k,p,sqrtS,bins;minBiasEvents=minBiasEvents)
   # ex = DifferentiationInterface.gradient(x->log(generatingfunction(batches[1],x,r_grid,step,CoM[1])), AutoForwardDiff(), zeros(eltype(r_grid),length(r_grid),2))
   # finalRes = zeros(Float64,length(batches)-1, length(r_grid), 2)
    finalRes=map(1:length(batches)-1) do i
        grad=zeros(eltype(r_grid),length(r_grid))
        #DifferentiationInterface.gradient(x->log(generatingfunction(batches[i],x,r_grid,step,CoM[i],f,norm)), AutoForwardDiff(), zeros(length(r_grid),2))
        DifferentiationInterface.gradient!(x->log(generatingfunction(batches[i],x,r_grid,step,CoM[i])),grad, AutoForwardDiff(), zeros(eltype(r_grid),length(r_grid),2))
        return grad
    end
    return map(x->f.(-norm .*x[:,1]),finalRes)
end

function generate_background(f,norm,batches,CoM,cc;r_grid=0:1:10,step=2pi/10)
    #batches, CoM=batched_events(Projectile1,Projectile2,w,k,p,sqrtS,bins;minBiasEvents=minBiasEvents)
    finalRes=DifferentiationInterface.gradient(x->log(generatingfunction(batches[cc],x,r_grid,step,CoM[cc])), AutoForwardDiff(), zeros(length(r_grid),2))
    return f.(-norm .*finalRes[:,1])
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
    finalRes=stack(map(1:length(batches)-1) do i
        hess=zeros(eltype(r_grid),2*length(r_grid),2*length(r_grid))
        map(m->DifferentiationInterface.hessian!(x->log(generatingfunction_h(batches[i],x,r_grid,step,CoM[i],m,norm)),hess, AutoForwardDiff(), zeros(eltype(r_grid),length(r_grid),2)),mList)
        return hess
    end)####From here on it is not type stable anymore
    hessianTransform=map(m->map(cc->reshape(finalRes[m,cc],(length(r_grid),2,length(r_grid),2)),1:length(batches)-1),1:length(mList))
    rIndeces=Iterators.product(1:length(r_grid),1:length(r_grid))
    twoPtFct=map(m->map(cc->map(rr->FT(hessianTransform[m][cc][rr[1],:,rr[2],:]),rIndeces),1:length(batches)-1),1:length(mList))
    return twoPtFct
end

function generate_2ptfct_all(norm,batches,CoM,mList;r_grid=0:1:10,step=2pi/50)
    res=zeros(ComplexF64,length(mList),length(CoM)-1,length(r_grid),length(r_grid))
    map(1:length(batches)-1) do cc
        map(1:length(mList)) do m
            res[m,cc,:,:]=generate_2ptfct(norm,batches,CoM,mList[m],cc;r_grid=r_grid,step=step)
        end
    end
    return res
end

function generate_2ptfct(norm,batches,CoM,m::Int,cc::Int;r_grid=0:1:10,step=2pi/50)
    hess=zeros(eltype(r_grid),2*length(r_grid),2*length(r_grid))
    finalRes=DifferentiationInterface.hessian!(x->log(generatingfunction_h(batches[cc],x,r_grid,step,CoM[cc],m,norm)),hess, AutoForwardDiff(), zeros(eltype(r_grid),length(r_grid),2))
    hessianTransform=reshape(finalRes,(length(r_grid),2,length(r_grid),2))
    rIndeces=Iterators.product(1:length(r_grid),1:length(r_grid))
    twoPtFct=map(rr->FT(hessianTransform[rr[1],:,rr[2],:]),rIndeces)
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
    twoPtFct_entropy=generate_2ptfct_all(norm,batches, CoM,mList;r_grid=r_grid,step=step)
    #@show  size(twoPtFct_entropy[1][1][1]),size(twoPtFct_entropy[1]),size(twoPtFct_entropy[1]),size(twoPtFct_entropy)
    twoPtFct=map(m->map(cc->map(r1->map(r2->twoPtFct_entropy[m][cc][r1,r2]*delta_factor(bg[cc][r1])*delta_factor(bg[cc][r2]),1:length(r_grid)),1:length(r_grid)),1:length(bg)),1:length(mList))
    #twoPtFct=map(m->map(cc->map(r1->map(r2->twoPtFct_entropy[m][cc][r1][r2],1:length(r_grid)),1:length(r_grid)),1:length(bg)),1:length(mList))
    nGrid=max(length(r_grid),n_ext_Grid)
    finalCorr=zeros(length(bg),2,nFields,nFields,length(mList),nGrid,nGrid)
    for cc in 1:length(bg)
        for m in 1:length(mList)
            for r1 in 1:length(r_grid)
                for r2 in 1:length(r_grid)
                    finalCorr[cc,1,1,1,m,r1,r2]=real(twoPtFct[m][cc][r1][r2])
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
    twoPtFct_entropy=generate_2ptfct_all(norm,batches, CoM,mList;r_grid=r_grid,step=step)
    twoPtFct=map(m->map(cc->map(r1->map(r2->twoPtFct_entropy[m][cc][r1,r2]*delta_factor(bg[cc][r1])*delta_factor(bg[cc][r2]),1:length(r_grid)),1:length(r_grid)),1:length(bg)),1:length(mList))
    bgString,twoptString=construct_trento_names(participants;extensionString="dat",mMode="2")
    writedlm(bgString,bg)
    for i in 1:length(mList)
        bgString,twoptString=construct_trento_names(participants;extensionString="dat",mMode=mList[i])
        writedlm(twoptString,real.(twoPtFct[i]))
    end

end


function generate_bg_two_pt_fct_save(f,delta_factor,norm,Projectile1,Projectile2,w,k,p,sqrtS,bins,mList;minBiasEvents=1000,r_grid=0:1:10,step=2pi/20,Threaded=true,n_ext_Grid=0,nFields=10,
    extensionString="dat",path="./",override_files=false,selected_bins=nothing)
   
    if (length(bins)+1)*100>minBiasEvents
       error("Not enough events for number of bins, increase minBiasEvents")
    end
    participants=Participants(Projectile1,Projectile2,w,sqrtS,k,p)
    nGrid=max(length(r_grid),n_ext_Grid)
    finalCorr=zeros(length(bins),2,nFields,nFields,length(mList),nGrid,nGrid)
    bg=zeros(length(bins),nGrid)    
    
    evt_gen_counter = 0 

    # Determine which bins to process
    bin_indices = eachindex(bins)
    if selected_bins !== nothing
        # selected_bins can be a single value or a vector of bin ranges, e.g. [2] or [2,3]
        bin_indices = selected_bins
    end

    for cc in bin_indices
        if cc==1
            lb = 0
            rb = bins[cc]
        else
            lb = bins[cc-1]
            rb = bins[cc]
        end
        bg_file = check_for_config(participants;extensionString=extensionString,path=path,cc=string(lb) * "-" * string(rb))[1]
        bgString=construct_trento_names(participants;extensionString=extensionString,cc=string(lb) * "-" * string(rb),path=path)[1]   
        if bg_file ==true && override_files==false
            bg[cc,:]=collect(readdlm(bgString))
        elseif (bg_file == false || override_files==true) && evt_gen_counter == 0
            evt_gen_counter += 1    
            events=rand(threaded(participants),minBiasEvents)
            batches, CoM=centralities_selection_CoM(events,bins;Threaded=Threaded)
        end
        for i in eachindex(mList)
            bg_file, twpt_file = check_for_config(participants;extensionString=extensionString,mMode=mList[i],path=path,cc=string(lb) * "-" * string(rb))
            bgString,twoptString=construct_trento_names(participants;extensionString=extensionString,mMode=mList[i],cc=string(lb) * "-" * string(rb),path=path)  
            if twpt_file == true && override_files==false
                finalCorr[cc,1,1,1,i,:,:]=collect(readdlm(twoptString))
            elseif (twpt_file == false || override_files==true) && evt_gen_counter == 0
                evt_gen_counter += 1    
                events=rand(threaded(participants),minBiasEvents)
                batches, CoM=centralities_selection_CoM(events,bins;Threaded=Threaded)
            end
        end

            
        if bg_file == false || override_files==true
            bg_small_grid=generate_background(f,norm,batches,CoM,cc,r_grid=r_grid,step=step)
            for r1 in 1:length(r_grid)
                bg[cc,r1]=bg_small_grid[r1]
            end
            writedlm(bgString,bg)
        end
               
        for i in eachindex(mList)
            bg_file, twpt_file = check_for_config(participants;extensionString=extensionString,mMode=mList[i],path=path,cc=string(lb) * "-" * string(rb))
            bgString,twoptString=construct_trento_names(participants;extensionString=extensionString,mMode=mList[i],cc=string(lb) * "-" * string(rb),path=path)   
            if twpt_file == false || override_files==true
                twoPtFct_entropy=generate_2ptfct(norm,batches,CoM,mList[i],cc;r_grid=r_grid,step=step)
                twoPtFct=map(r1->map(r2->twoPtFct_entropy[r1,r2]*delta_factor(bg[r1])*delta_factor(bg[r2]),1:length(r_grid)),1:length(r_grid))
                for r1 in 1:length(r_grid)
                    for r2 in 1:length(r_grid)
                        finalCorr[cc,1,1,1,i,r1,r2]=real(twoPtFct[r1][r2])
                    end
                end
                writedlm(twoptString,finalCorr[cc,1,1,1,i,:,:])
            end
        end
    end  
    return bg,finalCorr;
end

function generate_bg_two_pt_fct_save_faster(f,delta_factor,norm::Float64,Projectile1::NucleiWoodSaxon3D,Projectile2::NucleiWoodSaxon3D,w::Float64,k::Float64,p::Float64,sqrtS::Float64,bins::Vector{Int64},mList::Vector{Int64};minBiasEvents=1000,r_grid=0:1.:10,step=2pi/20,Threaded=false,n_ext_Grid=0,nFields=10,
    extensionString="dat",path="./",override_files=false)
#=
    # Basic validation
    if (length(bins) + 1) * 100 > minBiasEvents
        error("Not enough events for number of bins, increase minBiasEvents")
    end

    # Initialize participant system and arrays
    participants = Participants(Projectile1, Projectile2, w, sqrtS, k, p)
    nGrid = max(length(r_grid), n_ext_Grid)
    finalCorr = zeros(eltype(r_grid),length(bins), 2, nFields, nFields, length(mList), nGrid, nGrid)
    bg = zeros(eltype(r_grid),length(bins), nGrid)

    # Loop over all bins
    for cc in eachindex(bins)
        # Define lower and upper bin edges
        lb = cc == 1 ? 0 : bins[cc - 1]
        rb = bins[cc]

        # Background file path
        bgString = construct_trento_names(
            participants;
            extensionString = extensionString,
            cc = string(lb) * "-" * string(rb),
            path = path
        )[1]

        if override_files
            # =====================================
            # REGENERATE AND OVERWRITE ALL FILES
            # =====================================
            events = rand(threaded(participants), minBiasEvents)
            batches, CoM = centralities_selection_CoM(events, bins; Threaded = Threaded)

            # Generate background
            bg_small_grid=generate_background(f,norm,batches,CoM,cc,r_grid=r_grid,step=step)
            for r1 in 1:length(r_grid)
                bg[cc,r1]=bg_small_grid[r1]
            end
            writedlm(bgString,bg[cc,:])

            # Generate two-point correlation functions
          @inbounds  for i in eachindex(mList)
                _, twoptString = construct_trento_names(
                    participants;
                    extensionString = extensionString,
                    mMode = mList[i],
                    cc = string(lb) * "-" * string(rb),
                    path = path
                )

                 twoPtFct_entropy=generate_2ptfct(norm,batches,CoM,mList[i],cc;r_grid=r_grid,step=step)
                twoPtFct=map(r1->map(r2->twoPtFct_entropy[r1,r2]*delta_factor(bg_small_grid[r1])*delta_factor(bg_small_grid[r2]),1:length(r_grid)),1:length(r_grid))
                @inbounds for r1 in 1:length(r_grid)
                    @inbounds for r2 in 1:length(r_grid)
                        finalCorr[cc,1,1,1,i,r1,r2]=real(twoPtFct[r1][r2])
                    end
                end
                writedlm(twoptString,finalCorr[cc,1,1,1,i,:,:])
                
            end

        else
            # =====================================
            # READ EXISTING FILES ONLY
            # =====================================
            bg[cc, :] = collect(readdlm(bgString,eltype(r_grid)))

            for i in eachindex(mList)
                _, twoptString = construct_trento_names(
                    participants;
                    extensionString = extensionString,
                    mMode = mList[i],
                    cc = string(lb) * "-" * string(rb),
                    path = path
                )
                finalCorr[cc, 1, 1, 1, i, :, :] = collect(readdlm(twoptString,eltype(r_grid)))
            end
        end
    end

    return bg, finalCorr
    =#
end


function change_norm(bg,finalCorr,new_norm,old_norm,eos)
    rescaling_factor=new_norm/old_norm
    entropy(T)=pressure_derivative(T,Val(1),eos) #entropy as function of temperature
    entropyToTemp(T)=InverseFunction(entropy)(T) #inverse, i.e. T(s)
    dSdT(T)=pressure_derivative(T,Val(2),eos) #function to convert perturbations as function of bg temp, i.e. dT/ds(T_0)
    dSdTinverse(T) = 1/dSdT(T)
    rescaled_bg=entropyToTemp.(entropy.(bg) .* rescaling_factor)
    rescaled_finalCorr=finalCorr .* (rescaling_factor^2)
    return rescaled_bg, rescaled_finalCorr
end