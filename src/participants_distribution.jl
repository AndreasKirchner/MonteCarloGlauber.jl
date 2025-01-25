


struct Participant{T,S,V,M,C,D,F} <:ValueSupport
    part1::Vector{SVector{2}{T}}
    part2::Vector{SVector{2}{S}}
    shape1::Vector{V}
    shape2::Vector{M}
    n_coll::Int64
    sub_nucleon_width::C
    shape_parameter::D
    p::F
    R1::Float64
    R2::Float64
    b::Float64
end

impactParameter(x::Participant{T,S,V,M,C,D,F} ) where {T,S,V,M,C,D,F} = x.b

Base.eltype(::Participant{T,S,V,M,C,D,F} ) where {T,S,V,M,C,D,F}  = promote_type(T,S)

function Tp(x,y,w)
    w2=2*w^2
    1/(pi*w2)*exp(-(x^2+y^2)/w2)
end

function pmeanpos(a,b,p)
    return (0.5*(a^p+b^p))^(1/p)
end

function pmeanzer(a,b)
    return sqrt(a*b)
end


function pmean(a::Float64,b::Float64,p)
    if p>0.0 || p<0.0
        return pmeanpos(a,b,p)
    else 
        return pmeanzer(a,b)
    end
end

function (f::Participant{T,S,V,M,C,D,F} )(x,y) where {T,S,V,M,C,D,F}

    part1=f.part1
    part2=f.part2
    shape1=f.shape1
    shape2=f.shape2
    w= f.sub_nucleon_width
    p=f.p

    start=zero(eltype(f))
    ta=sum(zip(part1,shape1);init=start) do (pa ,ga)
        ga*Tp(x-pa[1],y-pa[2],w)
    end 
    start=zero(eltype(f))
    tb=sum(zip(part2,shape2);init=start) do (pa ,ga)
        ga*Tp(x-pa[1],y-pa[2],w)
    end 

    return pmean(ta,tb,p)#norm((ta,tb),p)

end

function fluctuating_thickness(x::Num1,y::Num2,part::Participant{T,S,V,M,C,D,F} ) where {Num1<:Real,Num2<:Real,T,S,V,M,C,D,F}

    part1=part.part1
    part2=part.part2
    shape1=part.shape1
    shape2=part.shape2
    w= part.sub_nucleon_width
    p=part.p

    start=zero(eltype(part))
    ta=sum(zip(part1,shape1);init=start) do (pa ,ga)
        ga*Tp(x-pa[1],y-pa[2],w)
    end 
    start=zero(eltype(part))
    tb=sum(zip(part2,shape2);init=start) do (pa ,ga)
        ga*Tp(x-pa[1],y-pa[2],w)
    end 

    return pmean(ta,tb,p)#norm((ta,tb),p)
end 


struct Participants{A,B,C,D,E,F,L} <:Sampleable{Univariate,Participant}
    nucl1::A
    nucl2::B
    sub_nucleon_width::C
    inpact_parameter::D
    sigma_gg::E
    shape_parameter::F
    total_cross_section::Float64
    p::L
end 

function Participants(n1,n2,w,b,sigma_NN,k,p)
    f(sigmagg,p)=totalcross_section(w,sigmagg,sigma_NN) 
    u0 = 1.
    prob = NonlinearProblem(f, u0)
    sol=solve(prob,SimpleNewtonRaphson())

    sigg=sol.u

    Participants(n1,n2,w,b,sigg,k,Float64(sigma_NN),p)
end 

function Participants(n1,n2,w,sigma_NN,k,p)
    if hasfield(typeof(n1), :R)
        R1=n1.R   
    else
        R1=mean(map(i->maximum(map(eachcol(n1[:,:,i])) do ll
            hypot(ll...)
            end),1:size(n1)[3]))
    end
    
    if hasfield(typeof(n2), :R)
        R2=n2.R   
    else
        R2=mean(map(i->maximum(map(eachcol(n2[:,:,i])) do ll
            hypot(ll...)
            end),1:size(n2)[3]))
    end
    b=(0,R1+R2+6*w)
    Participants(n1,n2,w,b,sigma_NN,k,p)
end

function Base.eltype(::Participants{A,B,C,D,E,F,L}) where {A,B,C,D,E,F,L} 
    T=promote_type(C,E,F,L) 
    Participant{T,T,T,T,C,F,L}
end 


function fluctuating_thickness!(out::V,part::T,grid_x::S,grid_y::R) where {V<:AbstractArray,T<:Participant,S<:AbstractVector,R<:AbstractVector}
    
    #part1=part.part1
    #part2=part.part2
    #shape1=part.shape1
    #shape2=part.shape2
    #w= part.sub_nucleon_width
    #p=part.p

    #function Tp(x,y)
    #    w2=2*w^2
    #    1/(pi*w2)*exp(-(x^2+y^2)/w2)
    #end

    @inbounds for i in eachindex(grid_x)
        x=grid_x[i]
        @inbounds for j in eachindex(grid_y)
            y=grid_y[j]
            
            #ta=sum(zip(part1,shape1)) do (pa ,ga)
            #    ga*Tp(x-pa[1],y-pa[2])
            #end 

            #tb=sum(zip(part2,shape2)) do (pa ,ga)
            #    ga*Tp(x-pa[1],y-pa[2])
            #end 
            
            #@show typeof(part2 )
            #tb=sum(zip(part2,shape2)) do pb ,gb
            #    gb*Tp(x-pb[1],y-pb[2])
            #end
            out[i,j]= part(x,y)#norm((ta,tb),p) #((ta^p+tb^p)/2)^(1/p)
        end 
    end

    return nothing

end

function fluctuating_thickness_threads!(out::V,part::T,grid_x::S,grid_y::R) where {V<:AbstractArray,T<:Participant,S<:AbstractVector,R<:AbstractVector}
    
    #part1=part.part1
    #part2=part.part2
    #shape1=part.shape1
    #shape2=part.shape2
    #w= part.sub_nucleon_width
    #p=part.p

    #function Tp(x,y)
    #    w2=2*w^2
    #    1/(pi*w2)*exp(-(x^2+y^2)/w2)
    #end

    cart = CartesianIndices(out)
    tmap!(out, eachindex(out)) do idx
        c = cart[idx]
        x=grid_x[c[1]]
        y=grid_x[c[2]]
        part(x,y)
    end
    

    return nothing

end



function fluctuating_thickness!(out::V,part::T,grid_x::S) where {V<:AbstractArray,T<:Participant,S<:AbstractVector}
    
    fluctuating_thickness!(out,part,grid_x,grid_x) 
end

function fluctuating_thickness_threads!(out::V,part::T,grid_x::S) where {V<:AbstractArray,T<:Participant,S<:AbstractVector}
    
    fluctuating_thickness_threads!(out,part,grid_x,grid_x) 
end


function fluctuating_thickness(part::T,grid_x::S,grid_y::R) where {T<:Participant,S<:AbstractVector,R<:AbstractVector}
    out=Matrix{eltype(part)}(undef,length(grid_x),length(grid_y))
    
    fluctuating_thickness!(out,part,grid_x,grid_y)
    return out 
end

function fluctuating_thickness_threads(part::T,grid_x::S,grid_y::R) where {T<:Participant,S<:AbstractVector,R<:AbstractVector}
    out=Matrix{eltype(part)}(undef,length(grid_x),length(grid_y))
    
    fluctuating_thickness_threads!(out,part,grid_x,grid_y)
    return out 
end

function fluctuating_thickness(part::T,grid_x::S) where {T<:Participant,S<:AbstractVector}
    
    fluctuating_thickness(part,grid_x,grid_x)
    
end

function fluctuating_thickness_threads(part::T,grid_x::S) where {T<:Participant,S<:AbstractVector}
    
    fluctuating_thickness_threads(part,grid_x,grid_x)
    
end

function fluctuating_thickness(part::Vector{T},grid_x::S,grid_y::R)where {T<:Participant,S<:AbstractVector,R<:AbstractVector}
    map(x->fluctuating_thickness(x,grid_x,grid_y),part)
end 

function fluctuating_thickness_threads(part::Vector{T},grid_x::S,grid_y::R)where {T<:Participant,S<:AbstractVector,R<:AbstractVector}
    tmap(x->fluctuating_thickness(x,grid_x,grid_y),part)
end 

function fluctuating_thickness(part::Vector{T},grid_x::S)where {T<:Participant,S<:AbstractVector}
    map(x->fluctuating_thickness(x,grid_x,grid_x),part)
end

function fluctuating_thickness_threads(part::Vector{T},grid_x::S)where {T<:Participant,S<:AbstractVector}
    tmap(x->fluctuating_thickness(x,grid_x,grid_x),part)
end


@inline @fastmath function totalcross_section(w,sigmaGG,sigmaNN) #double check

    #(4*pi*w^2* (Base.MathConstants.γ - expinti(-(sigmagg/(4pi* w^2))) + log(sigmagg/(4pi* w^2))))

    C=6 #truncation for integral
    return C^2/4 + expinti(-exp(-C^2/4)*sigmaGG/(4*pi*w^2)) - expinti(-sigmaGG/(4*pi*w^2)) -sigmaNN/(4*pi*w^2)

end 


@inline @fastmath function binary_impact_parameter_probability(b,nucleos::Participants{A,B,C,D,E,F,L}) where {A,B,C,D,E,F,L}
    w=nucleos.sub_nucleon_width
   # 2pi*b/nucleos.total_cross_section*(1- exp(
   #     -nucleos.sigma_gg/(4pi*w^2)*exp(-b^2/(4*w^2))
    #))
    Tnn=1/(4*pi*w^2)*exp(-b^2/(4*w^2))
    return 1-exp(-nucleos.sigma_gg*Tnn)

   #res=(1- exp(-nucleos.sigma_gg/(4pi*w^2)*exp(-b^2/(4*w^2))))*2pi*b/nucleos.total_cross_section

   #ifelse(res>1,one(res),res)

end 



function sampleb(rng,b::T) where {T<:Number}
    b
end

function sampleb(rng,b) 
    bmin=b[1]
    bmax=b[2]
    bmin+(bmax-bmin) *sqrt(rand(rng))
end

#2d case 

#3d case 

function Distributions._rand!(rng::AbstractRNG, nucleos::Participants{NUCL1, NUCL2, C, D, E, F, L}, A::AbstractVector{T}) where {NUCL1, NUCL2, C, D, E, F, L, T}

    R1=nucleos.nucl1.R
    R2=nucleos.nucl1.R
    nevent=length(A)#nucleos.n_event
    
    

    n1=sample_nevent_with_rotation_threaded(rng,nucleos.nucl1,nevent)
    n2=sample_nevent_with_rotation_threaded(rng,nucleos.nucl2,nevent)


    map!(A,eachindex(A)) do event 
    b=sampleb(rng,nucleos.inpact_parameter)
        event1=event 
        event2=event
        re1=fill(SVector{2,eltype(n1)}(zero(eltype(n1)),zero(eltype(n1))),0)
        re2=fill(SVector{2,eltype(n2)}(zero(eltype(n2)),zero(eltype(n2))),0)

        ncoll=1
        rot1=rand(rng,RotMatrix2)
        rot2=rand(rng,RotMatrix2)
         @inbounds for nucl1 in axes(n1,2)
            
            pos1=SVector{2}(n1[event1,nucl1,1],n1[event1,nucl1,2])
            pos1rot=rot1*pos1
            pos1rotshift= pos1rot +SVector{2}(-b/2,0)

            x_1=pos1rotshift[1]
            y_1=pos1rotshift[2]

            @inbounds for nucl2 in axes(n2,2)
                
                pos2=SVector{2}(n2[event2,nucl2,1],n2[event2,nucl2,2])
                pos2rot=rot2*pos2
                pos2rotshift= pos2rot +SVector{2}(b/2,0)

                x_2=pos2rotshift[1]
                y_2=pos2rotshift[2]

                impact_par=hypot(x_1-x_2,y_1-y_2)

                probability=binary_impact_parameter_probability(impact_par,nucleos)
 
                #accepted 
                if rand(rng,Bernoulli(probability))
 
                    push!(re1,SVector{2}(x_1,y_1))
                    push!(re2,SVector{2}(x_2,y_2))
             
                    ncoll+=1
                end 
         
                
            end 
        end 
        

        r1=re1#unique(re1)#this might be wrong
        r2=re2#unique(re2)
        k=nucleos.shape_parameter
        shape_1=rand(rng,Gamma(k,1/k),length(r1))
        shape_2=rand(rng,Gamma(k,1/k),length(r2))
  
        Participant(r1,r2,shape_1,shape_2,ncoll,nucleos.sub_nucleon_width,nucleos.shape_parameter,nucleos.p,R1,R2,b)
    end
    return A
end

function Distributions.rand(rng::AbstractRNG, nucleos::Participants{NUCL1, NUCL2, C, D, E, F, L}) where {NUCL1, NUCL2, C, D, E, F, L}

 #nucleos.n_event
 R1=nucleos.nucl1.R
 R2=nucleos.nucl1.R
    totsize1=nucleos.nucl1.N_nucleon
    totsize2=nucleos.nucl2.N_nucleon

    b=sampleb(rng,nucleos.inpact_parameter)



        n1_=transpose(rand(rng,nucleos.nucl1,totsize1))
        
            r = rand(rng,RotMatrix3)
            @inbounds for n in axes(n1_,1)
                pos=SVector{3}(n1_[n,1],n1_[n,2],n1_[n,3])
                rot=r*pos
                n1_[n,:].=rot
            end 

        
        n1= @views n1_[:,1:2]


        n2_=transpose(rand(rng,nucleos.nucl2,totsize2))
        
            r = rand(rng,RotMatrix3)
            @inbounds for n in axes(n2_,1)
                pos=SVector{3}(n2_[n,1],n2_[n,2],n2_[n,3])
                rot=r*pos
                n2_[n,:].=rot
            end 
        
        n2= @views n2_[:,1:2]

        re1=fill(SVector{2,eltype(n1)}(zero(eltype(n1)),zero(eltype(n1))),0)
        re2=fill(SVector{2,eltype(n2)}(zero(eltype(n2)),zero(eltype(n2))),0)

        ncoll=1
        rot1=rand(rng,RotMatrix2)
        rot2=rand(rng,RotMatrix2)
         @inbounds for nucl1 in axes(n1,1)


            pos1=SVector{2}(n1[nucl1,1],n1[nucl1,2])
            pos1rot=rot1*pos1
            pos1rotshift= pos1rot +SVector{2}(-b/2,0)

            x_1=pos1rotshift[1]
            y_1=pos1rotshift[2]
            

            @inbounds for nucl2 in axes(n2,1)
                
                pos2=SVector{2}(n1[nucl2,1],n2[nucl2,2])
                pos2rot=rot2*pos2
                pos2rotshift= pos2rot +SVector{2}(b/2,0)

                x_2=pos2rotshift[1]
                y_2=pos2rotshift[2]

                impact_par=hypot(x_1-x_2,y_1-y_2)

                probability=binary_impact_parameter_probability(impact_par,nucleos)
 
                #accepted 
                if rand(rng,Bernoulli(probability))
                
                    push!(re1,SVector{2}(x_1,y_1))
                    push!(re2,SVector{2}(x_2,y_2))
                    
                    ncoll+=1
                end 
         
                
            end 
        end 
        
        r1=re1
        r2=re2
        k=nucleos.shape_parameter
        shape_1=rand(rng,Gamma(k,1/k),length(r1))
        shape_2=rand(rng,Gamma(k,1/k),length(r2))
        return Participant(r1,r2,shape_1,shape_2,ncoll,nucleos.sub_nucleon_width,nucleos.shape_parameter,nucleos.p,R1,R2,b)
    
end
