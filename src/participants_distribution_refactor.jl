


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

@inline function Tp(x,y,w)
    w2=2*w^2
    1/(pi*w2)*exp(-(x^2+y^2)/w2)
end

@inline function pmeanpos(a,b,p)
    return (0.5*(a^p+b^p))^(1/p)
end

@inline function pmeanzer(a,b)
    return sqrt(a*b)
end


@inline function pmean(a,b,p)
    
    if p==0
        return pmeanzer(a,b)
    end 
    return pmeanpos(a,b,p)
end


function cross_section_from_energy(sqrtS) #returns the cross section for sqrtS in GeV
    return (0.073491826*log(sqrtS)-.19313457)*log(sqrtS)+3.123737545
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


struct Participants{A,B,C,D,E,F,G,L} <:Sampleable{Univariate,Participant}
    nucl1::A
    nucl2::B
    sub_nucleon_width::C
    inpact_parameter_magitude::D
    inpact_parameter_angle::G
    sigma_gg::E
    shape_parameter::F
    total_cross_section::Float64
    p::L
end 

"""
    Participants(n1,n2,w,s_NN,k,p[,b])

Constructor for the Participants distribution. It takes the following arguments:
    n1: Nucleus 1
    n2: Nucleus 2
    w: Subnucleon width
    s_NN: Energy in GeV
    k: Shape parameter
    p: Norm parameter
    b: Impact parameter distribution (optional) 
"""
function Participants

end 

function Participants(n1,n2,w,s_NN,k,p,b)

    sigma_NN=cross_section_from_energy(s_NN)
    f(sigmagg,p)=totalcross_section(w,sigmagg,sigma_NN) 
    u0 = 1.
    prob = NonlinearProblem{false}(f, u0)
    sol=solve(prob,SimpleNewtonRaphson())

    sigg=sol.u
    
    if length(b) ==1 
    return Participants(n1,n2,w,TriangularDist(b,b,b),Uniform(0,2pi),sigg,k,sigma_NN,p)
    end 

    return Participants(n1,n2,w,truncated(TriangularDist(0,b[2],b[2]),b[1],b[2]),Uniform(0,2pi),sigg,k,Float64(sigma_NN),p)
end 


dimension(s::Participants{A,B,C,D,E,F,G,L}) where {A,B,C,D,E,F,G,L}  = 1


function Participants(n1,n2,w,s_NN,k,p)
    
    R1=n1.R    
    R2=n2.R   

    b=(0,R1+R2+6*w)
    Participants(n1,n2,w,s_NN,k,p,b)
end

function Base.eltype(::Participants{A,B,C,D,E,F,G,L}) where {A,B,C,D,E,F,G,L} 
    T=promote_type(C,E,F,L) 
    Participant{T,T,T,T,C,F,L}
end 



function Base.copy(s::Participants{A,B,C,D,E,F,G,L}) where {A,B,C,D,E,F,G,L} 
b=unique(Distributions.params(s.inpact_parameter_magitude))

if length(b)==1
    
    return Participants(copy(s.nucl1),copy(s.nucl2),s.sub_nucleon_width,TriangularDist(first(b),first(b),first(b)),Uniform(0,2pi),s.sigma_gg,s.shape_parameter,s.total_cross_section,s.p)

end

if length(b)==2

    return Participants(copy(s.nucl1),copy(s.nucl2),s.sub_nucleon_width,truncated(TriangularDist(0,b[2],b[2]),b[1],b[2]),Uniform(0,2pi),s.sigma_gg,s.shape_parameter,s.total_cross_section,s.p)

end

end 




@inline @fastmath function totalcross_section(w,sigmaGG,sigmaNN) #double check

    #(4*pi*w^2* (Base.MathConstants.γ - expinti(-(sigmagg/(4pi* w^2))) + log(sigmagg/(4pi* w^2))))

    C=6 #truncation for integral
    return C^2/4 + expinti(-exp(-C^2/4)*sigmaGG/(4*pi*w^2)) - expinti(-sigmaGG/(4*pi*w^2)) -sigmaNN/(4*pi*w^2)

end 


@inline @fastmath function binary_impact_parameter_probability(b,nucleos::Participants{A,B,C,D,E,F,G,L}) where {A,B,C,D,E,F,G,L}
    w=nucleos.sub_nucleon_width

    Tnn=1/(4*pi*w^2)*exp(-b^2/(4*w^2))
    return 1-exp(-nucleos.sigma_gg*Tnn)

end 


#3d case 


function Distributions.rand(rng::AbstractRNG, nucleos::Participants{NUCL1, NUCL2, C, D, E, F,G, L}) where {NUCL1, NUCL2, C, D, E, F,G, L}

    
    R1=nucleos.nucl1.R
    R2=nucleos.nucl2.R 
    n1=rand(rng,nucleos.nucl1)
    n2=rand(rng,nucleos.nucl2)
    θ_b=rand(rng,nucleos.inpact_parameter_angle)
    s_th,c_th=sincos(θ_b)
    re1=Vector{SVector{2,eltype(nucleos.nucl1)}}()
    re2=Vector{SVector{2,eltype(nucleos.nucl2)}}()

    while true 
    b=rand(rng,nucleos.inpact_parameter_magitude)
    b_vec= SVector{2}(b*c_th/2,b*s_th/2)
    ncoll=0
    @inbounds for nucl1 in axes(n1,1)
        pos1=SVector{2}(n1[nucl1,1],n1[nucl1,2])
        pos1rotshift= pos1 -b_vec
        x_1=pos1rotshift[1]
        y_1=pos1rotshift[2]
        @inbounds for nucl2 in axes(n2,1)    
            pos2=SVector{2}(n1[nucl2,1],n2[nucl2,2])   
            pos2rotshift= pos2 +b_vec
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

    if ncoll>0 
        r1=unique(re1)
        r2=unique(re2)
        k=nucleos.shape_parameter
        distribution=Gamma(k,1/k)
        shape_1=rand(rng,distribution,length(r1))
        shape_2=rand(rng,distribution,length(r2))
        return Participant(r1,r2,shape_1,shape_2,ncoll,nucleos.sub_nucleon_width,nucleos.shape_parameter,nucleos.p,R1,R2,b)
    end 
    end 
end



function Distributions._rand!(rng::AbstractRNG, nucleos::Participants{NUCL1, NUCL2, C, D, E, F, L}, A::AbstractVector{T}) where {NUCL1, NUCL2, C, D, E, F, L, T}

    for i in eachindex(A)
        A[i]= rand(rng, nucleos)
    end 

    return A

end 
