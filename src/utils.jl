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
function centralities_selection(events::Vector{T} ;threaded=true) where {T<:Participant}
    
    if threaded 
        multiplicity=tmap(events)   do x 
            multiplicity(x)
        end
    else 
        multiplicity=map(events)   do x 
            multiplicity(x)
        end
    end 

    event_perm=sortperm(multiplicity,rev=true)
    multiplicity_sorted=multiplicity[event_perm]
 
    min_bias=length(events)
    n_event_per_bin=min_bias ÷ 100

    centrality_borders=map(x->x[1],multiplicity_sorted[1:n_event_per_bin:min_bias])
    return centrality_borders
end

