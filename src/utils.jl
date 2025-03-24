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