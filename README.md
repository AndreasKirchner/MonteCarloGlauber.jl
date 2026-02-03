# MonteCarloGlauber

[![Build Status](https://github.com/AndreasKirchner/MonteCarloGlauber.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/AndreasKirchner/MonteCarloGlauber.jl/actions/workflows/CI.yml?query=branch%3Amain)

The MonteCarloGlauber.jl package is a tool to simulate initial conditions for heavy-ion collisions based on the Glauber model.

##Installation 

To install the package and its dependencies one can use 
```julia
import Pkg
Pkg.add(url = "https://github.com/AndreasKirchner/MonteCarloGlauber.jl")
```
After this one can use the package simply by using 
```julia
using MonteCarloGlauber
```

## Basic Usage 
```julia
using Statistics
using MonteCarloGlauber
using Plots

# Create lead nuclei
n1 = Lead()
n2 = Lead()
# Collision parameters
w = 0.5          # nucleon width
s_NN = 2760      # Center-of-mass energy (GeV)
k = 1            # fluctuation parameter
p = 0.             # reduced thickness parameter
b = (1,2)          # impact parameter (optional)
```
This defines the two nuclei and the collision parameters. Asymmetric collision, e.g. Lead and Uranium are also possible.
To sample the event, we first create a struct 
```julia
participants = Participants(n1,n2,w,s_NN,k,p,b)
```
, in case the impact parameter range is specified. If b is not specified, one creates minimum bias events. The minimum bias events are then sampled via
```julia
Nev = 1000
participants = threaded(Participants(n1,n2,w,s_NN,k,p))
evt = rand(participants,Nev)
```
One can extract basic collision properties such as total entropy (in a.u.), number of binary collisions, number of participants and impact parameter with
```julia
mult = map(evt) do x
    x.multiplicity
end

n_coll = map(evt) do x
    x.n_coll
end

Npart = map(evt) do x
    MonteCarloGlauber.n_part(x)
end

b = map(evt) do x
    impactParameter(x)
end

# Plot distributions
histogram(b, bins = 10)
histogram(mult, bins = 10, yscale = :log10)
histogram(n_coll, bins = 10)

```

With the events sampled we can now compute things such as the impact parameter histogram
```julia
b_event=map(event) do x
    impactParameter(x) 
end 
histogram(b_event,nbins=100)
```
or the histogram for the number of collisions
```julia
ncoll_event=map(event) do x
    x.n_coll
end
histogram(ncoll_event,nbins=100,yscale=:log10)
```
We can also plot the events by evaluating them on a grid
```julia
profile=map(event)   do x 
    map(Iterators.product(-10:0.5:10,-10:0.5:10)) do y
        x(y...)
    end
end
heatmap(-10:0.5:10,-10:0.5:10,profile[1])
```
It is also possible to generate the background configuration and two point correlation function used in the FluiduM code via
```julia
fmGeV=5 # conversion factor for units
entropy(T)=47.5*4*T^3*pi^2/90*fmGeV^3 #entropy as function of temperature
entropyToTemp(T)=InverseFunction(entropy)(T) #inverse, i.e. T(s)
dSdT(T)=1/(3*47.5*4*T^2*pi^2/90*fmGeV^2) #function to convert perturbations as function of bg temp, i.e. dT/ds(T_0)
bg,twpt=generate_bg_two_pt_fct(entropyToTemp,dSdT,1,Lead(),Lead(),w,k,p,s_NN,[10,20],[2];minBiasEvents=1000,r_grid=0:1:10,nFields=10,n_ext_Grid=100)

```
We can also plot the background
```julia
plot(0:1:10,bg[1])
plot!(0:1:10,bg[2])
```
and the two point function
```julia
heatmap(twpt[1,1,1,1,1,:,:])
```
The usage of light ions is also possible with the only needed modification in the projectile definition
```julia
n1=TabulatedEvent("NLEFT_dmin_0.5fm_negativeweights_Ne.h5")
```
with TabulatedEvent being the file of possible neon configurations in this case.
