# MonteCarloGlauber

[![Build Status](https://github.com/AndreasKirchner/MonteCarloGlauber.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/AndreasKirchner/MonteCarloGlauber.jl/actions/workflows/CI.yml?query=branch%3Amain)

The MonteCarloGlauber.jl package is a tool to simulate initial conditions for heavy-ion collisions based on the Glauber model.

## Installation

To install the package and its dependencies, run:
```julia
import Pkg
Pkg.add(url = "https://github.com/AndreasKirchner/MonteCarloGlauber.jl")
```
Then load the package with:
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
This defines the two nuclei and the collision parameters. Asymmetric collisions, e.g. Lead and Uranium, are also possible.
To sample events, we first create a struct:
```julia
participants = Participants(n1,n2,w,s_NN,k,p,b)
```
If `b` is specified, the impact parameter range is used. If `b` is not specified, minimum-bias events are generated. The minimum-bias events are then sampled via:
```julia
Nev = 1000
participants = threaded(Participants(n1,n2,w,s_NN,k,p))
evt = rand(participants,Nev)
```
You can extract basic collision properties such as total entropy (a.u.), number of binary collisions, number of participants, and impact parameter with:
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

We can also plot the events by evaluating them on a grid:
```julia
profile=map(evt)   do x 
    map(Iterators.product(-10:0.5:10,-10:0.5:10)) do y
        x(y...)
    end
end
heatmap(-10:0.5:10,-10:0.5:10,profile[1])
```
It is also possible to generate the background configuration and two-point correlation function used in the FluiduM code via:
```julia
fmGeV=5 # conversion factor for units
entropy(T)=47.5*4*T^3*pi^2/90*fmGeV^3 #entropy as function of temperature
entropyToTemp(T)=InverseFunction(entropy)(T) #inverse, i.e. T(s)
dSdT(T)=1/(3*47.5*4*T^2*pi^2/90*fmGeV^2) #function to convert perturbations as function of background temp, i.e. dT/ds(T_0)
bg,twpt=generate_bg_two_pt_fct(entropyToTemp,dSdT,1,Lead(),Lead(),w,k,p,s_NN,[10,20],[2];minBiasEvents=1000,r_grid=0:1:10,nFields=10,n_ext_Grid=100)

```
We can also plot the background:
```julia
plot(0:1:10,bg[1])
plot!(0:1:10,bg[2])
```
and the two-point function:
```julia
heatmap(twpt[1,1,1,1,1,:,:])
```
The usage of light ions is also possible with only a small modification in the projectile definition.
This package ships light-ion configurations as an artifact, so you can use the built-in helpers:
```julia
n1 = Neon()
n2 = Oxigen()
```
If you want to load a specific file directly, you can still do:
```julia
n1 = TabulatedEvent(joinpath(MonteCarloGlauber.root_light_ion, "NLEFT_dmin_0.5fm_negativeweights_Ne.h5"))
```
where `root_light_ion` points to the installed artifact directory.

## Centrality selection
We create two nuclei:
```julia
n1 = Oxigen()
n2 = Oxigen()
```
We set some parameters:
```julia
w = 0.5          # nucleon width
s_NN = 2760      # Center-of-mass energy (GeV)
k = 1            # fluctuation parameter
p = 0.0             # reduced thickness parameter
```
We generate 10^6 events in a multithreaded way:
```julia
Nev = 1_000_000
participants = threaded(Participants(n1, n2, w, s_NN, k, p))
evt = rand(participants, Nev)
```
We compute some observables as before:
```julia
mult = multiplicity.(evt)
n_coll = map(evt) do x
    x.n_coll
end
Npart = MonteCarloGlauber.n_part.(evt)
b = impactParameter.(evt)
```
This function selects events for a given percentile:
```julia
function centrality_indices_by_rank(mult, bins)
    n = length(mult)
    order = sortperm(mult, rev = true)  # highest mult = most central
    out = Vector{Vector{Int}}(undef, length(bins) - 1)
    for (i, (lo, hi)) in enumerate(zip(bins[1:(end - 1)], bins[2:end]))
        lo_i = floor(Int, n * lo / 100) + 1
        hi_i = floor(Int, n * hi / 100)
        out[i] = order[lo_i:hi_i]
    end
    return out
end
```
Select the bin edges and plot the first centrality class:
```julia
bins = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
idx_bins = centrality_indices_by_rank(mult, bins)
histogram(multiplicity.(evt[idx_bins[1]]), bins = 100)
```
