# MonteCarloGlauber

[![Build Status](https://github.com/AndreasKirchner/MonteCarloGlauber.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/AndreasKirchner/MonteCarloGlauber.jl/actions/workflows/CI.yml?query=branch%3Amaster)

To install the package and its dependencies one can use 
```julia
import Pkg
Pkg.add(url = "https://github.com/fafafrens/Fluidum.jl](https://github.com/AndreasKirchner/MonteCarloGlauber.jl")
```
After this one can use the package simply by using 
```julia
using MonteCarloGlauber
```

## Basic Expample 
```julia
using Statistics
using MonteCarloGlauber
using Plots

n1= Lead()
n2= Uranium()
```
This defines two nuclei, Lead and Uranium. To sample for it just call rand 
```julia
rand(n1)
rand(n2)
```
will sample a matrix $(A,2)$ with $A$. 
To sample more call
```julia
rand(n1,100)
```
now you will have 100 event. 

rand(threarded(n1),100)

w= 1
s_NN=5000
k=1
p=0.
b=(1,2)

participants=Participants(n1,n2,2,s_NN,k,p,b)

event=rand(participants,1000)

b_event=map(event) do x
    impactParameter(x) 
end 

ncoll_event=map(event) do x
    x.n_coll
end

profile=map(event)   do x 
    map(Iterators.product(-10:10,-10:10)) do y
        x(y...)
    end
end 

heatmap(mean(profile))


Plots.histogram(b_event,nbins=100)
Plots.histogram(ncoll_event,nbins=100)


participants=Participants(n1,n2,w,s_NN,k,p)

event=rand(participants,1000)

b_event=map(event) do x
    impactParameter(x) 
end 

ncoll_event=map(event) do x
    x.n_coll
end 

using Plots
Plots.histogram(b_event,nbins=100)
Plots.histogram(ncoll_event,nbins=100)

