# Quickstart

## Installation
```julia
import Pkg
Pkg.add(url = "https://github.com/AndreasKirchner/MonteCarloGlauber.jl")
```

Then load the package:
```julia
using MonteCarloGlauber
```

## Sample Events
```julia
n1 = Lead()
n2 = Lead()

w = 0.5
s_NN = 2760
k = 1
p = 0.0

participants = Participants(n1, n2, w, s_NN, k, p)
events = rand(participants, 10_000)
```

## Multithreaded Sampling
To use all threads when sampling many events:
```julia
events = rand(threaded(participants), 100_000)
```

## Basic Observables
```julia
mult = multiplicity.(events)
coll = getfield.(events, :n_coll)
b = impactParameter.(events)
```
![Histograms](assets/histograms.png)

## Evaluate an Event on a Grid
`Participant` objects are callable and return the fluctuating thickness at `(x, y)`:
```julia
evt = events[1]
profile = [evt(x, y) for x in -10:0.5:10, y in -10:0.5:10]
```
![Example plot](assets/example_plot.png)


