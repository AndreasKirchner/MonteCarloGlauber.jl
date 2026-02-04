# MonteCarloGlauber.jl

MonteCarloGlauber.jl generates initial conditions for heavy-ion collisions using a Monte Carlo Glauber model with sub-nucleonic fluctuations. It provides 3D Woods-Saxon sampling, participant and binary-collision counting, and utilities for centrality selection and two-point correlators.

## Key Features
- Deformed 3D Woods-Saxon nuclei (`NucleiWoodSaxon3D`).
- Integrated Woods-Saxon sampling with interpolation (`IntegratedWoodSaxon`, `IntegratedWoodSaxonInterp`).
- Event-by-event participants and multiplicities (`Participants`, `Participant`).
- Light-ion configurations via artifacts (`Oxigen`, `Neon`).
- Background and two-point correlator helpers (`generate_bg_two_pt_fct`, `generate_bg_twpt_fct`).

## Quickstart
```julia
using MonteCarloGlauber

n1 = Lead()
n2 = Lead()

w = 0.5
s_NN = 2760
k = 1
p = 0.0

participants = Participants(n1, n2, w, s_NN, k, p)
events = rand(participants, 1000)

mult = multiplicity.(events)
b = impactParameter.(events)
```

See the Quickstart and Examples pages for more complete workflows.
