# Background

## Model Parameters
The Glauber model here exposes the following common parameters:
- `w`: Sub-nucleon width (fm). Controls the smearing of participant thickness.
- `k`: Fluctuation (shape) parameter for Gamma-distributed participant weights.
- `p`: Reduced thickness parameter in the generalized mean.
- `s_NN`: Center-of-mass energy in GeV; used to compute the inelastic cross section.

## Impact Parameter
When you call `Participants(n1, n2, w, s_NN, k, p)` without an explicit impact parameter range, the code generates minimum-bias events with an internal default range. To restrict the impact parameter distribution, pass a tuple `(bmin, bmax)` as the last positional argument:

```julia
participants = Participants(n1, n2, w, s_NN, k, p, (1.0, 2.0))
```

## Nucleus Choices
Use one of the built-in nucleus helpers:
- `Lead()`, `Gold()`, `Uranium()`, `Copper()`, `Xenon()` for Woods-Saxon nuclei.
- `Oxygen()` and `Neon()` for tabulated light-ion configurations.

## Event Profiles
`Participant` objects are callable: `evt(x, y)` returns the fluctuating thickness at coordinates `(x, y)`. This enables direct integration for eccentricities or sampling on grids for visualization.
