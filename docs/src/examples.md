# Examples

## Centrality Selection
`centralities_selection_events` partitions events by multiplicity (most central first). The `bins` vector defines the upper edges of the centrality classes in percent.

```julia
n1 = Oxygen()
n2 = Oxygen()

participants = Participants(n1, n2, 0.5, 2760, 1.0, 0.0)
events = rand(threaded(participants), 200_000)

bins = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
batches = MonteCarloGlauber.centralities_selection_events(events, bins)

central_0_10 = batches[1]
peripheral_50_60 = batches[6]
```

## Background and Two-Point Correlator
This is a high-level workflow that builds the background profile and correlators. It can be expensive for large `minBiasEvents`.

```julia
fmGeV = 5.0
entropy(T) = 47.5 * 4 * T^3 * pi^2 / 90 * fmGeV^3
entropyToTemp = InverseFunction(entropy)
dSdT(T) = 1 / (3 * 47.5 * 4 * T^2 * pi^2 / 90 * fmGeV^2)

bg, twopt = generate_bg_two_pt_fct(
    entropyToTemp,
    dSdT,
    1.0,
    Lead(),
    Lead(),
    0.5,
    1.0,
    0.0,
    2760,
    [10, 20],
    [2];
    minBiasEvents = 10_000,
    r_grid = 0:1:10,
    nFields = 10,
)
```
