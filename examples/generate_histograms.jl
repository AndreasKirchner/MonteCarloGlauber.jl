using MonteCarloGlauber
using Plots

# Generate fresh sample of events for histogram inspection
n1 = Lead()
n2 = Lead()

w = 0.5
s_NN = 2760
k = 1
p = 0.0

participants = Participants(n1, n2, w, s_NN, k, p)
events = rand(threaded(participants), 10_000)

# Extract event properties
mult = multiplicity.(events)
coll = getfield.(events, :n_coll)
b = impactParameter.(events)

# Impact Parameter Histogram
p_b = histogram(
    b,
    bins = 50,
    xlabel = "Impact Parameter b [fm]",
    ylabel = "Number of Events",
    title = "Impact Parameter Distribution",
    legend = false,
    grid = false,
    fillalpha = 0.7,
    linewidth = 0,
    size = (600, 400),
    guidefontsize = 12,
    tickfontsize = 10,
    dpi = 200
)

# Number of Binary Collisions Histogram
p_coll = histogram(
    coll,
    bins = 30,
    xlabel = "Number of Collisions (N_coll)",
    ylabel = "Number of Events",
    title = "Binary Collision Distribution",
    legend = false,
    grid = false,
    fillalpha = 0.7,
    linewidth = 0,
    yscale = :log10,
    size = (600, 400),
    guidefontsize = 12,
    tickfontsize = 10,
    dpi = 200
)

# Multiplicity Histogram (log scale)
p_mult = histogram(
    mult,
    bins = 50,
    xlabel = "Multiplicity",
    ylabel = "Number of Events",
    title = "Event Multiplicity Distribution",
    legend = false,
    grid = false,
    fillalpha = 0.7,
    linewidth = 0,
    yscale = :log10,
    size = (600, 400),
    guidefontsize = 12,
    tickfontsize = 10,
    dpi = 200
)

# Combined plot
plot(
    p_b, p_coll, p_mult, layout = (1, 3), size = (1800, 400), dpi = 200, left_margin = 10Plots.mm,
    right_margin = 5Plots.mm,
    top_margin = 0.0Plots.mm,
    bottom_margin = 10Plots.mm
)
