using MonteCarloGlauber
using Plots

# Generate fresh sample of events for histogram inspection
n1 = Lead()
n2 = Lead()

w = 0.5
s_NN = 2760
k = 1
p = 0.0

participants = Participants(n1, n2, w, s_NN, k, p, (0, 0))
events = rand(participants, 1)

evt = events[1]
grid = -10:0.1:10
profile = [evt(x, y) for x in grid, y in grid]

p0 = heatmap(
    grid, grid,
    profile,
    colorbar = :none,
    aspect_ratio = :equal,
    size = (600, 600),
    dpi = 200, axis = false,
    guidefontsize = 14, tickfontsize = 12, title = "Pb - Pb", titlefontsize = 20
)


n1 = Oxygen()
n2 = Oxygen()
participants = Participants(n1, n1, w, s_NN, k, p, (0, 0))
events = rand(participants, 1)

evt = events[1]
profile = [evt(x, y) for x in grid, y in grid]

p1 = heatmap(
    grid, grid,
    profile,
    colorbar = :none,
    aspect_ratio = :equal,
    size = (600, 600),
    dpi = 200, axis = false,
    guidefontsize = 14, tickfontsize = 12, title = "O - O", titlefontsize = 20
)


n1 = Neon()
n2 = Neon()
participants = Participants(n1, n1, w, s_NN, k, p, (0, 0))
events = rand(participants, 1)

evt = events[1]
profile = [evt(x, y) for x in grid, y in grid]

pm1 = heatmap(
    grid, grid,
    profile,
    colorbar = :none,
    aspect_ratio = :equal,
    size = (600, 600),
    dpi = 200, axis = false,
    guidefontsize = 14, tickfontsize = 12, title = "Ne - Ne", titlefontsize = 20
)

plot(
    p0, p1, pm1,
    layout = (1, 3),
    size = (900, 300),
    dpi = 200, grid = false

)
