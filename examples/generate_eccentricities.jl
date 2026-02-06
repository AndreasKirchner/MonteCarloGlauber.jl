using Statistics
using MonteCarloGlauber
using Plots

n1 = Lead()
n2 = Lead()
w = 0.5
s_NN = 5500
k = 1

p = 0.0
participants = Participants(n1, n1, w, s_NN, k, p)
evt = rand(threaded(participants), 50_000)

bins = [10, 20, 30, 40, 50]
batches = centralities_selection_events(evt, bins)

function avg_eccentricity(batch::Vector, n)
    ecc = epsilon_n_psi_n.(batch, n)
    return mean(getindex.(ecc, 1)), mean(getindex.(ecc, 2))
end

# ============================================================================
# Plot eccentricity vs centrality
# ============================================================================

# Apply avg_eccentricity to each batch and extract magnitude (first element)
ecc2 = [avg_eccentricity(batches[i], 2) for i in eachindex(batches)]
ecc3 = [avg_eccentricity(batches[i], 3) for i in eachindex(batches)]

ecc2_magnitude = getindex.(ecc2, 1)
ecc3_magnitude = getindex.(ecc3, 1)
ecc2_angle = getindex.(ecc2, 2)
ecc3_angle = getindex.(ecc3, 2)
# Create centrality labels
centrality_labels = ["0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-100%"]
centrality_indices = 1:length(batches)

# Plot eccentricity magnitudes across centrality bins
p_ecc = plot(
    centrality_indices,
    ecc2_magnitude,
    label = "ε₂",
    marker = :circle,
    linewidth = 2,
    markersize = 8,
    grid = false,
    xlabel = "Centrality Bin",
    ylabel = "Eccentricity Magnitude",
    legend = :topleft,
    size = (700, 450),
    guidefontsize = 12,
    tickfontsize = 10,
    dpi = 150
)

plot!(
    centrality_indices,
    ecc3_magnitude,
    label = "ε₃",
    marker = :square,
    linewidth = 2,
    markersize = 8
    , left_margin = 5Plots.mm,
    right_margin = 5Plots.mm,
    top_margin = 6Plots.mm,
    bottom_margin = 6Plots.mm
)


# Customize x-axis ticks to show centrality ranges
xticks!(p_ecc, centrality_indices, centrality_labels)

# Plot eccentricity magnitudes across centrality bins
p_ecc = plot(
    centrality_indices,
    ecc2_angle,
    label = "ε₂ angle",
    marker = :circle,
    linewidth = 2,
    markersize = 8,
    grid = false,
    xlabel = "Centrality Bin",
    ylabel = "Eccentricity Angle",
    legend = :topleft,
    size = (700, 450),
    guidefontsize = 12,
    tickfontsize = 10,
    dpi = 150
)

plot!(
    centrality_indices,
    ecc3_angle,
    label = "ε₃ angle",
    marker = :square,
    linewidth = 2,
    markersize = 8
)

xticks!(p_ecc, centrality_indices, centrality_labels)
