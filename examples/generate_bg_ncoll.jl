using MonteCarloGlauber
using Plots
const tau0 = 0.4
fmGeV = 5.0
n1 = Lead()
n2 = Lead()
w = 0.4
s_NN = 2760.0
k = 10.0
p = 0.0
Norm = 50.0
centrality_bins = [10, 20]

entropy(T) = 47.5 * 4 * T^3 * pi^2 / 90 * fmGeV^3
entropyToTemp = InverseFunction(entropy)
dSdT(T) = 1 / (3 * 47.5 * 4 * T^2 * pi^2 / 90 * fmGeV^2)


T = generate_bg(entropyToTemp, Norm, n1, n2, w, k, p, s_NN, centrality_bins)

ncoll = generate_ncoll(n1, n2, w, k, p, s_NN, centrality_bins)

T, ncoll = generate_bg_ncoll(entropyToTemp, Norm, n1, n2, w, k, p, s_NN, centrality_bins)

plot(transpose(T), label = "T [GeV]")
plot(transpose(ncoll), label = "ncoll [fm-2]")

T, ncoll = generate_bg_ncoll_save(entropyToTemp, Norm, n1, n2, w, k, p, s_NN, centrality_bins)
