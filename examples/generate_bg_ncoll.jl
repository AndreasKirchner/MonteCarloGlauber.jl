using MonteCarloGlauber

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

#entropyToTemp(T) = T
#dSdT(T) = 1.
#dSdTinverse(T) = 1.

T = generate_bg(entropyToTemp, Norm, n1, n2, w, k, p, s_NN, centrality_bins)

ncoll = generate_ncoll(n1, n2, w, k, p, s_NN, centrality_bins)
