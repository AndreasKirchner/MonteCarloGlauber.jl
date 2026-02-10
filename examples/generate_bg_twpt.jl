using MonteCarloGlauber
using Plots

n1 = Neon()
n2 = Neon()
w = 0.5
s_NN = 5500
k = 1
p = 0.0

fmGeV = 5.0
entropy(T) = 47.5 * 4 * T^3 * pi^2 / 90 * fmGeV^3
entropyToTemp = InverseFunction(entropy)
dSdT(T) = 1 / (3 * 47.5 * 4 * T^2 * pi^2 / 90 * fmGeV^2)

Norm = 100
m_list = [2, 3]
bins = [10, 20]

bg, twpt = generate_bg_twpt_fct(
    entropyToTemp, dSdT,
    Norm,
    n1, n2, w, k, p, s_NN,
    bins, m_list;
    minBiasEvents = 1000, r_grid = 0:1.0:5, NumPhiPoints = 20, Threaded = true, Nfields = 10
)

bg, twpt = generate_bg_twpt_fct_save(
    entropyToTemp, dSdT,
    Norm,
    n1, n2, w, k, p, s_NN,
    bins, m_list;
    minBiasEvents = 1000,
    r_grid = 0:1.0:5,
    NumPhiPoints = 20,
    Threaded = true,
    Nfields = 10,
    override_files = false
)


heatmap(twpt[1, 1, 1, 1, 1, :, :])
plot(bg[1, :])
