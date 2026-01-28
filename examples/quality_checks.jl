using Statistics
using MonteCarloGlauber
using Plots

n1 = Lead()
n2 = Lead()
w = 0.5
s_NN = 2760
k = 1


p = 0.0
participants = Participants(n1, n1, w, s_NN, k, p)
evt = rand(participants, 10000)

p = 1.0
participants = Participants(n1, n1, w, s_NN, k, p)
evt_p1 = rand(participants, 1000)

p = -1.0
participants = Participants(n1, n1, w, s_NN, k, p)
evt_pm1 = rand(participants, 1000)


ecc = MonteCarloGlauber.eccentricities_gl.(evt)
ecc_p1 = MonteCarloGlauber.eccentricities_gl.(evt_p1)
ecc_pm1 = MonteCarloGlauber.eccentricities_gl.(evt_pm1)

eps2_vec = MonteCarloGlauber.eps2.(ecc)
eps2_vec_p1 = MonteCarloGlauber.eps2.(ecc_p1)
eps2_vec_pm1 = MonteCarloGlauber.eps2.(ecc_pm1)


eps3_vec = MonteCarloGlauber.eps3.(ecc)
eps3_vec_p1 = MonteCarloGlauber.eps3.(ecc_p1)
eps3_vec_pm1 = MonteCarloGlauber.eps3.(ecc_pm1)


function select_cc_events_glauber(events, cc_fraction)
    glauber_vec = map(events) do x
        x.multiplicity
    end
    sorted_indices = sortperm(glauber_vec, rev = true)
    Nev = length(glauber_vec)
    cc_ev_num = div(cc_fraction * Nev, 100)
    selected_data = [events[sorted_indices[(1 + i * cc_ev_num):(cc_ev_num + i * cc_ev_num)]] for i in 0:(cc_fraction - 1)]
    return selected_data, sorted_indices
end

data, sorted_indices = select_cc_events_glauber(evt, 10)
data_p1, sorted_indices_p1 = select_cc_events_glauber(evt_p1, 10)
data_pm1, sorted_indices_pm1 = select_cc_events_glauber(evt_pm1, 10)


begin
    len_batch = length(data[1])
    eps2_avg = [mean(getindex.(eps2_vec[sorted_indices[(1 + i * len_batch):(len_batch + i * len_batch)]], 1)) for i in 0:9]
    scatter(eps2_avg, label = "p = 0")
    len_batch = length(data_p1[1])
    eps2_p1_avg = [mean(getindex.(eps2_vec_p1[sorted_indices_p1[(1 + i * len_batch):(len_batch + i * len_batch)]], 1)) for i in 0:9]
    scatter!(eps2_p1_avg, label = "p = 1")
    eps2_pm1_avg = [mean(getindex.(eps2_vec_pm1[sorted_indices_pm1[(1 + i * len_batch):(len_batch + i * len_batch)]], 1)) for i in 0:9]
    scatter!(eps2_pm1_avg, label = "p = -1", xlabel = "centrality", ylabel = "ϵ2")
end

begin
    len_batch = length(data[1])
    eps2_avg = [mean(getindex.(eps3_vec[sorted_indices[(1 + i * len_batch):(len_batch + i * len_batch)]], 1)) for i in 0:9]
    scatter(eps2_avg, label = "p = 0")
    len_batch = length(data_p1[1])
    eps2_p1_avg = [mean(getindex.(eps3_vec_p1[sorted_indices_p1[(1 + i * len_batch):(len_batch + i * len_batch)]], 1)) for i in 0:9]
    scatter!(eps2_p1_avg, label = "p = 1")
    eps2_pm1_avg = [mean(getindex.(eps3_vec_pm1[sorted_indices_pm1[(1 + i * len_batch):(len_batch + i * len_batch)]], 1)) for i in 0:9]
    scatter!(eps2_pm1_avg, label = "p = -1", xlabel = "centrality", ylabel = "ϵ3")
end


mult = map(evt) do x
    x.multiplicity
end
n_coll = map(evt) do x
    x.n_coll
end
b = map(evt) do x
    impactParameter(x)
end

histogram(b, bins = 10)
histogram(mult, bins = 10, yscale = :log10, normalize = :pdf)


begin
    len_batch = length(data[1])
    ncoll_avg = [mean(n_coll[sorted_indices[(1 + i * len_batch):(len_batch + i * len_batch)]]) for i in 0:9]
    mult_avg = [mean(mult[sorted_indices[(1 + i * len_batch):(len_batch + i * len_batch)]]) for i in 0:9]
    scatter(mult_avg, yscale = :log10, label = "multiplicity")
    scatter!(ncoll_avg, yscale = :log10, label = "Ncoll", xlabel = "centrality")
end
