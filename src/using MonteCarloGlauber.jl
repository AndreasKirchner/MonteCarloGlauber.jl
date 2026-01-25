using Revise
using Statistics
using MonteCarloGlauber
using StatsBase
using Plots
using StaticArrays
using OhMyThreads
using DifferentiationInterface
using IterTools
using BenchmarkTools
using DelimitedFiles
#using Cuba
#using FastGaussQuadrature
n1 = Lead()
n2 = Lead()
w = 0.5
s_NN = 2760
k = 1
p = 0.0


n1 = TabulatedEvent("src/NLEFT_dmin_0.5fm_negativeweights_Ne.h5")
using HDF5

open("src/NLEFT_dmin_0.5fm_negativeweights_Ne.h5", "r")


participants = Participants(n1, n1, w, s_NN, k, p)

e(T) = T^3pi
e(2)

energy2(0.6)

fmGeV = 5
energy2(T) = 47.5 * T^3 * pi^2 / 90 / 4 * fmGeV^3 #substitute entropy of the ideal gas of quarks (u,d,s) and gluons
entropyToEnergy(T) = InverseFunction(energy2)(T)
norm = 40
bg = MonteCarloGlauber.generate_bg(entropyToEnergy, norm, n1, n1, w, k, p, s_NN, [5]; minBiasEvents = 100, r_grid = 0:0.1:10)

plot(bg[1])

entropyToEnergy(1)

energy2(0.1)

rand(threaded(n1), 100)
participants = Participants(n1, n2, w, s_NN, k, p)
threaded(Participants(n1, n1, w, s_NN, k, p))

#INITIAL CONDITIONS
n1 = Lead()
n2 = Lead()
w = 1
s_NN = 5000
k = 1
p = 0.0
Norm = 50


participants = Participants(n2, n1, w, s_NN, k, p, 0.0)

participantsN = Participants(n1, n1, w, s_NN, k, p)

participantsN
participants

participants.nucl1
participantsN.nucl1
participants = Participants(n2, n1, w, s_NN, k, p, 0.0)
eve1 = rand(threaded(participants), 5)
eve2 = rand(threaded(participantsN), 5)

heatmap(eve1[1])


rand(participants.nucl1)
eve1.n_coll
eve2.n_coll

eve1.part1
eve1.part2
eve2.part1


event = rand(participantsN, 2)

event[1].n_coll
event[2].n_coll

event[1].part1
event[2].part2

participants = Participants(n2, n1, w, s_NN, k, p, (0.0, 1.0))
eve1 = rand(threaded(participants), 10)
profile = map(eve1) do x
    map(Iterators.product(-20:0.2:20, -20:0.2:20)) do y
        x(y...)
    end
end

heatmap(-20:0.2:20, -20:0.2:20, log.(profile[9]))
heatmap(-20:0.2:20, -20:0.2:20, profile[1])

com = center_of_mass(eve1[3])

com[2] / com[1]
com[3] / com[1]

threaded

heatmap(-10:0.1:10, -10:0.1:10, profile[6])
heatmap(-10:0.1:10, -10:0.1:10, profile[3])
heatmap(-10:0.1:10, -10:0.1:10, profile[4])


b_event = map(event) do x
    impactParameter(x)
end

ncoll_event = map(event) do x
    x.n_coll
end

profile = map(event) do x
    map(Iterators.product(-10:0.5:10, -10:0.5:10)) do y
        x(y...)
    end
end
using Plots

heatmap(mean(profile))

heatmap(profile[1])
heatmap(profile[2])
heatmap(profile[3])
heatmap(profile[4])

generate_bg_two_pt_fct()

Plots.histogram(b_event, nbins = 100)
Plots.histogram(ncoll_event, nbins = 100, normalize = true, yaxis = :log)
energy2(T) = T
bg, twpt = MonteCarloGlauber.generate_bg_two_pt_fct(energy2, energy2, 1, Lead(), Lead(), w, k, p, s_NN, [10, 20], [2]; minBiasEvents = 5000, r_grid = 0:0.2:10)

size(twpt)
heatmap(twpt[1, 1, 1, 1, 1, :, :])

using MonteCarloGlauber
bg2, twpt2 = MonteCarloGlauber.generate_bg_two_pt_fct(energy2, energy2, 1, Lead(), Lead(), w, k, p, s_NN, [10, 20], [2]; minBiasEvents = 100, r_grid = 0:1:10, nFields = 3, n_ext_Grid = 100)
heatmap(twpt2[1, 1, 1, :, :])
# the correlator has m, centra, r, r
# we want n complex, n field, nm, nr , nr
30 / 0.2
newCorr = zeros(2, 3, 1, 150, 150)
for cc in 1:2
    for m in 1:1
        for ri in 1:51
            for rj in 1:51
                newCorr[1, 1, m, ri, rj] = real(twpt[m][cc][ri][rj])
            end
        end
    end
end

heatmap(newCorr[1, 1, 1, :, :])

size(twpt)
size(twpt[1])
size(twpt[1][1])
size(twpt[1][1][1])
twpt[1][1][1][1]

using MonteCarloGlauber
using BenchmarkTools

bg2 = MonteCarloGlauber.generate_bg(energy2, 1, Lead(), Lead(), w, k, p, s_NN, [10, 20, 30, 40]; minBiasEvents = 5000, r_grid = 0:0.1:10)

plot(bg[1])
plot!(bg2[2])
plot!(bg2[3])
plot!(bg2[4])

MonteCarloGlauber.save_bg_two_pt_fct(energy2, energy2, 1, Lead(), Lead(), w, k, p, s_NN, [10, 20], [2, 3, 4]; minBiasEvents = 1000, r_grid = 0:1:4)
using DelimitedFiles
aa = readdlm("trento_two_pt_fct_m_2_Pb_Pb_w_0_5_sNN_6_206_p_1_0_k_1.dat")

aa

aa[1]

plot(bg[1], label = "0-10")
plot!(bg[2], label = "10-20")
heatmap(real.(permutedims(hcat(twpt[1][1]...))))

MonteCarloGlauber.cross_section_from_energy(s_NN)

participants = Participants(n1, n2, w, s_NN, k, p, (0.1, 6.8))
event = rand(threaded(participants), 500)
aba, acom = MonteCarloGlauber.centralities_selection_CoM(event, [10])
energy(T) = 3 * 8 * pi^2 / 15 * (3^2 - 1 + 7 / 4 * 3 * 3) * T^4 / (2 * 3 * 4)
Te4 = MonteCarloGlauber.InverseFunction(energy)

bb = MonteCarloGlauber.generate_background(Te4, 1, aba, acom; r_grid = 0:0.1:10)
plot(bb)
corr = MonteCarloGlauber.generate_2ptfct(Te4, 1, aba, acom, [2, 3, 4]; r_grid = 0:0.5:5)
heatmap(real.(corr[1][1]))
bg

writedlm("bg.dat", bg)
aaa = readdlm("bg.dat")


size(corr)
size(corr[1])
size(corr[1][1])
corr[1][1][1][1]
length(corr)

shape(corr)

heatmap(imag.(corr[1][1]))
bb[1][1]

using MonteCarloGlauber
bg, twpt = MonteCarloGlauber.generate_bg_two_pt_fct(Te, Te, 1, Lead(), Lead(), w, k, p, s_NN, [10, 20], [2, 3, 4]; minBiasEvents = 5000, r_grid = 0:0.5:5)

heatmap(real.(permutedims(hcat(twpt[1][1]...))))


b_event = map(event) do x
    impactParameter(x)
end
histogram(b_event)

multi = tmap(event) do x
    multiplicity(x, Nr = 50, Nth = 15)[1]
end

multi = tmap(event) do x
    multiplicity(x, Nr = 50, Nth = 15)
end

using DelimitedFiles
dat = readdlm("BKG276TeVAvgCC1.txt")

plot(dat[:, 2])

aba, acom = MonteCarloGlauber.centralities_selection_CoM(event, [10])

energy(T) = 3 * 8 * pi^2 / 15 * (3^2 - 1 + 7 / 4 * 3 * 3) * T^4 / (2 * 3 * 4)
Te4 = MonteCarloGlauber.InverseFunction(energy)
energy(1)
energy2(T) = T
Te = InverseFunction(energy2)

energy.([1, 2, 3])
Te4(1)
Te.([1, 2, 3])

bb = MonteCarloGlauber.generate_background(Te4, 1, aba, acom; r_grid = 0:1:10)
plot(bb, label = "new")
plot!(dat[:, 2], label = "old")
plot!(1.5 .* dat[:, 2], label = "old*1.5")

sum(dat[:, 2]) * 0.1

mean(i -> aba[1][i](0, 0), 1:length(aba[1]))

bb
bbb = acom[1]

bbb[1][2]

for (index, value) in mean(enumerate(aba[1]))
    index
end

aba[1][1](0, -1)

tl = ("a", "b", "c", "d")

sum(enumerate(aba[1])) do (i, conf)
    i
end

mean(enumerate(aba[1])) do (i, conf)
    i
end


1 + 2 + 3 + 4

enumerate(tl) do i
    i
end

bs = for (i, conf) in enumerate(tl)
    i
end

bs

grid = 0:1:4
MonteCarloGlauber.generatingfunction(aba[1], 1, 0:1:3, 0.5, acom[1], Te4, 1)
bb = MonteCarloGlauber.generate_background(Te4, 1, aba, acom; r_grid = 0:1:3)
aba[1]
acom
####this is the plot
histogram(multi, nbins = 100, normalize = true, yscale = :log10)


obj = fit(Histogram, multi; nbins = 100)

dump(obj)


profile = mapreduce(+, 1:10000) do _
    x = rand(participants)
    map(Iterators.product(-10:10, -10:10)) do y
        x(y...)
    end
end


rand(threaded(participants))
rand(threaded(participants), 10)


ncoll_event = map(event) do x
    x.n_coll
end

rand()

using Distributions

a1 = 10 * sqrt.(rand(100000))
histogram(a1, bins = 100)
aa = truncated(TriangularDist(0, 10, 10), 0, 10)
a2 = rand(aa, 100000)
histogram(a1, bins = 100)
histogram!(a2, bins = 100)


x1 = rand(Uniform(-10, 10), 100000)
x2 = rand(Uniform(-10, 10), 100000)
b = sqrt.(x1 .^ 2 .+ x2 .^ 2)
histogram!(b, normalize = true)
histogram(a1, bins = 100, normalize = true)
histogram!(a2, bins = 100, normalize = true)

dd = sqrt.(Uniform(0, 1))
dd = Uniform(0, 1)

rand(dd, 12)

n1.R

profile = map(event) do x
    map(Iterators.product(-10:10, -10:10)) do y
        x(y...)
    end
end


using BenchmarkTools

one_profile = rand(participants)
multiplicity(one_profile)
@code_warntype multiplicity(one_profile)
@benchmark multiplicity($one_profile)

using Integrals
function multiplicity_hcubature(con::T) where {T <: Participant}

    f(u, p) = SVector{3}(1, u[1], u[2]) .* con(u[1], u[2])
    domain = (SVector{2}(-Inf, -Inf), SVector{2}(Inf, Inf))
    prob = IntegralProblem{false}(f, domain)
    sol = solve(prob, GaussLegendre(); reltol = 1.0e-3, abstol = 1.0e-3)
    return sol.u


end

@benchmark  multiplicity_hcubature($one_profile)


heatmap(profile)


Plots.histogram(b_event, nbins = 100)
Plots.histogram(ncoll_event, nbins = 100)


participants = Participants(n1, n2, w, s_NN, k, p)

event = rand(participants, 1000)

b_event = map(event) do x
    impactParameter(x)
end

ncoll_event = map(event) do x
    x.n_coll
end

using Plots
Plots.histogram(b_event, nbins = 100)
Plots.histogram(ncoll_event, nbins = 100)

participants = (Participants(n1, n2, w, s_NN, k, p))

event = rand(participants, 1000)

b_event = map(event) do x
    impactParameter(x)
end

ncoll_event = map(event) do x
    x.n_coll
end

using Plots
Plots.histogram(b_event, nbins = 100)
Plots.histogram(ncoll_event, nbins = 100)


using Profile

using PProf
Profile.clear()
Profile.Allocs.clear()
@profile rand(aa, 10_000_000)
Profile.Allocs.@profile rand(aa, 100)
pprof()
rand(threarded(aa), 1000)
prof = Profile.Allocs.fetch()
PProf.Allocs.pprof(prof; from_c = false)

@benchmark rand($aa, 10)

aa = NucleiWoodSaxon3D(10, 1.5, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
aa = IntegratedWoodSaxon(10, 1, 5.0, 1.0, 0.0)
aa = IntegratedWoodSaxonInterp(10, 1, 5.0, 1.0, 0.0)
bb = threarded(aa)
cc = Threarded(aa, Threads.nthreads())
dd = Threarded(aa, 10 * Threads.nthreads())

using BenchmarkTools
@benchmark rand($aa, 10000)
@benchmark rand($bb, 10000)
@benchmark rand($cc, 10000)
@benchmark rand($dd, 10000)

Normal(0, 1)

TriangularDist(0, 1)
@code_warntype rand(Random.default_rng(), aa)
histogram(rand(TriangularDist(1, 1, 1), 100000))


rand(aa, 100)


rand(mcm, 10000)

c = stack(rand(mcm, 100000))
using Plots

plot(c)
histogram(c[1, :], nbins = 100)
histogram(c[2, :], nbins = 100)
histogram(c[3, :], nbins = 100)

plot(c[1, :] .^ 2 .+ c[2, :] .^ 2)
plot(c[2, :])

histogram

using StatsBase

h = fit(Histogram, (c[1, :], c[2, :]); nbins = 100)

plot(h)
using BenchmarkTools
rng = Random.default_rng()

@benchmark rand($rng, $mcm)

@code_warntype rand(rng, mcm)

@code_warntype mcmcaccept(rng, 0.1, mcm)

@benchmark mcmcaccept($rng, 0.1, $mcm)

using Profile

Profile.Allocs.clear()

@time Profile.Allocs.@profile sample_rate = 1 mcmcaccept(rng, 0.1, mcm)

using PProf

PProf.Allocs.pprof(from_c = false)


using Distributions

d = product_distribution((x = Uniform(), y = Uniform(), z = Uniform()))


using Random
Random.default_rng()


aa = NucleiWoodSaxon3D(100, 1.5, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
Participants(aa, aa, 1, 1, 6.4, 1, 0)

@code_warntype NucleiWoodSaxon3D(10, 1.5, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
@code_warntype Participants(aa, aa, 1, 1, 6.4, 1, 0)

pp = Participants(aa, aa, 1, 6.4, 1, 0)
@benchmark rand(pp, 2)
sdsd = Random.default_rng()
@code_warntype rand(sdsd, pp)

rand(pp)


rand(sdsd, pp)

asas = threarded(pp)


asdassdd = rand(asas, 1_000_000)
Participant

ncoll = map(asdassdd) do x
    x.n_coll
end

minimum(ncoll)

dump(asdassdd[findmin(ncoll)[2]])


function fluctuating_thickness!(out::V, part::T, grid_x::S, grid_y::R) where {V <: AbstractArray, T <: Participant, S <: AbstractVector, R <: AbstractVector}

    #part1=part.part1
    #part2=part.part2
    #shape1=part.shape1
    #shape2=part.shape2
    #w= part.sub_nucleon_width
    #p=part.p

    #function Tp(x,y)
    #    w2=2*w^2
    #    1/(pi*w2)*exp(-(x^2+y^2)/w2)
    #end

    @inbounds for i in eachindex(grid_x)
        x = grid_x[i]
        @inbounds for j in eachindex(grid_y)
            y = grid_y[j]

            #ta=sum(zip(part1,shape1)) do (pa ,ga)
            #    ga*Tp(x-pa[1],y-pa[2])
            #end

            #tb=sum(zip(part2,shape2)) do (pa ,ga)
            #    ga*Tp(x-pa[1],y-pa[2])
            #end

            #@show typeof(part2 )
            #tb=sum(zip(part2,shape2)) do pb ,gb
            #    gb*Tp(x-pb[1],y-pb[2])
            #end
            out[i, j] = part(x, y) #norm((ta,tb),p) #((ta^p+tb^p)/2)^(1/p)
        end
    end

    return nothing

end

function fluctuating_thickness_threads!(out::V, part::T, grid_x::S, grid_y::R) where {V <: AbstractArray, T <: Participant, S <: AbstractVector, R <: AbstractVector}

    #part1=part.part1
    #part2=part.part2
    #shape1=part.shape1
    #shape2=part.shape2
    #w= part.sub_nucleon_width
    #p=part.p

    #function Tp(x,y)
    #    w2=2*w^2
    #    1/(pi*w2)*exp(-(x^2+y^2)/w2)
    #end

    cart = CartesianIndices(out)
    tmap!(out, eachindex(out)) do idx
        c = cart[idx]
        x = grid_x[c[1]]
        y = grid_x[c[2]]
        part(x, y)
    end


    return nothing

end


function fluctuating_thickness!(out::V, part::T, grid_x::S) where {V <: AbstractArray, T <: Participant, S <: AbstractVector}

    return fluctuating_thickness!(out, part, grid_x, grid_x)
end

function fluctuating_thickness_threads!(out::V, part::T, grid_x::S) where {V <: AbstractArray, T <: Participant, S <: AbstractVector}

    return fluctuating_thickness_threads!(out, part, grid_x, grid_x)
end


function fluctuating_thickness(part::T, grid_x::S, grid_y::R) where {T <: Participant, S <: AbstractVector, R <: AbstractVector}
    out = Matrix{eltype(part)}(undef, length(grid_x), length(grid_y))

    fluctuating_thickness!(out, part, grid_x, grid_y)
    return out
end

function fluctuating_thickness_threads(part::T, grid_x::S, grid_y::R) where {T <: Participant, S <: AbstractVector, R <: AbstractVector}
    out = Matrix{eltype(part)}(undef, length(grid_x), length(grid_y))

    fluctuating_thickness_threads!(out, part, grid_x, grid_y)
    return out
end

function fluctuating_thickness(part::T, grid_x::S) where {T <: Participant, S <: AbstractVector}

    return fluctuating_thickness(part, grid_x, grid_x)

end

function fluctuating_thickness_threads(part::T, grid_x::S) where {T <: Participant, S <: AbstractVector}

    return fluctuating_thickness_threads(part, grid_x, grid_x)

end

function fluctuating_thickness(part::Vector{T}, grid_x::S, grid_y::R) where {T <: Participant, S <: AbstractVector, R <: AbstractVector}
    return map(x -> fluctuating_thickness(x, grid_x, grid_y), part)
end

function fluctuating_thickness_threads(part::Vector{T}, grid_x::S, grid_y::R) where {T <: Participant, S <: AbstractVector, R <: AbstractVector}
    return tmap(x -> fluctuating_thickness(x, grid_x, grid_y), part)
end

function fluctuating_thickness(part::Vector{T}, grid_x::S) where {T <: Participant, S <: AbstractVector}
    return map(x -> fluctuating_thickness(x, grid_x, grid_x), part)
end

function fluctuating_thickness_threads(part::Vector{T}, grid_x::S) where {T <: Participant, S <: AbstractVector}
    return tmap(x -> fluctuating_thickness(x, grid_x, grid_x), part)
end

mm = fluctuating_thickness(aba[1][2], -10:0.1:10, -10:0.1:10)
mm4 = fluctuating_thickness(aba[2][1], -10:0.1:10, -10:0.1:10)

heatmap(mm)
heatmap(mm4)

using HDF5
close(fid)
fid = h5open("test_profs.h5", "r")
tprof = read(fid)["event_1"]
heatmap(tprof)
p1 = heatmap(read(fid)["event_0"])
p2 = heatmap(read(fid)["event_2"])
p3 = heatmap(read(fid)["event_3"])
p4 = heatmap(read(fid)["event_4"])
p5 = heatmap(read(fid)["event_5"])
p6 = heatmap(read(fid)["event_6"])
p7 = heatmap(read(fid)["event_7"])
p8 = heatmap(read(fid)["event_8"])
plot(p1, p2, p3, p4, p5, p6, p7, p8, size = (900, 900))
keys(fid)

tprof["Event_1"]


using MonteCarloGlauber
using Statistics
using MonteCarloGlauber
using Plots

n1 = Lead()
n2 = Lead()
rand(n1)
rand(n2)

rand(n1, 100)

rand(threaded(n1), 100)

w = 1
s_NN = 5000
k = 1
p = 0.0
b = (1, 2)

participants = Participants(n1, n2, w, s_NN, k, p, b)

participants = threaded(Participants(n1, n2, w, s_NN, k, p))
event = rand(participants, 1000)

b_event = map(event) do x
    impactParameter(x)
end
histogram(b_event, nbins = 100)

ncoll_event = map(event) do x
    x.n_coll
end
histogram(ncoll_event, nbins = 100, yscale = :log10)

profile = map(event) do x
    map(Iterators.product(-10:0.5:10, -10:0.5:10)) do y
        x(y...)
    end
end
heatmap(-10:0.5:10, -10:0.5:10, profile[4])

b_event = map(1:1000) do x
    impactParameter(rand(participants))
end

energy2(T) = T
bg, twpt = generate_bg_two_pt_fct(energy2, energy2, 1, Lead(), Lead(), w, k, p, s_NN, [10, 20], [2]; minBiasEvents = 1000, r_grid = 0:1:10, nFields = 3, n_ext_Grid = 20)

plot(0:1:10, bg[1])
plot!(0:1:10, bg[2])

heatmap(twpt[1, 1, 1, :, :])

using HDF5

test = h5open("src/NLEFT_dmin_0.5fm_negativeweights_Ne.h5", "r")
configs = test["configs"]
size(configs)
close(test)
size(configs)

size(configs)

a = configs[:, :, 10]

using MonteCarloGlauber
aa = TabulatedEvent("NLEFT_dmin_0.5fm_negativeweights_Ne.h5")
aa.array
pp = Participants(aa, aa, 1, 1, 6, 1, 0)
event = rand(pp, 10)
profile = map(event) do x
    map(Iterators.product(-10:0.5:10, -10:0.5:10)) do y
        x(y...)
    end
end
heatmap(-10:0.5:10, -10:0.5:10, profile[1])

ab = Participants(Lead(), Lead(), 1, 1, 6, 1, 0)
rand(ab, 10)

rand(Lead(), 10)
MonteCarloGlauber.drawIndeces(1, aa)
rand(pp, 10)

using Rotations
rot = rand(RotMatrix{3})
ww = aa.array[:, :, 1]

rr = rot * ww

rr[2]
using MonteCarloGlauber
aa = NucleiWoodSaxon3D(10, 1.5, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
bb = IntegratedWoodSaxon(10, 1, 5.0, 1.0, 0.0)
cc = IntegratedWoodSaxonInterp(10, 1, 5.0, 1.0, 0.0)
dd = threaded(aa)
ee = Threaded(aa, Threads.nthreads())
ff = Threaded(aa, 10 * Threads.nthreads())

ll = Threaded(aa, 5)
ll.N_nucleon
aa.N_nucleon
ll.copy_buff[1].N_nucleon


using Test
size(aa) == (10, 2)
eltype(aa) == Float64
size(bb) == (10, 2)
eltype(bb) == Float64
size(cc) == (10, 2)
eltype(cc) == Float64
size(dd) == (10, 2)

rand(aa, 100)

event = Participants(aa, aa, 1, 1, 6.4, 1, 0)


typeof(rand(event))
ab = rand(event, 10)
ab[1].part1

rand(event, 100)

@test rand(event, 100) isa Vector{Array{Float64, 2}}
@test rand(event).part1 isa Vector{SVector{2, Float64}}

@test eltype(dd) == Float64
@test size(ee) == (10, 2)
@test eltype(ee) == Float64
@test size(ff) == (10, 2)
@test eltype(ff) == Float64


using Fluidum
using MonteCarloGlauber
using Plots

#INITIAL CONDITIONS
n1 = Lead()
n2 = Lead()
w = 0.5
s_NN = 5000.0
k = 1.0
p = 0.0
Norm = 50.0
eos = FluiduMEoS()
centrality_bins = [5, 10, 20, 30]
mlist = [2, 3]

entropy(T) = pressure_derivative(T, Val(1), FluiduMEoS()) #entropy as function of temperature
entropyToTemp(T) = InverseFunction(entropy)(T) #inverse, i.e. T(s)
dSdT(T) = pressure_derivative(T, Val(2), FluiduMEoS()) #function to convert perturbations as function of bg temp, i.e. dT/ds(T_0)
dSdTinverse(T) = 1 / dSdT(T)
#bg,twpt=MonteCarloGlauber.generate_bg_two_pt_fct_save_faster(entropyToTemp,dSdTinverse,Norm,n1,n2,w,k,p,s_NN,centrality_bins,mlist;minBiasEvents=1000,r_grid=0:0.2:15,n_ext_Grid=100,override_files=true,path="./")
@code_warntype entropy(1)
@code_warntype entropyToTemp(1)

entropy(1)

using QuadGK
using Statistics

function variance(configuration, r_1, r_2, m)
    return mean(configuration) do conf
        integral, er = quadgk(0, 2pi, rtol = 1.0e-3) do θ_1
            result, error = quadgk(0, 2pi, rtol = 1.0e-3) do θ_2
                x_1 = r_1 * cos(θ_1)
                y_1 = r_1 * sin(θ_1)
                x_2 = r_2 * cos(θ_2)
                y_2 = r_2 * sin(θ_2)
                conf(x_1, y_1) * conf(x_2, y_2) * exp(im * m * (θ_1 - θ_2))
            end
            return result
        end
        return integral
    end
end

function variance_at(configuration, r_1, r_2, m, len)
    return mean(configuration) do conf
        mean(range(0, 2pi, len)) do θ_1
            mean(range(0, 2pi, len)) do θ_2
                x_1 = r_1 * cos(θ_1)
                y_1 = r_1 * sin(θ_1)
                x_2 = r_2 * cos(θ_2)
                y_2 = r_2 * sin(θ_2)
                conf(x_1, y_1) * conf(x_2, y_2) * exp(im * m * (θ_1 - θ_2))
            end
        end
    end
end

function variance_at_2(configuration, r_1, r_2, m, len)
    result = zero(first(configuration)(r_1, r_2) * im)
    nevnet = length(configuration)
    θrange = range(0, 2pi, len)
    for i_conf in eachindex(configuration)
        sums = zero(result)
        for θ_1 in θrange
            x_1 = r_1 * cos(θ_1)
            y_1 = r_1 * sin(θ_1)
            c_1 = configuration[i_conf](x_1, y_1)
            for θ_2 in θrange
                x_2 = r_2 * cos(θ_2)
                y_2 = r_2 * sin(θ_2)
                c_2 = configuration[i_conf](x_2, y_2)
                sums += c_1 * c_2 * exp(im * m * (θ_1 - θ_2))
            end
        end
        result += sums
    end
    return result / (nevnet * len * len)
end

using Statistics

function mean_at(configuration, r_1, m, len)
    return mean(configuration) do conf
        mean(range(0, 2pi, len)) do θ_1
            x_1 = r_1 * cos(θ_1)
            y_1 = r_1 * sin(θ_1)
            conf(x_1, y_1)exp(im * m * (θ_1))
        end
    end
end


function cumulant_at(configuration, r_1, r_2, m, len)
    return variance_at(configuration, r_1, r_2, m, len) - mean_at(configuration, r_1, m, len) * conj(mean_at(configuration, r_2, m, len))
end

participants = Participants(n1, n2, w, s_NN, k, p)

ddd = rand(participants, 100)

variance_at(ddd, 1, 1, 0, 20)

variance_at(ddd, 2, 1, 0, 20)
variance_at(ddd, 1, 2, 0, 20)

variance_at(ddd, 1, 1, 0, 20)
variance_at_2(ddd, 1, 1, 0, 20)

rgrid = 0:1:10

resListBG = mean_at.(Ref(ddd), rgrid, Ref(0), Ref(20))

res = map(rgrid) do r1
    map(rgrid) do r2
        cumulant_at(ddd, r1, r2, 0, 20)
    end
end

heatmap(real.(res))

resi = stack(real.(res))
heatmap(resi)

size(resi)

plot(rgrid, real.(resListBG))

mean_at(ddd, 1, 0, 10)
cumulant_at(ddd, 1, 1, 0, 15)


using BenchmarkTools

@benchmark variance_at_2($ddd, 1, 1, 0, 15)

@code_warntype variance_at_2(ddd, 1, 1, 0, 15)

100 * 23 * 10_000 / 2 / 60 / 60

function second_cumulant(configuration, r_1, r_2, m, len)
    result = zero(first(configuration)(r_1, r_2)im)
    result_av = zero(result)
    nevnet = length(configuration)
    θrange = range(0, 2pi, len)
    for i_conf in eachindex(configuration)
        @inbounds conf = configuration[i_conf]
        sums = zero(result)
        sumaverege = zero(result)
        for θ_1 in θrange
            s1, c1 = sincos(θ_1)
            x_1 = r_1 * s1
            y_1 = r_1 * c1
            c_1 = conf(x_1, y_1)
            sfirst, cfirst = sincos(m * θ_1)
            #firstfactor=c_1*(sfirst+im*cfirst)
            firstfactor = c_1 * (cfirst + im * sfirst)
            sumaverege += firstfactor
            for θ_2 in θrange
                s2, c2 = sincos(θ_2)
                x_2 = r_2 * s2
                y_2 = r_2 * c2
                c_2 = conf(x_2, y_2)
                sdiff, cdiff = sincos(m * (θ_2))
                sums += c_2 * (cdiff - im * sdiff) * firstfactor
            end
        end
        result += sums
        result_av += sumaverege
    end
    return (result - abs2(result_av) / nevnet) / (nevnet * len * len)
end

ddd = rand(participants, 500)
second_cumulant(ddd, 1.2, 1, 0, 20)
cumulant_at(ddd, 1.2, 1, 0, 20)

using BenchmarkTools
@benchmark second_cumulant(ddd, 2, 1, 0, 15)
@benchmark cumulant_at(ddd, 2, 1, 0, 15)

rgrid = 0:1:10

res = map(rgrid) do r1
    map(rgrid) do r2
        second_cumulant(ddd, r1, r2, 3, 20)
    end
end

second_cumulant(ddd, 2, 1, 2, 5)
second_cumulant(ddd, 2, 1, 2, 10)
second_cumulant(ddd, 2, 1, 2, 15)
second_cumulant(ddd, 2, 1, 2, 20)
second_cumulant(ddd, 2, 1, 2, 25)
second_cumulant(ddd, 2, 1, 2, 30)
second_cumulant(ddd, 2, 1, 2, 60)
second_cumulant(ddd, 2, 1, 2, 100)


second_cumulant(ddd, 1, 1, 2, 5)
second_cumulant(ddd, 1, 1, 2, 10)
second_cumulant(ddd, 1, 1, 2, 15)
second_cumulant(ddd, 1, 1, 2, 20)
second_cumulant(ddd, 1, 1, 2, 25)
second_cumulant(ddd, 1, 1, 2, 30)
second_cumulant(ddd, 1, 1, 2, 60)
second_cumulant(ddd, 1, 1, 2, 400)

0.02039 / 0.0212916

resi = stack(real.(res))
heatmap(resi)

function second_cumulant_old(configuration, r_1, r_2, m, norm, len)
    result = zero(first(configuration)(r_1, r_2)im)
    result_av = zero(result)
    nevnet = length(configuration)
    θrange = range(0, 2pi, len)
    for i_conf in eachindex(configuration)
        @inbounds conf = configuration[i_conf]
        sums = zero(result)
        sumaverege = zero(result)
        for θ_1 in θrange
            s1, c1 = sincos(θ_1)
            x_1 = r_1 * c1
            y_1 = r_1 * s1
            c_1 = norm * conf(x_1, y_1)
            sfirst, cfirst = sincos(m * θ_1)
            #firstfactor=c_1*(sfirst+im*cfirst)
            firstfactor = c_1 * (cfirst + im * sfirst)
            sumaverege += firstfactor #<s^m(r)>
            for θ_2 in θrange
                s2, c2 = sincos(θ_2)
                x_2 = r_2 * c2
                y_2 = r_2 * s2
                c_2 = norm * conf(x_2, y_2)
                sdiff, cdiff = sincos(m * (θ_2))
                sums += c_2 * (cdiff - im * sdiff) * firstfactor # int phi_i phi_j entropy(ri,phii)*entropy(r_j,phij_)*exp(i*m*(phi_i-phi_j))
            end
        end
        result += sums
        result_av += sumaverege
    end
    return (result - abs2(result_av) / nevnet) / (nevnet * len * len)
    #return (result-result_av*result_av_cc/nevnet)/(nevnet*len*len)
end

function second_cumulant_new(configuration, r_1, r_2, m, norm, len)
    result = zero(first(configuration)(r_1, r_2)im)
    result_av = zero(result)
    result_av_cc = zero(result)
    nevnet = length(configuration)
    θrange = range(0, 2pi, len)
    for i_conf in eachindex(configuration)
        @inbounds conf = configuration[i_conf]
        sums = zero(result)
        sumaverege = zero(result)
        sumaverage_cc = zero(result)
        for θ_1 in θrange
            s1, c1 = sincos(θ_1)
            x_1 = r_1 * c1
            y_1 = r_1 * s1
            c_1 = norm * conf(x_1, y_1)
            sfirst, cfirst = sincos(m * θ_1)
            #firstfactor=c_1*(sfirst+im*cfirst)
            firstfactor = c_1 * (cfirst + im * sfirst)
            sumaverege += firstfactor #<s^m(r)>
            for θ_2 in θrange
                s2, c2 = sincos(θ_2)
                x_2 = r_2 * c2
                y_2 = r_2 * s2
                c_2 = norm * conf(x_2, y_2)
                sdiff, cdiff = sincos(m * (θ_2))
                sums += c_2 * (cdiff - im * sdiff) * firstfactor # int phi_i phi_j entropy(ri,phii)*entropy(r_j,phij_)*exp(i*m*(phi_i-phi_j))
                sumaverage_cc = c_2 * (cdiff - im * sdiff)
            end
        end
        result += sums
        result_av += sumaverege
        result_av_cc += sumaverage_cc
    end
    #return (result-abs2(result_av)/nevnet)/(nevnet*len*len)
    return (result - result_av * result_av_cc / nevnet) / (nevnet * len * len)
end

#<ss>-<s><s>= sum ss /N- sum

r_grid = 0:1:20

oldList = second_cumulant_old.(Ref(ddd), Ref(1), r_grid, Ref(2), Ref(1), Ref(10))
newList = second_cumulant_new.(Ref(ddd), Ref(1), r_grid, Ref(2), Ref(1), Ref(10))

plot(real.(oldList), label = "old")
plot!(real.(newList), label = "new")

second_cumulant_new(ddd, 1.2, 1, 2, 1, 20)
second_cumulant_old(ddd, 1.2, 1, 2, 1, 20)

second_cumulant_new(ddd, 1, 1, 2, 1, 20)
second_cumulant_old(ddd, 1, 1, 2, 1, 20)

second_cumulant_new(ddd, 15, 0.7, 2, 1, 20)
second_cumulant_old(ddd, 15, 0.7, 2, 1, 20)


second_cumulant(ddd, 1.2, 1, 0, 1, 10)
second_cumulant(ddd, 1.2, 1, 0, 1, 200)
second_cumulant(ddd, 1, 1, 0, 1, 20)

using MonteCarloGlauber

participants = Participants(n1, n2, w, s_NN, k, p)


ddd = rand(participants, 10)

MonteCarloGlauber.center_of_mass_gl(ddd[1])
MonteCarloGlauber.center_of_mass_gl(ddd[2])
MonteCarloGlauber.center_of_mass_gl(ddd[3])
MonteCarloGlauber.center_of_mass_gl(ddd[4])

aa = MonteCarloGlauber.prepare_accumulation(ddd[2])
aa2 = MonteCarloGlauber.prepare_accumulation2(ddd[1])
MonteCarloGlauber.center_of_mass_gl!(ddd[2], aa)

0.029 / 0.104

MonteCarloGlauber.center_of_mass(ddd[2])
0.032 / 0.108

zeros(eltype(rgrid), length(rgrid), length(rgrid))

function generate_tw_pt_fct(configuration, r_grid, m, norm, len)
    correlator = zeros(eltype(r_grid), length(r_grid), length(r_grid))
    for i1 in 1:length(r_grid)
        r1 = r_grid[i1]
        for i2 in i1:length(r_grid)
            r2 = r_grid[i2]
            correlator[i1, i2] = real(second_cumulant_new(configuration, r1, r2, m, norm, len))
            correlator[i2, i1] = correlator[i1, i2]
        end
    end
    return correlator
end

#dims final correlator
#cc, reim, field, field, m ,r1,r2


#1 sort profiles in cc
#2 get bg and two fct with norm
#3 apply Eos(bg) and dt/ds*twoPtFct
#4 put in matrix

using MonteCarloGlauber
using Plots


participants = Participants(n1, n2, w, s_NN, k, p)

ddd = rand(participants, 100)
batches = MonteCarloGlauber.centralities_selection_events(ddd, [10, 20])

function generate_bg(batches, bins, r_grid, InverseEoS, Norm; NumPhiPoints = 20, Threaded = true)
    bg = zeros(eltype(r_grid), length(bins), length(r_grid))
    #batches=MonteCarloGlauber.centralities_selection_events(events,bins;Threaded=Threaded)
    for cc_batches in 1:(length(batches) - 1)
        for r_i in eachindex(r_grid)
            r = r_grid[r_i]
            #@show cc_batches r_i
            bg[cc_batches, r_i] = real.(mean_at(batches[cc_batches], r, 0, NumPhiPoints))

        end
    end
    return InverseEoS.(Norm .* bg)

end

bg_new = generate_bg(batches, [10, 20], 0.0:1:20, entropyToTemp, 10)

@code_warntype generate_bg(batches, [10, 20], 0.0:1:20, entropyToTemp, 10)
@code_warntype entropyToTemp(1.0)

plot(bg_new[1, :])
plot!(bg_new[2, :])

last(bg_new[1, :])
last(bg_new[2, :])

aaa = zeros(3, 5, 5, 6, 7)
aaa
size(aaa)[4]

function generate_tw_pt_fct_entropy(batches, bins, r_grid, m_list, Norm; NumPhiPoints = 20, Threaded = true, Nfields = 10)
    finalCorrelator = zeros(eltype(r_grid), length(bins), 2, Nfields, Nfields, length(m_list), length(r_grid), length(r_grid))
    #batches=MonteCarloGlauber.centralities_selection_events(events,bins;Threaded=Threaded)
    for cc in 1:(length(batches) - 1)
        for m in 1:length(m_list)
            for r1 in 1:length(r_grid)
                for r2 in r1:length(r_grid)
                    finalCorrelator[cc, 1, 1, 1, m, r1, r2] = real(second_cumulant_new(batches[cc], r_grid[r1], r_grid[r2], m_list[m], Norm, NumPhiPoints))
                    finalCorrelator[cc, 1, 1, 1, m, r2, r1] = finalCorrelator[cc, 1, 1, 1, m, r1, r2]
                end
            end
        end
    end
    return finalCorrelator
end


function generate_bg_twpt_fct(f, delta_factor, norm, Projectile1, Projectile2, w, k, p, sqrtS, bins, mList; minBiasEvents = 1000000, r_grid = 0.0:1:10, NumPhiPoints = 20, Threaded = true, Nfields = 10)
    correlator = zeros(eltype(r_grid), length(bins), 2, Nfields, Nfields, length(mList), length(r_grid), length(r_grid))
    participants = Participants(Projectile1, Projectile2, w, sqrtS, k, p)
    if Threaded
        events = rand(threaded(participants), minBiasEvents)
    else
        events = rand(participants, minBiasEvents)
    end
    temp_shift_vector = MonteCarloGlauber.shift_and_compute.(events)
    new_event = getindex.(temp_shift_vector, 1)
    mult = getindex.(temp_shift_vector, 2)

    batches = MonteCarloGlauber.centralities_selection_events(new_event, bins, mult; Threaded = Threaded)
    bg = generate_bg(batches, bins, r_grid, f, norm; NumPhiPoints = NumPhiPoints, Threaded = Threaded)
    tw_pt_entropy = generate_tw_pt_fct_entropy(batches, bins, r_grid, mList, norm; NumPhiPoints = NumPhiPoints, Nfields = Nfields)

    for cc in 1:(length(batches) - 1) #TODO put this in own function
        for m in 1:length(mList)
            for r1 in 1:length(r_grid)
                for r2 in r1:length(r_grid)
                    correlator[cc, 1, 1, 1, m, r1, r2] = tw_pt_entropy[cc, 1, 1, 1, m, r1, r2] * delta_factor(bg[cc, r1]) * delta_factor(bg[cc, r2])
                    correlator[cc, 1, 1, 1, m, r2, r1] = correlator[cc, 1, 1, 1, m, r1, r2]
                end
            end
        end
    end

    return bg, correlator
end
@code_warntype generate_bg_twpt_fct(entropyToTemp, dSdTinverse, 30, n1, n2, w, 10, p, s_NN, [10], [2, 3]; minBiasEvents = 500, r_grid = 0.0:1:20)

using Statistics
k
using MonteCarloGlauber

using Fluidum
using MonteCarloGlauber
using Plots

#INITIAL CONDITIONS
n1 = Lead()
n2 = Lead()
w = 0.5
s_NN = 5000.0
k = 1.0
p = 0.0
Norm = 50.0
eos = FluiduMEoS()
centrality_bins = [5, 10, 20, 30]
mlist = [2, 3]

entropy(T) = T #pressure_derivative(T,Val(1),FluiduMEoS()) #entropy as function of temperature
entropyToTemp(T) = T #InverseFunction(entropy)(T) #inverse, i.e. T(s)
dSdT(T) = T #pressure_derivative(T,Val(2),FluiduMEoS()) #function to convert perturbations as function of bg temp, i.e. dT/ds(T_0)
dSdTinverse(T) = T #1/dSdT(T)


a10, b10 = generate_bg_twpt_fct(entropyToTemp, dSdTinverse, 30, n1, n2, w, 10, p, s_NN, [10], [2, 3]; minBiasEvents = 100, r_grid = 0.0:1:20)


a1, b1 = generate_bg_twpt_fct(entropyToTemp, dSdTinverse, 30, n1, n2, w, 1, p, s_NN, [10], [2, 3]; minBiasEvents = 500, r_grid = 0.0:1:20)
a01, b01 = generate_bg_twpt_fct(entropyToTemp, dSdTinverse, 30, n1, n2, w, 0.1, p, s_NN, [10], [2, 3]; minBiasEvents = 500, r_grid = 0.0:1:20)

a10

heatmap(b10[1, 1, 1, 1, 1, :, :])
plot(b10[1, 1, 1, 1, 1, 10, :])
heatmap(b10[1, 1, 1, 1, 1, 1:9, 1:9])

minimum(b10[1, 1, 1, 1, 1, :, :])

@code_warntype generate_bg_twpt_fct(entropyToTemp, dSdTinverse, 30, n1, n2, w, 10, p, s_NN, [10], [2, 3]; minBiasEvents = 500, r_grid = 0.0:1:20)
@benchmark generate_bg_twpt_fct($entropyToTemp, $dSdTinverse, 30, n1, n2, w, 10, p, s_NN, [10], [2, 3]; minBiasEvents = 500, r_grid = 0.0:1:20)

plot(dSdTinverse.(a10[1, :] .+ 0.01))

exponential_tail_pointlike(function_profile, x; xmax = 8, offset = 0.015)
temperature_funct = linear_interpolation(r, temperature_profile; extrapolation_bc = Flat()) #from 9.9 fm it's flat

temp_exp = exponential_tail_pointlike.(Ref(temperature_funct), radius; xmax = x_max_x, offset = offset_x)

using Interpolations
r_grid = range(0.0, 20)
a10
bgfun = Fluidum.linear_interpolation(r_grid, a10[1, :]; extrapolation_bc = Flat())


bgfun(2.23)

bg_good = Fluidum.exponential_tail_pointlike.(Ref(bgfun), r_grid; xmax = 10, offset = 0.015)

function addTailBG(bg, r_grid; xmax = 10, offset = 0.015)
    bgfun = Fluidum.linear_interpolation(r_grid, bg; extrapolation_bc = Flat())
    return Fluidum.exponential_tail_pointlike.(Ref(bgfun), r_grid; xmax = xmax, offset = offset)
end

addTailBG(a10[1, :], r_grid)

plot(r_grid, bg_good)
plot!(r_grid, a10[1, :])

plot(r_grid, dSdTinverse.(bg_good))
plot!(r_grid, dSdTinverse.(a10[1, :]))


plot(a10[1, :])
plot(dSdTinverse.(a10[1, :]))
minimum(a10[1, :])

# small k -> large flucts
sqrt(maximum(b10[1, 1, 1, 1, 1, :, :])) / maximum(a10[1, :])
sqrt(maximum(b1[1, 1, 1, 1, 1, :, :])) / maximum(a1[1, :])
sqrt(maximum(b01[1, 1, 1, 1, 1, :, :])) / maximum(a01[1, :])

sqrt(b10[1, 1, 1, 1, 2, 5, 5]) / a10[1, 5]
sqrt(b1[1, 1, 1, 1, 2, 5, 5]) / a1[1, 5]
sqrt(b01[1, 1, 1, 1, 2, 5, 5]) / a01[1, 5]


1 / sqrt(50)

sqrt(-minimum(b10[1, 1, 1, 1, 1, :, :])) / maximum(a10[1, :])
sqrt(-minimum(b1[1, 1, 1, 1, 1, :, :])) / maximum(a1[1, :])
sqrt(-minimum(b01[1, 1, 1, 1, 1, :, :])) / maximum(a01[1, :])

plot(a10[1, :])
plot!(a1[1, :])
plot!(a01[1, :])

plot(heatmap(b[1, 1, 1, 1, 1, 1:30, 1:30]), heatmap(b[2, 1, 1, 1, 1, 1:30, 1:30]), heatmap(b[1, 1, 1, 1, 2, 1:30, 1:30]), heatmap(b[2, 1, 1, 1, 2, 1:30, 1:30]))

heatmap(b10[1, 1, 1, 1, 1, :, :])


plot(b[1, 1, 1, 1, 1, 1, 1:20])
plot!(a[1, :])

plot(a[1, :])
plot(dSdTinverse.(a[1, :]))

sqrt(maximum(b[1, 1, 1, 1, 1, 1:30, 1:30]))
maximum(a[1, :])

sqrt(maximum(b[1, 1, 1, 1, 1, 1:30, 1:30])) / maximum(a[1, :])


batches
second_cumulant(ddd, 1, 1, 2, 100, 20)

fun(T) = T
length(0.0:0.1:10)
ddd = rand(participants, 100)
a = generate_bg(batches, [10, 20], 0.0:1:10, fun, 2)
a[1, 2]
plot(a[1, :])
plot!(a[2, :])
plot!(a[3, :])

b = generate_tw_pt_fct_entropy(ddd, [10, 20], 0.0:0.1:10, [2, 3], 10)

heatmap(b[1, 1, 1, 1, 2, :, :])

using BenchmarkTools
@code_warntype generate_bg(ddd, [10, 20], 0.0:0.1:10, 1)
@benchmark generate_bg($ddd, [10, 20], 0.0:0.1:10, 1)

@code_warntype generate_tw_pt_fct_entropy(ddd, [10, 20], 0.0:1:5, [2, 3], 10)
@benchmark generate_tw_pt_fct_entropy(ddd, [10, 20], 0.0:1:5, [2, 3], 10)


@code_warntype generate_tw_pt_fct(ddd, 0.0:1:3, 2, 1, 20)
@code_warntype generate_tw_pt_fct(ddd, 0.0:1:6, 2, 1, 20)
@benchmark generate_tw_pt_fct(ddd, 0.0:1:3, 2, 1, 20)
@benchmark generate_tw_pt_fct(ddd, 0.0:1:6, 2, 1, 20)

length(0:1:3)
length(0:1:6)

length(ddd)

4 * 5 / 2
7 * 8 / 2


ddd[1]

MonteCarloGlauber.center_of_mass(ddd[1])
0.03 / 0.108

a, b = MonteCarloGlauber.centralities_selection_CoM(ddd, [0, 10])

a
b[1]

b[1][1][2] / b[1][1][1]
b[1][1][3] / b[1][1][1]
ddd[1](1, 1)

MonteCarloGlauber.fluctuating_thickness(1, 1, ddd[1])

b[1][10][2] / b[1][10][1]
b[1][10][3] / b[1][10][1]

a = generate_tw_pt_fct(ddd, 0:0.2:10.0, 2, 1, 20)

heatmap(a)

event = rand(participants, 10)
profile = map(event) do x
    map(Iterators.product(-10:0.5:10, -10:0.5:10)) do y
        x(y...)
    end
end
heatmap(-10:0.5:10, -10:0.5:10, profile[7])
