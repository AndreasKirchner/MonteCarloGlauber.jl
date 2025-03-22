using QuadGK
using Distributions
using Random
using Base.Threads: nthreads
using OhMyThreads: @spawn, fetch
using OhMyThreads
using Rotations
using LinearAlgebra
using Interpolations
using StaticArrays
using Plots
using NonlinearSolve
using SpecialFunctions

include("Metropolis_Hastings.jl")
include("Integrated_wood_sakson_refactor.jl")
include("Integrated_wood_sakson_interp_refactor.jl")
include("Nucleos3D_refactor.jl")
include("participants_distribution_refactor.jl")

using MonteCarloGlauber


using Profile

using PProf
Profile.clear()
Profile.Allocs.clear()
@profile rand(aa,10_000_000)
Profile.Allocs.@profile rand(aa,100)
pprof()
rand(threarded(aa),1000)
prof = Profile.Allocs.fetch()
PProf.Allocs.pprof(prof; from_c=false)

@benchmark rand($aa,10)

aa=NucleiWoodSaxon3D(10, 1.5, 5., 1., 0., 0., 0., 0., 0.)
aa=IntegratedWoodSaxon(10,1,5.,1.,0.)
aa=IntegratedWoodSaxonInterp(10,1,5.,1.,0.)
bb=Threarded(aa,2*nthreads())
cc=Threarded(aa,nthreads())
dd=Threarded(aa,10*nthreads())

using BenchmarkTools
@benchmark rand($aa,10000)
@benchmark rand($bb,10000)
@benchmark rand($cc,10000)
@benchmark rand($dd,10000)

Normal(0,1)

TriangularDist(0, 1) 
@code_warntype rand(Random.default_rng() , aa)
histogram(rand(TriangularDist(1, 1,1),100000))


rand(aa,100)



rand(mcm,10000)

c=stack(rand(mcm,100000))
using Plots

plot(c)
histogram(c[1,:],nbins=100)
histogram(c[2,:],nbins=100)
histogram(c[3,:],nbins=100)

plot(c[1,:].^2 .+c[2,:].^2)
plot(c[2,:])

histogram

using StatsBase

h = fit(Histogram, (c[1,:],c[2,:]);nbins=100)

plot(h)
using BenchmarkTools
rng=Random.default_rng()

@benchmark rand($rng,$mcm)

@code_warntype rand(rng,mcm)

@code_warntype mcmcaccept(rng,0.1,mcm)

@benchmark mcmcaccept($rng,0.1,$mcm)

using Profile

Profile.Allocs.clear()

@time Profile.Allocs.@profile sample_rate=1 mcmcaccept(rng,0.1,mcm)

using PProf

PProf.Allocs.pprof(from_c=false)


using Distributions

d = product_distribution((x=Uniform(),y=Uniform(),z=Uniform()))




using Random
Random.default_rng()




aa=NucleiWoodSaxon3D(100, 1.5, 5., 1., 0., 0., 0., 0., 0.)
Participants(aa,aa,1,1,6.4,1,0)

@code_warntype NucleiWoodSaxon3D(10, 1.5, 5., 1., 0., 0., 0., 0., 0.)
@code_warntype Participants(aa,aa,1,1,6.4,1,0)

pp=Participants(aa,aa,1,6.4,1,0)
@benchmark rand(pp,2)
sdsd=Random.default_rng()
@code_warntype rand(sdsd,pp,)

rand(pp)


rand(sdsd,pp)

asas=threarded(pp)


asdassdd=rand(asas,1_000_000)
Participant

ncoll=map(asdassdd) do x
x.n_coll    
end 

minimum(ncoll)

dump(asdassdd[findmin(ncoll)[2]])