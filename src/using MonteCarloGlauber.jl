
using Statistics
using MonteCarloGlauber
using StatsBase
using Plots
using StaticArrays
using OhMyThreads
using Cuba
using FastGaussQuadrature
 


n1= Lead()
n2= Uranium()

rand(n1)

rand(n2)

rand(threaded(n1),100)

w= 1
s_NN=5000
k=1
p=0.
b=(10,11)

participants=Participants(n1,n2,2,s_NN,k,p)

event=rand(threaded(participants),50_000)

b_event=map(event) do x
    impactParameter(x) 
end 

ncoll_event=map(event) do x
    x.n_coll
end


profile=map(event)   do x 
    map(Iterators.product(-10:10,-10:10)) do y
        x(y...)
    end
end 

multi=tmap(event)   do x 
    multiplicity(x)
end


histogram(multi,nbins=100,normalize=true,yscale=:log10)

obj=fit(Histogram,multi)

obj.edges


profile=mapreduce(+,1:10000)   do _ 
    x=rand(participants) 
    map(Iterators.product(-10:10,-10:10)) do y
        x(y...)
    end
end



rand(threaded(participants))
rand(threaded(participants),10)







using BenchmarkTools

one_profile=rand(participants)
multiplicity(one_profile)
@code_warntype multiplicity(one_profile)
@benchmark multiplicity($one_profile)

using Integrals
function multiplicity_hcubature(con::T) where {T<:Participant}

f(u, p) = SVector{3}(1,u[1],u[2]) .*con(u[1],u[2])
domain = (SVector{2}(-Inf, -Inf), SVector{2}(Inf, Inf))
    prob = IntegralProblem{false}(f, domain)
    sol = solve(prob, GaussLegendre(); reltol = 1e-3, abstol = 1e-3)
    sol.u


end 

@benchmark  multiplicity_hcubature($one_profile)



heatmap(profile)


Plots.histogram(b_event,nbins=100)
Plots.histogram(ncoll_event,nbins=100)


participants=Participants(n1,n2,w,s_NN,k,p)

event=rand(participants,1000)

b_event=map(event) do x
    impactParameter(x) 
end 

ncoll_event=map(event) do x
    x.n_coll
end 

using Plots
Plots.histogram(b_event,nbins=100)
Plots.histogram(ncoll_event,nbins=100)

participants=(Participants(n1,n2,w,s_NN,k,p))

event=rand(participants,1000)

b_event=map(event) do x
    impactParameter(x) 
end 

ncoll_event=map(event) do x
    x.n_coll
end 

using Plots
Plots.histogram(b_event,nbins=100)
Plots.histogram(ncoll_event,nbins=100)




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
bb=threarded(aa)
cc=Threarded(aa,Threads.nthreads())
dd=Threarded(aa,10*Threads.nthreads())

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