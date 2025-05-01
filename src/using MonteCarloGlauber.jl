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
n1= Lead()
n2= Lead()
w= 0.5
s_NN=2760
k=1
p=1.


bg,twpt=MonteCarloGlauber.generate_bg_two_pt_fct(energy2,energy2,1,Lead(),Lead(),w,k,p,s_NN,[10,20],[2,3,4];minBiasEvents=5000,r_grid=0:1:10)

MonteCarloGlauber.save_bg_two_pt_fct(energy2,energy2,1,Lead(),Lead(),w,k,p,s_NN,[10,20],[2,3,4];minBiasEvents=1000,r_grid=0:1:4)
using DelimitedFiles
aa=readdlm("trento_two_pt_fct_m_2_Pb_Pb_w_0_5_sNN_6_206_p_1_0_k_1.dat")

aa

aa[1]

plot(bg[1],label="0-10")
plot!(bg[2],label="10-20")
heatmap(real.(permutedims(hcat(twpt[1][1]...))))

MonteCarloGlauber.cross_section_from_energy(s_NN)

participants=Participants(n1,n2,w,s_NN,k,p,(0.1,6.8))
event=rand(threaded(participants),500)
aba, acom=MonteCarloGlauber.centralities_selection_CoM(event,[10])
energy(T)=3*8*pi^2/15 *(3^2 -1 +7/4* 3*3)*T^4/(2*3*4)
Te4=MonteCarloGlauber.InverseFunction(energy)

bb=MonteCarloGlauber.generate_background(Te4,1,aba,acom;r_grid=0:0.1:10)
plot(bb)
corr=MonteCarloGlauber.generate_2ptfct(Te4,1,aba,acom,[2,3,4];r_grid=0:0.5:5)
heatmap(real.(corr[1][1]))
bg

writedlm("bg.dat",bg)
aaa=readdlm("bg.dat")


size(corr)
size(corr[1])
size(corr[1][1])
corr[1][1][1][1]
length(corr)

shape(corr)

heatmap(imag.(corr[1][1]))
bb[1][1]

using MonteCarloGlauber
bg,twpt=MonteCarloGlauber.generate_bg_two_pt_fct(Te,Te,1,Lead(),Lead(),w,k,p,s_NN,[10,20],[2,3,4];minBiasEvents=5000,r_grid=0:0.5:5)

heatmap(real.(permutedims(hcat(twpt[1][1]...))))



b_event=map(event) do x
    impactParameter(x) 
end 
histogram(b_event)

multi=tmap(event)   do x 
    multiplicity(x,Nr=50,Nth=15)[1]
end

multi=tmap(event)   do x 
    multiplicity(x,Nr=50,Nth=15)
end

using DelimitedFiles
dat=readdlm("BKG276TeVAvgCC1.txt")

plot(dat[:,2])

aba, acom=MonteCarloGlauber.centralities_selection_CoM(event,[10])

energy(T)=3*8*pi^2/15 *(3^2 -1 +7/4* 3*3)*T^4/(2*3*4)
Te4=MonteCarloGlauber.InverseFunction(energy)
energy(1)
energy2(T)=T
Te=InverseFunction(energy2)

energy.([1,2,3])
Te4(1)
Te.([1,2,3])

bb=MonteCarloGlauber.generate_background(Te4,1,aba,acom;r_grid=0:1:10)
plot(bb,label="new")
plot!(dat[:,2],label="old")
plot!(1.5 .*dat[:,2],label="old*1.5")

sum(dat[:,2])*0.1

mean(i->aba[1][i](0,0),1:length(aba[1]))

bb
bbb=acom[1]

bbb[1][2]

for (index, value) in mean(enumerate(aba[1]))
    index
end

aba[1][1](0,-1)

tl=("a","b","c","d")

sum(enumerate(aba[1])) do (i,conf)
    i
end

mean(enumerate(aba[1])) do (i,conf)
    i
end


1+2+3+4

enumerate(tl) do i 
    i
end

bs=for (i,conf) in enumerate(tl)
   i
end

bs

grid=0:1:4
MonteCarloGlauber.generatingfunction(aba[1],1,0:1:3,0.5,acom[1],Te4,1)
bb=MonteCarloGlauber.generate_background(Te4,1,aba,acom;r_grid=0:1:3)
aba[1]
acom
####this is the plot 
histogram(multi,nbins=100,normalize=true,yscale=:log10)



obj=fit(Histogram,multi;nbins=100)

dump(obj)


profile=mapreduce(+,1:10000)   do _ 
    x=rand(participants) 
    map(Iterators.product(-10:10,-10:10)) do y
        x(y...)
    end
end



rand(threaded(participants))
rand(threaded(participants),10)


ncoll_event=map(event) do x
    x.n_coll
end

rand()

using Distributions

a1=10*sqrt.(rand(100000))
histogram(a1,bins=100)
aa=truncated(TriangularDist(0,10,10),0,10)
a2=rand(aa,100000)
histogram(a1,bins=100)
histogram!(a2,bins=100)


x1=rand(Uniform(-10,10),100000)
x2=rand(Uniform(-10,10),100000)
b=sqrt.(x1 .^2 .+ x2 .^2)
histogram!(b,normalize=true)
histogram(a1,bins=100,normalize=true)
histogram!(a2,bins=100,normalize=true)

dd=sqrt.(Uniform(0,1))
dd=Uniform(0,1)

rand(dd,12)

n1.R

profile=map(event)   do x 
    map(Iterators.product(-10:10,-10:10)) do y
        x(y...)
    end
end 





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


function fluctuating_thickness!(out::V,part::T,grid_x::S,grid_y::R) where {V<:AbstractArray,T<:Participant,S<:AbstractVector,R<:AbstractVector}
    
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
        x=grid_x[i]
        @inbounds for j in eachindex(grid_y)
            y=grid_y[j]
            
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
            out[i,j]= part(x,y)#norm((ta,tb),p) #((ta^p+tb^p)/2)^(1/p)
        end 
    end

    return nothing

end

function fluctuating_thickness_threads!(out::V,part::T,grid_x::S,grid_y::R) where {V<:AbstractArray,T<:Participant,S<:AbstractVector,R<:AbstractVector}
    
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
        x=grid_x[c[1]]
        y=grid_x[c[2]]
        part(x,y)
    end
    

    return nothing

end



function fluctuating_thickness!(out::V,part::T,grid_x::S) where {V<:AbstractArray,T<:Participant,S<:AbstractVector}
    
    fluctuating_thickness!(out,part,grid_x,grid_x) 
end

function fluctuating_thickness_threads!(out::V,part::T,grid_x::S) where {V<:AbstractArray,T<:Participant,S<:AbstractVector}
    
    fluctuating_thickness_threads!(out,part,grid_x,grid_x) 
end


function fluctuating_thickness(part::T,grid_x::S,grid_y::R) where {T<:Participant,S<:AbstractVector,R<:AbstractVector}
    out=Matrix{eltype(part)}(undef,length(grid_x),length(grid_y))
    
    fluctuating_thickness!(out,part,grid_x,grid_y)
    return out 
end

function fluctuating_thickness_threads(part::T,grid_x::S,grid_y::R) where {T<:Participant,S<:AbstractVector,R<:AbstractVector}
    out=Matrix{eltype(part)}(undef,length(grid_x),length(grid_y))
    
    fluctuating_thickness_threads!(out,part,grid_x,grid_y)
    return out 
end

function fluctuating_thickness(part::T,grid_x::S) where {T<:Participant,S<:AbstractVector}
    
    fluctuating_thickness(part,grid_x,grid_x)
    
end

function fluctuating_thickness_threads(part::T,grid_x::S) where {T<:Participant,S<:AbstractVector}
    
    fluctuating_thickness_threads(part,grid_x,grid_x)
    
end

function fluctuating_thickness(part::Vector{T},grid_x::S,grid_y::R)where {T<:Participant,S<:AbstractVector,R<:AbstractVector}
    map(x->fluctuating_thickness(x,grid_x,grid_y),part)
end 

function fluctuating_thickness_threads(part::Vector{T},grid_x::S,grid_y::R)where {T<:Participant,S<:AbstractVector,R<:AbstractVector}
    tmap(x->fluctuating_thickness(x,grid_x,grid_y),part)
end 

function fluctuating_thickness(part::Vector{T},grid_x::S)where {T<:Participant,S<:AbstractVector}
    map(x->fluctuating_thickness(x,grid_x,grid_x),part)
end

function fluctuating_thickness_threads(part::Vector{T},grid_x::S)where {T<:Participant,S<:AbstractVector}
    tmap(x->fluctuating_thickness(x,grid_x,grid_x),part)
end

mm=fluctuating_thickness(aba[1][2],-10:0.1:10,-10:0.1:10)
mm4=fluctuating_thickness(aba[2][1],-10:0.1:10,-10:0.1:10)

heatmap(mm)
heatmap(mm4)

using HDF5
close(fid)
fid = h5open("test_profs.h5", "r")
tprof=read(fid)["event_1"]
heatmap(tprof)
p1=heatmap(read(fid)["event_0"])
p2=heatmap(read(fid)["event_2"])
p3=heatmap(read(fid)["event_3"])
p4=heatmap(read(fid)["event_4"])
p5=heatmap(read(fid)["event_5"])
p6=heatmap(read(fid)["event_6"])
p7=heatmap(read(fid)["event_7"])
p8=heatmap(read(fid)["event_8"])
plot(p1,p2,p3,p4,p5,p6,p7,p8,size=(900,900))
keys(fid)

tprof["Event_1"]