using MonteCarloGlauber
using Test
#using StaticArrays

@testset "MonteCarloGlauber.jl" begin
    aa=NucleiWoodSaxon3D(10, 1.5, 5., 1., 0., 0., 0., 0., 0.)
    bb=IntegratedWoodSaxon(10,1,5.,1.,0.)
    cc=IntegratedWoodSaxonInterp(10,1,5.,1.,0.)
    dd=threaded(aa)
    ee=Threaded(aa,Threads.nthreads())
    ff=Threaded(aa,10*Threads.nthreads())  
    
    @test size(aa) == (10,2)
    @test eltype(aa) == Float64
    @test size(bb) == (10,2)
    @test eltype(bb) == Float64
    @test size(cc) == (10,2)
    @test eltype(cc) == Float64
    @test size(dd) == (10,2)    

    @test eltype(dd) == Float64
    @test size(ee) == (10,2)
    @test eltype(ee) == Float64
    @test size(ff) == (10,2)
    @test eltype(ff) == Float64 
    
    #rng=Random.default_rng()
    #@test rand(rng,aa,100) isa Array{Float64,2}
    #@test rand(rng,bb,100) isa Array{Float64,2}
    #@test rand(rng,cc,100) isa Array{Float64,2}
    #@test rand(rng,dd,100) isa Array{Float64,2}
    #@test rand(rng,ee,100) isa Array{Float64,2}
    #@test rand(rng,ff,100) isa Array{Float64,2}  
    #
    #@test rand(rng,aa) isa Array{Float64,2}
    #@test rand(rng,bb) isa Array{Float64,2}     
    #@test rand(rng,cc) isa Array{Float64,2}
    #@test rand(rng,dd) isa Array{Float64,2}      
    #@test rand(rng,ee) isa Array{Float64,2}
    #@test rand(rng,ff) isa Array{Float64,2}        
#
    #event=Participants(aa,aa,1,1,6.4,1,0)
#
    #
    #
    #@test rand(rng,event,100) isa Vector{Array{Float64,2}}
    #@test rand(rng,event) isa Array{Float64,2}

    @test rand(aa,100) isa Array{Float64,2}
    @test rand(bb,100) isa Array{Float64,2}
    @test rand(cc,100) isa Array{Float64,2}
    @test rand(dd,100) isa Array{Float64,2}
    @test rand(ee,100) isa Array{Float64,2}
    @test rand(ff,100) isa Array{Float64,2}  
    
    @test rand(aa) isa Array{Float64,2}
    @test rand(bb) isa Array{Float64,2}     
    @test rand(cc) isa Array{Float64,2}
    @test rand(dd) isa Array{Float64,2}      
    @test rand(ee) isa Array{Float64,2}
    @test rand(ff) isa Array{Float64,2}        

    event=Participants(aa,aa,1,1,6.4,1,0)

    
    
    #@test rand(event,100)[4].part1 isa Vector{SVector{2,Float64}}
    #@test rand(event,100)[4].part2 isa Vector{SVector{2,Float64}}
    #@test rand(event).part1 isa Vector{SVector{2,Float64}}
    #@test rand(event).part2 isa Vector{SVector{2,Float64}}
#TODO how to add package to test?

end
