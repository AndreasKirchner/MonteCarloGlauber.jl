using MonteCarloGlauber
using Test

@testset "MonteCarloGlauber.jl" begin
    aa=NucleiWoodSaxon3D(10, 1.5, 5., 1., 0., 0., 0., 0., 0.)
    bb=IntegratedWoodSaxon(10,1,5.,1.,0.)
    cc=IntegratedWoodSaxonInterp(10,1,5.,1.,0.)
    dd=threarded(aa)
    ee=Threarded(aa,Threads.nthreads())
    ff=Threarded(aa,10*Threads.nthreads())  
    
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
    
    rng=Random.default_rng()
    @test rand(rng,aa,100) isa Array{Float64,2}
    @test rand(rng,bb,100) isa Array{Float64,2}
    @test rand(rng,cc,100) isa Array{Float64,2}
    @test rand(rng,dd,100) isa Array{Float64,2}
    @test rand(rng,ee,100) isa Array{Float64,2}
    @test rand(rng,ff,100) isa Array{Float64,2}  
    
    @test rand(rng,aa) isa Array{Float64,2}
    @test rand(rng,bb) isa Array{Float64,2}     
    @test rand(rng,cc) isa Array{Float64,2}
    @test rand(rng,dd) isa Array{Float64,2}      
    @test rand(rng,ee) isa Array{Float64,2}
    @test rand(rng,ff) isa Array{Float64,2}        

    Parti
end
