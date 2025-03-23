module MonteCarloGlauber

# Write your package code here.
using Distributions
using Random 
using SpecialFunctions
using SimpleNonlinearSolve
using StatsBase
using StaticArrays
using LinearAlgebra
using OhMyThreads
using Rotations
using DifferentiationInterface
using QuadGK
using Interpolations
using Base.Threads



include("Metropolis_Hastings.jl")
include("Integrated_wood_sakson_refactor.jl")
include("Integrated_wood_sakson_interp_refactor.jl")
include("Nucleos3D_refactor.jl")
include("participants_distribution_refactor.jl")
include("Nucleos_definition.jl")
include("utils.jl")



export NucleiWoodSaxon3D,IntegratedWoodSaxonInterp,IntegratedWoodSaxon, Participants, Participant, threaded, Threaded, rand 
export impactParameter 

export Uranium, Lead, Gold, Copper, Xenon


end