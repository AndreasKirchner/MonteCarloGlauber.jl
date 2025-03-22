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

#This file has all the functions to run the event by event initial conditions and to generate the fluidum intiail conditions

#using QuadGK
#using Interpolations



#using BenchmarkTools

#using Profile
#using PProf


export NucleiWoodSaxon3D,IntegratedWoodSaxonInterp,IntegratedWoodSaxon,Participants,Participant
threarded, Threarded, rand 





end