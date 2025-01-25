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


include("participants_distribution.jl")
include("nucleos_distribution.jl")
include("FluiduM_running.jl")

#This file has all the functions to run the event by event initial conditions and to generate the fluidum intiail conditions

#using QuadGK
#using Interpolations



#using BenchmarkTools

#using Profile
#using PProf







export trento_cmd_line

end