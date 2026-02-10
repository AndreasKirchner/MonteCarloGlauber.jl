module MonteCarloGlauber

# Write your package code here.
using Distributions
using Random
using SpecialFunctions
using NonlinearSolve
using StatsBase
using StaticArrays
using LinearAlgebra
using OhMyThreads
using Rotations
using DifferentiationInterface
using QuadGK
using Interpolations
using Base.Threads
using DelimitedFiles
using HDF5
using FastGaussQuadrature
using Statistics
using Artifacts


include("Metropolis_Hastings.jl")
include("Integrated_wood_saxson_refactor.jl")
include("Integrated_wood_saxson_interp_refactor.jl")
include("Nucleus3D_refactor.jl")
include("participants_distribution_refactor.jl")
include("Nucleus_definition.jl")
include("utils.jl")
include("Tabulated_Event.jl")


"""
    generate_bg_two_pt_fct(f, delta_factor, norm, Projectile1, Projectile2, w, k, p, sqrtS, bins, mList; 
                           minBiasEvents=1000000, r_grid=0:1:10, step=2pi/20, Threaded=true, n_ext_Grid=0, nFields=10)

Generates the background and two-point correlation function for a given set of parameters.

# Arguments
- `f`: Function to convert initial entropy/enerty to temperature
- `delta_factor`: Function to convert the correltator from entropy/energy to temperature perturbations as function of the background temperature (e.g. dT/ds (T_BG))
- `norm`: Normalization factor for the profiles
- `Projectile1`: The first projectile in the collision system
- `Projectile2`: The second projectile in the collision system
- `w`: Sub-nucleon width parameter
- `k`: Shape parameter for the participant distribution
- `p`: Reduced thickness parameter
- `sqrtS`: Center-of-mass energy of the collision system in GeV
- `bins`: A vector specifying the centrality bins starting implicitly from 0, i.e. [10,20] corresponds to the classes 0-10% and 10-20%
- `mList`: A list of modes for which the two-point correlation function is calculated

# Keyword Arguments
- `minBiasEvents` (default: `1000000`): The minimum number of events to generate for the simulation.
- `r_grid` (default: `0:1:10`): The radial grid used for the trento calculations
- `step` (default: `2pi/20`): The angular step size for the calculations
- `Threaded` (default: `true`): Whether to use multithreading for event generation and processing.
- `n_ext_Grid` (default: `0`): The size of the external grid for the two-point function. The calculation is only carried out for the values specified in `r_grid`. Extra points are padded with zeros.
- `nFields` (default: `10`): The number of fields in the final correlation output.

# Returns
- `bg`: A vector containing the generated background values for each centrality bin.
- `finalCorr`: A 5-dimensional array containing the real part of the two-point correlation function. 
  The dimensions are `[2, nFields, length(mList), nGrid, nGrid]`, where:

# Example
#```julia
bg, finalCorr = generate_bg_two_pt_fct(
    f, delta_factor, norm, Projectile1, Projectile2, 
    w, k, p, sqrtS, bins, mList; 
    minBiasEvents=1000000, r_grid=0:1:10, step=2pi/20, Threaded=true
)

"""
#
#function generate_bg_two_pt_fct
#end

const root_light_ion = artifact"light_ion_conf"


export NucleiWoodSaxon3D, IntegratedWoodSaxonInterp, IntegratedWoodSaxon, Participants, Participant, threaded, Threaded, rand, TabulatedEvent

export impactParameter, multiplicity, center_of_mass, centralities_selection_events, epsilon_n_psi_n

export InverseFunction

export Uranium, Lead, Gold, Copper, Xenon, Oxygen, Neon

export generate_bg, generate_tw_pt_fct_entropy, generate_bg_twpt_fct, change_norm, generate_bg_two_pt_fct_save


end
