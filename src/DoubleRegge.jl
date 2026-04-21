module DoubleRegge

using PartialWaveFunctions
using SpecialFunctions
using DelimitedFiles
using LinearAlgebra
using StaticArrays
using Measurements
using TypedTables
using Polynomials
using Statistics
using Parameters
using QuadGK
using UnROOT
using Optim
using Cuba
using GSL

export TwoParticleDiffraction
export ReactionSystem
export compass_ηπ, compass_η′π
include("masses.jl")

export modelS
export cosθ1, ϕTY
include("kinematics.jl")
include("modelS.jl")

export Psi
export pw_project
include("LMbasis.jl")

import Base: rand
export path_ηπ, path_η′π
export description_ηπ, description_η′π
export read_data
export TwoBodyPartialWaves
export TwoBodyPartialWaveAs, TwoBodyPartialWaveIϕs
export changerepresentation
include("dataIO.jl")

export broadcast_over_tree
export invmasssq, invariants, λ, q
include("roottrees.jl")

export recamp, constructamps
export dNdϕ, dNdcosθ
export integrate_dcosθdϕ
export arg, amplitude, Iϕ
include("reconstruction.jl")

export intensity, phase
export phi_asymmetry, phi_asymmetry_2d
include("observables.jl")

export α_a2, α_f2, α_ℙ, α_a2′, α_π1
export modelDR
export trajectory
export Vertex
export ReggeExchange
export EventKinematics
export DoubleReggeModel
export with_parameters
export six_exchanges, ten_exchanges
export v_ππ1η, v_πa2η, v_πa2′η, v_πf2π, v_πℙπ, v_pf2p, v_pℙp
include("modelD.jl")

export bootstrap_band, slice
include("bootstrap.jl")

export inlims
export .., reorder
export plotsfolder, fitsfolder
export shiftbyperiod, alignperiodicsequence
export meanshiftbyperiod
include("utils.jl")

export constrained_pw_projection
export constrained_pw_projection_with_derivative
export fold, unfold
include("constrainedpw.jl")

export bartlettambiguities
include("bartlett_ambiguities.jl")


using Plots.RecipesBase
include("plotting_recipes.jl")

end  # module DoubleRegge
