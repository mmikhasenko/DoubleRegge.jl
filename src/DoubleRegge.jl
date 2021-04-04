module DoubleRegge

using Parameters
using DelimitedFiles
using QuadGK
using PartialWaveFunctions
using GSL
using Statistics
using UpROOT
using StaticArrays
using LinearAlgebra
using Cuba
using SpecialFunctions
using Optim
using Measurements
using TypedTables
using Polynomials

export TwoParticleDiffraction
export setsystem!
export G
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

export phi_asymmetry, phi_asymmetry_2d
include("observables.jl")

export α_a2, α_f2, α_ℙ
export modelDR
export trajectory
export sixexchages
export build_model
include("modelD.jl")

export bootstrap_band, slice
include("bootstrap.jl")

export inlims
export ..
export plotsfolder, fitsfolder
export shiftbyperiod, alignperiodicsequence
include("utils.jl")

export constrained_pw_projection
export constrained_pw_projection_with_derivative
export fold, unfold
include("constrainedpw.jl")

export bartlettambiguities
include("bartlett_ambiguities.jl")

end  # module DoubleRegge
