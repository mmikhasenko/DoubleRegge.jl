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

export TwoParticleDiffraction
export setsystem!
export G
include("masses.jl")

export modelS
export cosθ1, ϕTY
include("kinematics.jl")
include("modelS.jl")

export arg, Psi
export LMs
include("LMbasis.jl")

export x_IδI_ϕδϕ_compass_ηπ
include("dataIO.jl")

export broadcast_over_tree
export invmasssq, invariants, λ
include("roottrees.jl")

export recamp, constructamps
export dNdϕ, dNdcosθ
export integrate_dcosθdϕ
include("reconstruction.jl")

export phi_asymmetry, phi_asymmetry_2d
include("observables.jl")

export α_a2, α_f2, α_ℙ
export modelDR
export trajectory
export sixexchages
include("modelD.jl")

export bootstrap_band, slice
include("bootstrap.jl")


end  # module DoubleRegge
