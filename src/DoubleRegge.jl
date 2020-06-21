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

export arg, Psi
export LMs
include("LMbasis.jl")

export get_intesity, get_phase
include("dataIO.jl")

export broadcast_over_tree
export invmasssq, invariants, λ
include("roottrees.jl")

export recamp, constructamps
export dNdϕ, dNdcosθ
include("reconstruction.jl")

export modelS
export cosθ1, ϕTY
include("kinematics.jl")
include("masses.jl")
include("modelS.jl")

export α_a2, α_f2, α_ℙ
export modelDR
export trajectory
include("modelD.jl")

export bootstrap_band, slice
include("bootstrap.jl")


end  # module DoubleRegge
