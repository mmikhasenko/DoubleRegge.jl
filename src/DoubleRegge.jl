module DoubleRegge

using Parameters
using DelimitedFiles
using QuadGK
using PartialWaveFunctions
using GSL
using Statistics

export arg, Psi
export LMs
include("LMbasis.jl")

export get_intesity, get_phase
include("dataIO.jl")

export broadcast_over_tree
export invmasssq, invariants, λ
export cosθη_of_s1t1, ϕTY
include("roottrees.jl")

export recamp, constructamps
export dNdϕ, dNdcosθ
include("reconstruction.jl")

export modelS
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
