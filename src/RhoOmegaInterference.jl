module RhoOmegaInterference

using QuadGK
using LinearAlgebra
using Parameters

export mπ, mπsq
export mρ, mω
export Γρ, Γω, Brωππ
include("constants.jl")

export ρ_2b, ρ2π_2b, ρ3π_2b
export ρ3π_ρπ
include("phasespace.jl")

export Nominal
include("settings.jl")

export ρ2π, ρ3π
export A_X2ππ, I_X2ππ
export I_ππ2ππ, I_3π23π
include("intensities.jl")

export k, Φ2
include("kinematics.jl")

export GS
export GS_disp
# 
export CM_Pwave, BW_Pwave, phiBW_Pwave
export I_ωBW
include("chewmandelstam.jl")

export K, Knoh²
export T, D
include("kmatrix.jl")

include("GKPY.jl")

end