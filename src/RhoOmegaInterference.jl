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

export k
export CM_Pwave, BW_Pwave, phiBW_Pwave
include("chewmandelstam.jl")

export K, Knoh²
export T, D
include("kmatrix.jl")

include("GKPY.jl")

end