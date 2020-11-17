module RhoOmegaInterference

using QuadGK
using LinearAlgebra

export mπ, mπsq
export mρ, mω
export Γρ, Γω, Brωππ
include("constants.jl")

export ρ, ρρ_2b, ρω_2b
export gρ, gω, gωππ
include("phasespace.jl")

export k
export CM_Pwave, BW_Pwave, phiBW_Pwave
include("chewmandelstam.jl")

export K, Knoh²
export T, D
include("kmatrix.jl")

include("GKPY.jl")

end