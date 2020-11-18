# nominal model
ρ2π(s) = ρ2π_2b(s; R = Nominal.R)
ρ3π(s) = ρ3π_ρπ(s; R = Nominal.R)
# 
I_ππ2ππ(s) = abs2(T(s+1e-8im; K=K, ρ2π=ρ2π, ρ3π=ρ3π)[1,1])
# 
I_3π23π(s) = abs2(T(s+1e-8im; K=K, ρ2π=ρ2π, ρ3π=ρ3π)[2,2])
I_ωBW(s) = abs2(1/(mω^2-s-1im*mω*Γω))

# production amplitudes
function A_X2ππ(s;α)
    Tv = T(s+1e-8im; K=K, ρ2π=ρ2π, ρ3π=ρ3π)
    return Tv[1,1] + α*Tv[1,2]/(mρ^2-s)
end
# intensity
I_X2ππ(s; α) = abs2(A_X2ππ(s; α=α))

# # no h^2 term
# Inoh²_ππ2ππ(s) = abs2(T(s+1e-8im; K=Knoh², ρ2π=ρ2π, ρ3π=ρ3π)[1,1])
# Inoh²_X2ππ(s; α) = abs2(Anoh²_X2ππ(s; α=α))
# Inoh²_3π23π(s) = abs2(T(s+1e-8im; K=Knoh², ρ2π=ρ2π, ρ3π=ρ3π)[2,2])
# function Anoh²_X2ππ(s;α)
#     Tv = T(s+1e-8im; K=Knoh², ρ2π=ρ2π, ρ3π=ρ3π)
#     return Tv[1,1] + α*Tv[1,2]/(mρ^2-s)
# end