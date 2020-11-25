# nominal model
ρ2π(s) = ρ2π_2b(s; R = Nominal.R)
ρ3π(s) = ρ3π_ρπ(s; R = Nominal.R)
# 
I_ππ2ππ(s) = abs2(T(s+1e-8im; K=K, ρ2π=ρ2π, ρ3π=ρ3π)[1,1])
I_3π23π(s) = abs2(T(s+1e-8im; K=K, ρ2π=ρ2π, ρ3π=ρ3π)[2,2])

# production amplitudes
function A_X2ππ(s;α)
    Tv = T(s+1e-8im; K=K, ρ2π=ρ2π, ρ3π=ρ3π)
    return Tv[1,1] + α*Tv[1,2]/(mρ^2-s)
end
# intensity
I_X2ππ(s; α) = abs2(A_X2ππ(s; α=α))*p(s; sth=4mπ^2)^2
