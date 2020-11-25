
function couplings(; ρ2π, ρ3π)
    gρ = sqrt(mρ*Γρ / ρ2π(mρ^2))
    gω = sqrt(mω*Γω / ρ3π(mω^2))
    gωππ = sqrt(mω*Γω*Brωππ / ρ2π(mω^2))
    return (gρ=gρ, gω=gω, gωππ=gωππ)
end

function K(s; ρ2π, ρ3π)
    @unpack gρ, gω, gωππ = couplings(; ρ2π=ρ2π, ρ3π=ρ3π)
    g = [gρ, 0]; h = [gωππ, gω]
    Kv = g' .* g ./ (mρ^2-s) + h' .* h ./ (mω^2-s)
end

function Knoh²(s; ρ2π, ρ3π)
    @unpack gρ, gω, gωππ = couplings(; ρ2π=ρ2π, ρ3π=ρ3π)
    g = [gρ, 0];
    Kv = g' .* g ./ (mρ^2-s) + [0 gωππ*gω; gωππ*gω gω^2] ./ (mω^2-s)
end

function T(s; K=error("K"), ρ2π, ρ3π)
    𝕀 = Matrix{Float64}(I,(2,2))
    Rho = [ρ2π(s) 0; 0 ρ3π(s)]
    Kv = K(s; ρ2π=ρ2π, ρ3π=ρ3π)
    return Kv * inv(𝕀 - 1im * Rho * Kv)
end

function D(s; K=error("K"), ρ2π, ρ3π)
    𝕀 = Matrix{Float64}(I,(2,2))
    Rho = [ρ2π_2b(s) 0; 0 ρ3π(s)]
    Kv = K(s; ρ2π=ρ2π, ρ3π=ρ3π)
    return inv(𝕀 - 1im * Rho * Kv)
end
