
function K(s)
    g = [gρ, 0]; h = [gωππ, gω]
    Kv = g' .* g ./ (mρ^2-s) + h' .* h ./ (mω^2-s)
end
function Knoh²(s)
    g = [gρ, 0];
    Kv = g' .* g ./ (mρ^2-s) + [0 gωππ*gω; gωππ*gω gω^2] ./ (mω^2-s)
end

function T(s; K=error("K"), ρρ = ρρ_2b, ρω = ρω_2b)
    𝕀 = Matrix{Float64}(I,(2,2))
    Rho = [ρρ(s) 0; 0 ρω(s)]
    Kv = K(s)
    return Kv * inv(𝕀 - 1im * Rho * Kv)
end
function D(s; K=error("K"), ρρ = ρρ_2b, ρω = ρω_2b)
    𝕀 = Matrix{Float64}(I,(2,2))
    Rho = [ρρ_2b(s) 0; 0 ρω(s)]
    Kv = K(s)
    return inv(𝕀 - 1im * Rho * Kv)
end
