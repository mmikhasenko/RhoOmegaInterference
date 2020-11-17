using QuadGK
using Plots

const mπ = 0.140;
# B terms
ρ(s) = sqrt(1-4mπ^2/s)/(8π)
BW(s; g=1.8, m=0.77) = g^2/(m^2-s-1im*g^2*ρ(s))

compact(s; B) = B(s) + BW(s)/π*quadgk(s′->B(s′)*ρ(s′)/(s′-s-1e-7im), 4mπ^2, Inf)[1]

let
    plot(xlab="m(pipi) (GeV)", color_palette=palette(:wong))
    plot!(e->e*abs2(BW(e^2))*ρ(e^2), 2mπ, 1.1, lab="cT")
    plot!(e->e*abs2(compact(e^2; B=s->1/(1+s)))*ρ(e^2)*1000, 2mπ, 1.1, lab="B+TdispB")
    plot!(e->e*abs2(compact(e^2; B=s->1/(1+s)^2))*ρ(e^2)*1600, 2mπ, 1.1, lab="B'+TdispB'")
end
