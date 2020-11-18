using RhoOmegaInterference
# 
using Plots
using LaTeXStrings
using Parameters
# 
theme(:wong)
#
ρ2π(s) = ρ2π_2b(s; R = Nominal.R)
ρ3π(s) = ρ3π_ρπ(s; R = Nominal.R)
# 
I_ππ2ππ(s) = abs2(T(s+1e-8im; K=K, ρ2π=ρ2π, ρ3π=ρ3π)[1,1])
Inoh²_ππ2ππ(s) = abs2(T(s+1e-8im; K=Knoh², ρ2π=ρ2π, ρ3π=ρ3π)[1,1])

let
    plot( sp=1, e->e*I_ππ2ππ(e^2), mω-2Γω, mω+2Γω, lab=L"h^2 \neq 0", xlab=L"m_{\pi\pi}")
    plot!(sp=1, e->e*Inoh²_ππ2ππ(e^2), mω-2Γω, mω+2Γω, lab=L"h^2 = 0")
    vline!(sp=1, [mω], lab=L"m_{\omega}", ls=:dash)
    plot!(sp=1, xlab=L"m_{\pi\pi}\,\,(\mathrm{GeV})",
        title=L"\pi\pi\to \pi\pi\,\,\textrm{cross section}")
    #
    eofΔe(Δe) = (1e-6*Δe+mω)
    ωxzoom = Γω*Brωππ/2*1e6
    plot!(inset=bbox(0.6,0.4,0.35,0.35))
    plot!(sp=2, Δe->eofΔe(Δe)*I_ππ2ππ(eofΔe(Δe)^2), -ωxzoom, +ωxzoom, lab="")
    plot!(sp=2, Δe->eofΔe(Δe)*Inoh²_ππ2ππ(eofΔe(Δe)^2), -ωxzoom, +ωxzoom, lab="")
    vline!(sp=2, [0], lab="", ls=:dash, xlab=L"m_{\pi\pi}-m_\omega\,(\mathrm{keV})", lw=1)
end
savefig(joinpath("plots","pipi_amplitude.pdf"))

# 
# let
#     plot(e->e*abs2(([gωππ gω] * D(e^2+1e-8im;  K=K))[1] / (mω^2-e^2)), 0.6, 1.0, lab="")
#     plot!(e->e*abs2(([gωππ gω] * D(e^2+1e-8im; K=Knoh²))[1] / (mω^2-e^2)), 0.6, 1.0, lab="")
#     vline!([mω], lab="")
# end


# let α = 230
#     P(s) = [gρ 0] ./ (mρ^2-s) + 1.5*cis(α/180*π) * [gωππ gω] ./ (mω^2-s)
#     plot(e->e*abs2(( P(e^2) * D(e^2+1e-8im; K=K))[1]), 0.6, 1.0, lab="")
#     plot!(e->e*abs2((P(e^2) * D(e^2+1e-8im; K=Knoh²))[1]), 0.6, 1.0, lab="")
#     vline!([mω], lab="")
# end
# #
# Tππ(s; k=1.8e-3) = 1 / (mρ^2-s-1im * gρ^2*ρρ(s)) * (1 + k * s / (mω^2-s-1im*mω*Γω))
# # plot(s->abs2(Tππ(s)), 0.5, 0.68, lab="")

# let α = 0, eN = 0.75
#     P(s) = [1; -12*cis(30/180*π)]
#     # f1(e) = 1/e^2 * abs2(( T(e^2+1e-8im; K=K1)*P(e^2))[1]); n1=f1(eN)
#     f2(e) = 1/e^2 * abs2(( T(e^2+1e-8im; K=K2)*P(e^2))[1]); n2=f2(eN)
#     f3(e) = 1/e^2 * abs2(( Tππ(e^2+1e-8im))[1]); n3=f3(eN)
#     # plot( e->f1(e)/n1, √0.5, √0.67, lab="")
#     plot(xlab=L"m(\pi\pi)\,\,(\mathrm{MeV})", title=L"e^+e^-\to\pi^+\pi^-", size=(500,350), color_palette=palette(:wong))
#     plot!(e->f2(e)/n2, √0.5, √0.67,  lab="K / [1-i Rho K] * F")
#     plot!(e->f3(e)/n3, √0.5, √0.67, lab="CLEO fit: BW(1+BW)")
#     vline!([mω], lab="", ann=(mω,0.7,text("m(omega)",7,rotation=90,:bottom)), leg=:topleft)
# end
# savefig(joinpath("plots","comparison_to_FF_complex.pdf"))

# let α = 0, eN = 0.75
#     P(s) = [1; 1.71e-1/(mρ^2-s)]
#     # f1(e) = 1/e^2 * abs2(( T(e^2+1e-8im; K=K1)*P(e^2))[1]); n1=f1(eN)
#     f2(e) = 1/e^2 * abs2(( T(e^2+1e-8im; K=K2)*P(e^2))[1]); n2=f2(eN)
#     f3(e) = 1/e^2 * abs2(( Tππ(e^2+1e-8im))[1]); n3=f3(eN)
#     # plot( e->f1(e)/n1, √0.5, √0.67, lab="")
#     plot(xlab=L"m(\pi\pi)\,\,(\mathrm{MeV})", title=L"e^+e^-\to\pi^+\pi^-", size=(500,350), color_palette=palette(:wong))
#     plot!(e->f2(e)/n2, √0.5, √0.67,  lab="K / [1-i Rho K] * F")
#     plot!(e->f3(e)/n3, √0.5, √0.67, lab="CLEO fit: BW(1+BW)")
#     vline!([mω], lab="", ann=(mω,0.7,text("m(omega)",7,rotation=90,:bottom)), leg=:topleft)
# end
# savefig(joinpath("plots","comarison_to_FF_0.17overpole.pdf"))

# let f2(e) = 1/e^2 * abs2(( T(e^2+1e-8im; K=K2)*P(e^2))[1]); n2=f2(mρ)
#     f3(e) = 1/e^2 * abs2(( Tππ(e^2+1e-8im))[1]); n3=f3(mρ)
# f(p) = sum(x->abs2(f2(x)-f3(x)), )


# let
#     plot(e->e*abs2(T(e^2+1e-8im; K=K1)[1,2]), mω-2Γω*Brωππ, mω+2Γω*Brωππ, lab="")
#     vline!([mω], lab="")
# end




# function K(s; pars)
#     @unpack g1sq, m1sq, g2sq, m2sq, k = pars
#     return [g1sq / (m1sq-s) k; k g2sq / (m2sq-s)]
# end

# function T(s; pars)
#     𝕀 = Matrix{Float64}(I,(2,2))
#     Rho = [ρρ(s) 0; 0 ρω(s)]
#     Kv = K(s; pars=pars)
#     return Kv * inv(𝕀 - 1im * Rho * Kv)
# end

# # function TA(s; pars)
# #     𝕀 = Matrix{Float64}(I,(2,2))
# #     Rho = [ρρ(s) 0; 0 ρω(s)]
# #     #
# #     g = [, ]
# #     Kv = g' .* g ./ (mρ^2-s) + h' .* h ./ (mρ^2-s)
# #     return Kv * inv(𝕀 - 1im * Rho * Kv)
# # end


# using SymPy

# let
#     s, = @vars s positive=true
#     ρ1, ρ2, m1, m2, k = @vars ρ_1 ρ_2 m_1 m_2 k positive=true
#     g1, g2 = @vars g_1 g_2 positive=true
#     #
#     𝕀 = Matrix{Int}(I,(2,2))
#     Rho = [ρ1 0; 0 ρ2]
#     Kv = K(s; pars=(g1sq=g1^2, m1sq=m1^2, g2sq=g2^2, m2sq=m2^2, k=k))
#     v = Kv * inv(𝕀 - 1im * Rho * Kv)
#     v = simplify.(v)
#     v[1,2]
# end

# let
#     s, = @vars s positive=true
#     ρ1, ρ2, m1, m2, k = @vars ρ_1 ρ_2 m_1 m_2 k positive=true
#     g1, g2 = @vars g_1 g_2 positive=true
#     α1ρ, α1ω, α2ρ, α2ω = @vars α1ρ α1ω α2ρ α2ω positive=true
#     #
#     𝕀 = Matrix{Int}(I,(2,2))
#     Rho = [ρ1 0; 0 ρ2]
#     Kv = K(s; pars=(g1sq=g1^2, m1sq=m1^2, g2sq=g2^2, m2sq=m2^2, k=k))
#     #

#     F = [α1ρ*g1/(m1^2-s) + α1ω*g2/(m2^2-s)
#          α2ρ*g1/(m1^2-s) + α2ω*g2/(m2^2-s)]
#     #
#     D = 𝕀 - 1im * Kv * Rho
#     v = inv(D) * det(D) * F
#     v = simplify.(v)
#     # print(v[1], )
#     v[1]
#     # LaTeXString(sympy.latex(v[1]))
# end

# using LaTeXStrings

# plot()
# annotate!([(0.5,0.5,LaTeXString("α"))])

# print
