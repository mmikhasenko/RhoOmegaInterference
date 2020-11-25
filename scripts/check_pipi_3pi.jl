using RhoOmegaInterference
# 
using LaTeXStrings
using Parameters
# 
using Plots
import Plots.PlotMeasures:mm
theme(:wong2, lw=2)
#
p1 = let
    plot(xlab=L"m_{\pi\pi}\,\,(\mathrm{GeV})",
        title=L"|\hat{T}(\pi\pi \to \pi\pi)|^2", leg=:topleft)
    f1(s) = I_ππ2ππ(s)
    plot!( e->f1(e^2)/f1(mρ^2-1e-4), 2mπ, 0.8, lab="K-matrix")
    # 
    f2(s) = abs2(BW_Pwave(s; R=Nominal.R, endep=:rho))
    plot!(e->f2(e^2)/f2(mρ^2-1e-4), 2mπ, 0.8, lab="BW")
end
p2 = let
    plot(xlab=L"m_{3\pi}\,\,(\mathrm{GeV})",
        title=L"|\hat{T}(3\pi \to 3\pi)|^2", leg=:topleft)
    plot!(e->I_3π23π(e^2)/I_3π23π(mω^2), mω-2Γω, mω+2Γω, lab="K-matrix")
    plot!(e->I_ωBW(e^2)/I_ωBW(mω^2), mω-2Γω, mω+2Γω, lab="BW")
end

plot(p1,p2, size=(900,350), bottom_margin=2mm)
savefig(joinpath("plots","BW_vs_Kmatrix.pdf"))

# let
#     plot(xlab=L"m_{\pi\pi}\,\,(\mathrm{GeV})",
#         title=L"|A(\rho/\omega \to \pi\pi)|^2", leg=:topleft)
#     f(s) = I_X2ππ(s; α=1000)
#     plot!( e->f(e^2)/f(mρ^2-1e-4), mω-2Γω, mω+2Γω, lab="K-matrix")
#     #
#     f(s) = abs2(BW_Pwave(s; R=Nominal.R, endep=:rho)/(mω^2-s-1im*mω*Γω))
#     plot!(e->f(e^2)/f(mρ^2-1e-4), mω-2Γω, mω+2Γω, lab="BW*BW")
# end
