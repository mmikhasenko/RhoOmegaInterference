using RhoOmegaInterference
using RhoOmegaInterference.GKPY
using Plots
using LaTeXStrings

annlab(v) = text(string(round(v,digits=2))*" deg.",:bottom,8)
let ulim = 0.8, R0 = 1.5
    plot(ylab="scattering phase (deg.)", xlab=L"m_{\pi\pi}\,(\mathrm{GeV})",
        color_palette=theme_palette(:wong), leg=:left)
    # 
    v = δ3(ulim^2)*180/π
    plot!(e->δ3(e^2)*180/π, 2mπ, ulim, lab="GKPY F-wave", ann=(ulim-0.02,v,annlab(v)))
    # 
    v = δ1(ulim^2)*180/π
    plot!(e->δ1(e^2)*180/π, 2mπ, ulim, lab="GKPY P-wave", ann=(ulim-0.02,v,annlab(v)))
    #
    plot!(e->phiBW_Pwave(e^2; R=R0, endep=:CM)*180/π, 2mπ, ulim, lab="BW P-wave (R=$(R0)/GeV)",
        ann=(ulim-0.02,v,annlab(v)))
    #
    plot!(e->phiBW_Pwave(e^2; R=R0, endep=:rho)*180/π, 2mπ, ulim, lab="BW P-wave (R=$(R0)/GeV)",
        ann=(ulim-0.02,v,annlab(v)))
    #
    plot!(e->phiBW_Pwave(e^2; R=5)*180/π, 2mπ, ulim, lab="BW P-wave (R=5/GeV)",
        ann=(ulim-0.02,v,annlab(v)))
end
savefig(joinpath("plots","Pwave_phaseshift.pdf"))
# (sin(113.1/180*π)/sin(0.02/180*π))^2
