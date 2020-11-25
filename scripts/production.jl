using RhoOmegaInterference
# 
using Plots
import Plots.PlotMeasures:mm
using LaTeXStrings
using Parameters
# 
theme(:wong2, lw=2, size=(500,350))
const λ = RhoOmegaInterference.λ

gr()

const mDstar0 = 2.00685
const mD0 = 1.86483
const mJψ = 3.0969
pq(s) = sqrt(λ(s,mπ^2,mπ^2)*λ(s,(mDstar0+mD0)^2,mJψ^2))

p1 = let α = 0.15
    plot(sp=1, xlab=L"m_{\pi\pi}\,\,(\mathrm{GeV})",
        title=L"p^3\,q\,|\hat{A}_{\pi\pi}\,|^2", leg=:topleft)
    ev = range(2mπ+1e-9, mDstar0+mD0-mJψ-1e-9, length=100)
    # 
    f1(s) = I_X2ππ(s; α=0.0)
    calv = map(e->pq(e^2)*f1(e^2)/f1(0.73^2), ev)
    plot!(ev, calv, lab=L"\alpha = 0")
    # 
    f2(s) = I_X2ππ(s; α=-α)
    calv = map(e->pq(e^2)*f2(e^2)/f2(0.73^2), ev)
    plot!(ev, calv, lab=latexstring("\\alpha = -$(α)"))
    
    f3(s) = I_X2ππ(s; α=α)
    calv = map(e->pq(e^2)*f3(e^2)/f3(0.73^2), ev)
    plot!(ev, calv, lab=latexstring("\\alpha = $(α)"))
    plot!()
end
savefig(joinpath("plots","X2pipi_intensity.pdf"))

p2 = let α = 0.15
    plot(sp=1, xlab=L"m_{\pi\pi}\,\,(\mathrm{GeV})",
        title=L"p^2\,|\hat{A}_{\pi\pi}\,|^2", leg=:topleft)
    f1(s) = I_X2ππ(s; α=0.0)
    plot!( e->f1(e^2)/f1(0.73^2), 2mπ, 0.8, lab=L"\alpha = 0")
    f2(s) = I_X2ππ(s; α=-α)
    plot!(e->f2(e^2)/f2(0.73^2), 2mπ, 0.8, lab=latexstring("\\alpha = -$(α)"))
    f3(s) = I_X2ππ(s; α=α)
    plot!(e->f3(e^2)/f3(0.73^2), 2mπ, 0.8, lab=latexstring("\\alpha = $(α)"))
    #
    vspan!([mDstar0+mD0-mJψ, 0.8], lab="", l=(:gray, :dash), α=0.1)
end
savefig(joinpath("plots","X2pipi_intensity_nophasespace.pdf"))

plot(p1, p2, layout=grid(1,2), size=(900,350), bottom_margin=3mm)
savefig(joinpath("plots","X2pipi_intensity_combined.pdf"))
