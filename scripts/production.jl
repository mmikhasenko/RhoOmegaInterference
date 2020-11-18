using RhoOmegaInterference
# 
using Plots
using LaTeXStrings
using Parameters
# 
theme(:wong2, lw=2, size=(500,350))
const λ = RhoOmegaInterference.λ

const mDstar0 = 2.00685
const mD0 = 1.86483
const mJψ = 3.0969
pq(s) = sqrt(λ(s,mπ^2,mπ^2)*λ(s,(mDstar0+mD0)^2,mJψ^2))

let
    plot(sp=1, xlab=L"m_{\pi\pi}\,\,(\mathrm{GeV})",
        title=L"|X\to \pi\pi\,|^2", leg=:topleft)
    f(s) = I_X2ππ(s; α=0.0)
    plot!( e->pq(e^2)*f(e^2)/f(0.73^2), 2mπ, mDstar0+mD0-mJψ, lab=L"\alpha = 0")
    f(s) = I_X2ππ(s; α=-0.05)
    plot!(e->pq(e^2)*f(e^2)/f(0.73^2), 2mπ, mDstar0+mD0-mJψ, lab=L"\alpha = -0.05")
    f(s) = I_X2ππ(s; α=0.05)
    plot!(e->pq(e^2)*f(e^2)/f(0.73^2), 2mπ, mDstar0+mD0-mJψ, lab=L"\alpha = 0.05")
end
savefig(joinpath("plots","X2pipi_intensity.pdf"))
