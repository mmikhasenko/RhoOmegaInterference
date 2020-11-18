using RhoOmegaInterference
using Plots
using LaTeXStrings

theme(:wong, size=(500,350), lw=2)

Nominal.R = 1.5

let
    plot(xlab=L"m_{3\pi}\,(\textrm{GeV})", ylab=L"\textrm{phase-space factor}", leg=:topleft)
    plot!(e->ρ_2b(e^2;sth=(mρ+mπ)^2, R=Nominal.R), mρ+mπ, 3.6, lab=L"\rho\pi\,\,\textrm{two-body}")
    plot!(e->ρ3π_ρπ(e^2; R=Nominal.R), 3mπ, 3.6, lab=L"\rho\pi\,\,\textrm{with}\,\,\rho\to\pi\pi")
end
savefig(joinpath("plots","phase_space_tb_vs_qtb.pdf"))


let xlims = (2mπ,0.8)
    plot(xlab=L"m_{3\pi}\,(\textrm{GeV})", ylab=L"\textrm{phase-space factor}", leg=:topleft)
    f(s) = ρ_2b(s;sth=(3mπ)^2,R=Nominal.R)
    plot!(e->f(e^2)/f(mω^2), xlims..., lab=L"\rho\pi\,\,\textrm{two-body,}\,\,s_\textrm{th}=(3m_\pi)^2")
    f(s) = ρ3π_2b(s; R=Nominal.R)
    plot!(e->f(e^2)/f(mω^2), xlims..., lab=L"\rho\pi\,\,\textrm{two-body,}\,\,s_\textrm{th}=0")
    f(s) = ρ3π_ρπ(s; R=Nominal.R)
    plot!(e->f(e^2)/f(mω^2), xlims..., lab=L"\rho\pi\,\,\textrm{with}\,\,\rho\to\pi\pi")
end
savefig(joinpath("plots","phase_space_tb_vs_qtb.pdf"))



using QuadGK
f(s) = ρ3π_ρπ(s; R=Nominal.R)
g(s) = s/π*quadgk(s′->f(s′)/s′/(s′-s-1e-6im), (3mπ)^2, Inf)[1]
const g0 = g(mω^2)
gsub(s) = g(s)- real(g0)

let
    ev = range(2mπ,2.8,length=50)
    calv = gsub.(ev.^2)
    plot(xlab=L"m_{3\pi}\,(\textrm{GeV})", ylab=L"\textrm{phase-space factor}", leg=:topright,
        title="Chew Mandelstam (dispesion integral)")
    plot!(ev, [real(calv) imag(calv)], lab=["re CM" "im CM"], lw=2)
    lens!([2mπ, 0.8], [-0.09, 0.02], inset = (1, bbox(0.12, 0.0, 0.4, 0.4)), frame=:box)
end
savefig(joinpath("plots","qtb_dispersive_term.pdf"))