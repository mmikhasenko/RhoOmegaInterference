using RhoOmegaInterference
using RhoOmegaInterference.GKPY
using Plots
using Plots.PlotMeasures:mm
using LaTeXStrings
theme(:wong2, frame=:box, minorticks=true, grid=false, ylims=(0,:auto), lw=1.5)

list_of_functions = let R0 = 1.45
    [
    L"\mathrm{GKPY\,\,P}\textrm{-}\mathrm{wave}" => e->δ1(e^2)*180/π,
    latexstring("\\mathrm{GS\\,\\,P}\\textrm{-}\\mathrm{wave}") => e->phiBW_Pwave(e^2; R=R0, endep=:GS)*180/π,
    latexstring("\\mathrm{BW\\,\\,P}\\textrm{-}\\mathrm{wave}\\,\\,(R=$(R0)/\\mathrm{GeV})") => 
        e->phiBW_Pwave(e^2; R=R0, endep=:CM)*180/π,
    # latexstring("\\mathrm{CM\\,\\,P}\\textrm{-}\\mathrm{wave}\\,\\,(R=$(R0)/\\mathrm{GeV})") =>
    #     e->phiBW_Pwave(e^2; R=R0, endep=:rho)*180/π,
    L"\mathrm{BW\,\,P}\textrm{-}\mathrm{wave\,\,(R=5/GeV)}" => e->phiBW_Pwave(e^2; R=5)*180/π,
    ]
end;

let ulim = 0.8
    plot(layout=grid(2,1, heights=(0.70,0.20)), link=:x, xlim=(2mπ, ulim), size=(500,500))
    plot!(sp=1, ylab=L"\mathrm{scattering\,\,phase},\,\delta\,\,\mathrm{(deg)}",
        bottom_margin=-3mm, leg=:left)
    for (l,f) in list_of_functions
        plot!(sp=1, f, 2mπ, ulim, lab=l)
    end
    for (l,f) in list_of_functions
        plot!(sp=2, e->f(e)-list_of_functions[1][2](e), 2mπ, ulim, lab="")
    end
    plot!(sp=2, ylim=(-2,2), xlab=L"m_{\pi\pi}\,(\mathrm{GeV})",
        ylab=L"\delta-\delta_{\mathrm{GKPY}}",
        bottom_margin=-15mm)
end
plot!(sp=1, e->δ3(e^2)*180/π, 2mπ, 0.8, lab= L"\mathrm{GKPY\,\,F}\textrm{-}\mathrm{wave}", c=7)

savefig(joinpath("plots","Pwave_phaseshift.pdf"))
