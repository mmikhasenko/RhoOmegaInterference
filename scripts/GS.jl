using RhoOmegaInterference
# using ForwardDiff
using Plots


function gs1(s)
    r = GS_disp(mρ^2)
    return (GS_disp(s) - real(r)) / (r - real(r))
end
gs2(s) = GS(s,mρ)/1im


#
let
    plot(layout=grid(1,2), size=(700,400), frame=:origin)
    plot!(sp=1,e->real(gs1(e^2+1e-5im)), -0.1, 1.2, lab="GS DISP")
    plot!(sp=2,e->imag(gs1(e^2+1e-5im)), -0.1, 1.2, lab="GS DISP", ylim=(:auto,:auto))
    plot!(sp=1,e->real(gs2(e^2+1e-5im)), -0.1, 1.2, lab="GS")
    plot!(sp=2,e->imag(gs2(e^2+1e-5im)), -0.1, 1.2, lab="GS", ylim=(:auto,:auto))
end