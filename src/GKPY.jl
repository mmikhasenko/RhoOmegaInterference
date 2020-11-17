module GKPY

# using Plots
using Parameters
using QuadGK

# masses from the paper
const mπ = 0.13957;
const mρ = 0.7736;
const Γρ = 0.146;
const mK = 0.496;

# fit parameters (constrained)
const Ppars = (B0 = 1.043, B1 = 0.19, λ1 = 1.39, λ2 = -1.70, ϵ1 = 0.00, ϵ2 = 0.07, e0 = 1.05)
const Gpars = (B0 = 1.09e5, B1 = 1.41e5, λ = 0.051e5)

export δ1, δ3
export cotδ1, cotδ3

# maps
k(s) = sqrt(s/4-mπ^2)
w(s; s0 = 1.45^2) = (sqrt(s) - sqrt(s0-s)) / (sqrt(s) + sqrt(s0-s))

# phase shifts
# P-wave
function cotδ1_less_1050(s; pars)
    s < 4mπ^2 && return 0.0
    @unpack e0, B0, B1 = pars
    return sqrt(s)/(2*k(s)^3)*(mρ^2-s)*(2mπ^3/(mρ^2*sqrt(s)) + B0 + B1*w(s; s0 = e0^2))
end
function δ1_more_1050(s; pars)
    s < 4mπ^2 && return 0.0
    @unpack  λ1, λ2, ϵ1, ϵ2 = pars
    λ0 = acot(cotδ1_less_1050(4mK^2; pars=pars))
    return λ0 + λ1*(sqrt(s)/(2mK)-1) + λ2*(sqrt(s)/(2mK)-1)^2
end
function δ1(s; pars = Ppars)
    @unpack e0 = pars
    if s < e0^2 
        v = acot(cotδ1_less_1050(s; pars=pars))
        return (v<0) ? v+π : v
    end
    return δ1_more_1050(s; pars=pars)
end
cotδ1(s; pars = Ppars) = cot(δ1(s; pars=pars))

# D-wave
function cotδ3(s; pars=Gpars)
    @unpack B0, B1, λ = pars
    return sqrt(s)/(2*k(s)^7)*mπ^6*(2λ*mπ/sqrt(s) + B0 + B1*w(s))
end
δ3(s; pars=Gpars) = acot(cotδ3(s; pars=pars))

end
