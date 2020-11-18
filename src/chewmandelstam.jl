
k(s) = sqrt(s/4-mπ^2)
CM_Pwave(s; R=5) = s/π * quadgk(s′->k(s′)^3/ (1+(R*k(s′))^2) / sqrt(s′) / (s′*(s′-s-1e-7im)), 4mπ^2, Inf)[1]

function BW_Pwave(s; R=Nominal.R, endep=:CM)
    Γ = Γρ;
    if endep == :rho
        p = k(s); p0 = k(mρ^2)
        Γ *= p^3/p0^3*(1+(R*p0)^2) / (1+(R*p)^2) * mρ/sqrt(s)
    end
    if endep == :CM
        r = CM_Pwave(mρ^2; R=R)
        Γ *= (CM_Pwave(s; R=R) - real(r)) / (r - real(r))
    end
    return 1/(mρ^2-s-1im*mρ*Γ)
end

function phiBW_Pwave(s; R=0.5, endep=:CM)
    a = BW_Pwave(s; R=R, endep=endep)
    return atan(imag(a), real(a))
end

I_ωBW(s) = abs2(1/(mω^2-s-1im*mω*Γω))
