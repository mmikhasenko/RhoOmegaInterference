
CM_Pwave(s; R=5) = s/π * quadgk(s′->k(s′)^3/ (1+(R*k(s′))^2) / sqrt(s′) / (s′*(s′-s-1e-7im)), 4mπ^2, Inf)[1]

"""
    GS_disp(s)

P-wave factor k^3/√s is used under the dispersion integral.
Two subtractions. The limear term is fixed by condition
that the first derivative is zero at the mass of rho.
"""
GS_disp(s) = 0.061s+s^2/π * quadgk(s′->k(s′)^3 / sqrt(s′) / (s′^2*(s′-s-1e-7im)), 4mπ^2, Inf)[1]

# GS(s) = 
hGS(s) = 2/π*k(s)/√s * log( (√s + 2k(s)) / (2mπ) )
hGS′(s) = (ϵ=1e-7; (hGS(s+ϵ)-hGS(s-ϵ)) / (2ϵ))
# 
realGS(s,m) = - m/k(m^2)^3*(k(s)^2*(hGS(s)-hGS(m^2)) - (s - m^2)*k(m^2)^2*hGS′(m^2) )
imagGS(s,m) = (k(s)/k(m^2))^3 * m / √s
GS(s,m) = realGS(s,m) + 1im*imagGS(s,m)


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
    if endep == :GS
        # r = GS_disp(mρ^2)
        # Γ *= (GS_disp(s) - real(r)) / (r - real(r))
        Γ *= GS(s,mρ)/1im
    end
    return 1/(mρ^2-s-1im*mρ*Γ)
end

function phiBW_Pwave(s; R=0.5, endep=:CM)
    a = BW_Pwave(s; R=R, endep=endep)
    return atan(imag(a), real(a))
end

I_ωBW(s) = abs2(1/(mω^2-s-1im*mω*Γω))
