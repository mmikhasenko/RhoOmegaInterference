
B1(z) = z/(1+z)
p(s; sth) = sqrt(s-sth)/2
# 
ρ_2b(s; sth, R) = 2*p(s; sth=sth) / sqrt(s) * B1(p(s; sth=sth)^2*R^2)
ρ2π_2b(s; R) = ρ_2b(s; sth=4mπsq, R=R)
ρ3π_2b(s; R) = ρ_2b(s; sth=0, R=R)
#
λ(x,y,z) = x^2+y^2+z^2-2x*y-2y*z-2z*x
function ρ3π_ρπ(s; R)
    function integrand(σ)
        q = sqrt(λ(s,σ,mπ^2))/(2sqrt(s))
        # ρ → ππ P-wave break up momentum is neglected in both numerator and energy-dep width
        (2q)/sqrt(s) * B1(q^2*R^2) * 2*mρ*Γρ / ((mρ^2-σ)^2+mρ^2*Γρ^2) / (2π)
    end
    return quadgk(integrand, 4mπ^2, (sqrt(s)-mπ)^2)[1]
end
#
