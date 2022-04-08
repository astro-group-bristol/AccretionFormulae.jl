using Plots
using QuadGK
using GeodesicBase

α_13 = 0
α_22 = 0
α_13 = 0
α_52 = 0

ϵ_3 = 0

M = 1.99e30
r = M
vϕ_vals = []
# a = 0.0

# Johannsen 2014 eq 28
function gsqrt(ϵ_3)
    gsqrt =   (r/abs(1 + (α_13 * M^3)/(r^3) -  
                        (α_22 * a^2 * M^2)/(r^4) + 
                        (α_13 * a^2 * M^3)/(r^5))) *
                √(((1 + (ϵ_3 * M^3)/(r^3))^3)/
                (1 + (α_52 * M^2)/(r^2)))
end

u = [0, init_radius, π/2, 0]
vs = [@SVector [0.0, 0.0, 0.0, -vϕ/100] for vϕ in vϕ_vals]

m = CarterMethodBL(M = 1.0, a = 0.998)

for v in vs
    En = GeodesicBase.E(m, u, v)
    L_z = GeodesicBase.Lz(m, u, v)

    values = (r, v[4], En, L_z)
    

Ω = 0.1
L_z = 1

pt1 = 0
pt2 = (M^2)/(gsqrt(ϵ_3)*(E - Ω*L_z))
pt3 = quadgk(x -> ((α*x^3)/(β*x*2 - (α - x)^2)^(1/2)), r_isco, r, rtol=1e-3)

f_disk = pt1 + pt2 + pt3
print(f_disk)