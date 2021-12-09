using Plots

function f(r)

    x = sqrt(r/M)
    x0 = sqrt(R_isco/M)
    x1 = 2*cos((1/3)*(acos(a_star))-(pi/3))
    x2 = 2*cos((1/3)*(acos(a_star))+(pi/3))
    x3 = -2*cos((1/3)*(acos(a_star)))

    flux = (3/(2*M))*(1/(x^(2)*(x^(3)-(3*x)+(2*a_star))))*(x-x0-((3/2)*a_star*log(x/x0)) 
            - (3*(x1-a_star)^2)/(x1*(x1-x2)*(x1-x3))*log((x-x1)/(x0-x1)) 
            - (3*(x2-a_star)^2)/(x2*(x2-x1)*(x2-x3))*log((x-x2)/(x0-x2)) 
            - (3*(x3-a_star)^2)/(x3*(x3-x1)*(x3-x2))*log((x-x3)/(x0-x3)))

    return flux

end

function D(R, f)
    D = mdot*f/(4*pi*R)
    return D
end

G = 6.67e-11
c = 3e8
L_sun = 3.8e26
M_sun = 1.99e30

eta = 0.1
M = 10*M_sun
a_star = 0.5

R_isco = 6*G*M/(c^2)
L_edd = 3e4*L_sun*(M/M_sun)
Mdot = L_edd/(c^2*eta)
mdot = 0.1*Mdot

r_vals = LinRange(R_isco, 10*R_isco, 10000)
D_vals = []

for r in r_vals
    push!(D_vals, D(r, f(r)))
end

plot(r_vals, D_vals)

