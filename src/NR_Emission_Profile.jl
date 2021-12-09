using Plots

function (D(R))
    D = (3*G*M*mdot)/(8*pi*R^3)*(1-(R_isco/R)^(1/2))
    return D
end

G = 6.67e-11
M = 2e30*10
c = 3e8
eta = 0.1

R_isco = 6*G*M/(c^2)

L_edd = 3e4 *3.8e26 * (2e30*10)/1.99e30

Mdot = L_edd/(c^2*eta)
mdot = 0.1*Mdot

r_vals = LinRange(R_isco,10*R_isco,10000)
D_vals = []

for R in r_vals #creates a series of D (dissipation) values corresponding to radial distances
    push!(D_vals, D(R)) #(cycles through R values, taken from r_vals, 
                        # and inputs them into function D(R), which is a function of R
                        # the resulting value is then added to D_vals)
end

sig_SB = 5.67e-8

T_vals = []
function (T(D))
    T = (D/sig_SB)^(1/4)
    return T
end

for D in D_vals
    push!(T_vals, T(D))

end

plot(r_vals, T_vals)