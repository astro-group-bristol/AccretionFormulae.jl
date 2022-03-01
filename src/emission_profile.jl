#functions
"""
From Reynolds (2020) (Eq2-4)
Calculates the radius of the inner most stable circular orbit of the a black hole,
with spin a_star, and mass M.
"""
function r_isco(a_star, M)
    # M *= 1.99e30
    a = a_star
    r_g = 2*G*M/(c^2)
    # r_g = 1
    z1 = 1+∛(1-a^2)*(∛(1+a)+∛(1-a))
    z2 = √(3*a^2+z1^2)
    if a >= 0
        r_isco = (3+z2-√((3-z1)*(3+z1+2*z2)))*r_g
    elseif a < 0
        r_isco = (3+z2+√((3-z1)*(3+z1+2*z2)))*r_g
    end
end

"""
From Page & Thorne (1974) (Eq15n)
Calculates the function f at radius r, for a black hole of spin a_star, and mass M, 
which is used to calculate the flux given by its accretion disk.
"""
function flux(r, a_star, M)

    R_isco = r_isco(a_star, M)
    x = √(r/M)
    x0 = √(R_isco/M)
    x1 = 2*cos((1/3)*(acos(a_star))-(π/3))
    x2 = 2*cos((1/3)*(acos(a_star))+(π/3))
    x3 = -2*cos((1/3)*(acos(a_star)))

    flux = (3/(2*M))*(1/(x^2*(x^3-(3*x)+(2*a_star))))*(x-x0-((3/2)*a_star*log(x/x0)) 
            - (3*(x1-a_star)^2)/(x1*(x1-x2)*(x1-x3))*log((x-x1)/(x0-x1)) 
            - (3*(x2-a_star)^2)/(x2*(x2-x1)*(x2-x3))*log((x-x2)/(x0-x2)) 
            - (3*(x3-a_star)^2)/(x3*(x3-x1)*(x3-x2))*log((x-x3)/(x0-x3)))
end

function mdot(M)
    L_edd = 3e4*L_☼*(M/M_☼)
    Mdot = -L_edd/(c^2*η)
    mdot = -0.1*Mdot
end

function diss(mdot, r, a_star, M)
    # diss = ((c^6)/(G^2))*mdot*f(r, a_star, M)/(4*π*r)
    diss = mdot*flux(r, a_star, M)/(4*π*r)
end

function temperature(r, a_star, M)
    m_dot = mdot(M)
    temperature = (diss(m_dot, r, a_star, M)/σ_SB)^(1/4)
end

"""
From Fanton et al. (1997) (Eq78)
"""
function observed_temperature(r, a_star, M, g)
    T = temperature(r, a_star, M)
    observed_temperature =  g*T
end

# constants
G = 6.67e-11
c = 3e8
L_☼ = 3.8e26
M_☼ = 1.99e30
σ_SB = 5.67e-8

η = 0.1
# M = 10*M_☼
# a_star = 0.1

export observed_temperature, r_isco