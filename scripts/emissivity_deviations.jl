using Plots
using QuadGK
using GeodesicBase
using CarterBoyerLindquist
using ComputedGeodesicEquations
using GeodesicTracer
using Statistics
using Interpolations

include("geodesic_search.jl")

α_13 = 0
α_22 = 0
α_13 = 0
α_52 = 0

ϵ_3 = 0

M = 1
R_isco = AccretionFormulae.r_isco(0.998, 1)

# getting solutions from the search
r_vals_long, _, r_dict_vϕ, r_dict_Lz, r_dict_En = search(
                                                        init_radii = R_isco+1:0.05:30,
                                                        vϕ_vals = 0.00:0.001:10.00
                                                        )

# setting up unique solution lists
vϕ_vals = []
Lz_vals = []
En_vals = []
r_vals = []

# sorting dictionaries
r_dict_vϕ = sort(r_dict_vϕ)
r_dict_Lz = sort(r_dict_Lz)
r_dict_En = sort(r_dict_En)

# averaging repeat solutions to find the gradients
for init_radius in keys(r_dict_vϕ)
    if length(r_dict_vϕ[init_radius]) != 0
        push!(vϕ_vals,  mean(r_dict_vϕ[init_radius]))
        push!(Lz_vals,  mean(r_dict_Lz[init_radius]))
        push!(En_vals,  mean(r_dict_En[init_radius]))
        push!(r_vals, init_radius)
    end
end

# finding gradients of relevant values
itp_vϕ = interpolate((r_vals,), vϕ_vals, Gridded(Linear()));
itp_Lz = interpolate((r_vals,), Lz_vals, Gridded(Linear()));

"""
Johannsen 2014 eq 28
"""
function gsqrt(ϵ_3, r, a)
    gsqrt =   (r/abs(1 + (α_13 * M^3)/(r^3) -  
                        (α_22 * a^2 * M^2)/(r^4) + 
                        (α_13 * a^2 * M^3)/(r^5))) *
                √(((1 + (ϵ_3 * M^3)/(r^3))^3)/
                (1 + (α_52 * M^2)/(r^2)))
end

f_vals = []
for (i, init_radius) in enumerate(r_vals)
    # vϕ_vals = r_dict_vϕ[init_radius]
    # u = [0, init_radius, π/2, 0]
    # vs = [@SVector [0.0, 0.0, 0.0, -vϕ/100] for vϕ in r_dict_vϕ[init_radius]]
    # m = BoyerLindquist(M = 1.0, a = 0.998)

    solutions = [] # all different combos of solutions for each r value

    # for v in vs
    #     @show(v)
        # En = GeodesicBase.E(m, u, v)
        # L_z = GeodesicBase.Lz(m, u, v)
        # push!(L_z_vals, L_z)
    En = En_vals[i]
    Lz = Lz_vals[i]
    vϕ = vϕ_vals[i]
    g_factor = gsqrt(0, init_radius, 0.998)
    M = 1
    dΩdr = only(Interpolations.gradient(itp_vϕ, init_radius))
    dLzdr = only(Interpolations.gradient(itp_Lz, init_radius))

    #     values = (init_radius, v[4], En, L_z)
    #     push!(solutions, values)
    # # end
    
    # @show(solutions)


    pt1 = dΩdr
    pt2 = (M^2)/(g_factor*(En - vϕ*Lz)^2)
    pt3 = quadgk(x -> ((En-vϕ*Lz)*dLzdr), init_radius, R_isco, rtol=1e-3)

    @show(pt1)
    @show(pt2)
    @show(pt3)

    f_disk = -(pt1 * pt2 * pt3[1])
    push!(f_vals, f_disk)
end
# f_vals .^= (1/4)
plt = plot(r_vals, f_vals)
@show(plt)
# L_z_vals_mean = []
# for init_radius in keys(r_dict_vϕ)
#     push!(L_z_vals_mean,  mean(r_dict_vϕ[init_radius]))
#     push!(r_vals, init_radius)
# end

# itp_L_z = interpolate((r_vals,), vϕ_vals, Gridded(Linear()));


# plt = plot(r_vals_long, L_z_vals)
# display(plt)