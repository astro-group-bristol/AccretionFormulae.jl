include("iron_line_profile.jl")
include("arbitrary_f/iron_line_profile_approx.jl")

# x_vals_exact, bins_exact = iron_line_profile(
#                                             tolerance=1e-12,
#                                             # size_multiplier=3,
#                                             # fov=12,
#                                             # dtmax=5,
#                                             obs_angle=30,
#                                             output = "data"
#                                             )


vars = [-2,0,2]

labels = []
x_vals_vals = []
bins_vals = []

for var in vars
    x_vals_approx, bins_approx = iron_line_profile_approx(
                                                tolerance=1e-12,
                                                size_multiplier=10,
                                                fov=12,
                                                dtmax=2,
                                                obs_angle=30,
                                                output = "data",
                                                ϵ_3 = var ,
                                                spin = 0.8,
                                                normalised = false
                                                )
    push!(bins_vals, bins_approx)
    push!(x_vals_vals, x_vals_approx)
    push!(labels, "(ϵ_3 = $var)")
end
fixed_labels=permutedims(labels)
plt1 = plot(x_vals_vals, bins_vals, label=fixed_labels)
display(plt1)
# plot!(x_vals_exact, bins_exact, label="Page Thorne f",legend=:topleft)


labels = []
x_vals_vals = []
bins_vals = []

for var in vars
    x_vals_approx, bins_approx = iron_line_profile_approx(
                                                tolerance=1e-12,
                                                size_multiplier=10,
                                                fov=12,
                                                dtmax=2,
                                                obs_angle=30,
                                                output = "data",
                                                α_22 = var ,
                                                spin = 0.8 ,
                                                )
    push!(bins_vals, bins_approx)
    push!(x_vals_vals, x_vals_approx)
    push!(labels, "(α_22 = $var)")
end
fixed_labels=permutedims(labels)
plt2 = plot(x_vals_vals, bins_vals, label=fixed_labels)
display(plt2)