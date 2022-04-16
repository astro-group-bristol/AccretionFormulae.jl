include("iron_line_profile.jl")
include("arbitrary_f/iron_line_profile_approx.jl")

x_vals_exact, bins_exact = iron_line_profile(
                                            tolerance=1e-12,
                                            # size_multiplier=3,
                                            # fov=12,
                                            # dtmax=5,
                                            obs_angle=30,
                                            output = "data"
                                            )


ϵ_3_vals = 1:1:10

labels = []
x_vals_vals = []
bins_vals = []

for ϵ_3 in ϵ_3_vals
    x_vals_approx, bins_approx = iron_line_profile_approx(
                                                tolerance=1e-12,
                                                # size_multiplier=3,
                                                # fov=12,
                                                # dtmax=5,
                                                obs_angle=30,
                                                output = "data",
                                                ϵ_3=ϵ_3
                                                )
    push!(bins_vals, bins_approx)
    push!(x_vals_vals, x_vals_approx)
    push!(labels, "Johannsen f (ϵ_3 = $ϵ_3)")
end
fixed_labels=permutedims(labels)
plot(x_vals_vals, bins_vals, label=fixed_labels)
plot!(x_vals_exact, bins_exact, label="Page Thorne f",legend=:topleft)
