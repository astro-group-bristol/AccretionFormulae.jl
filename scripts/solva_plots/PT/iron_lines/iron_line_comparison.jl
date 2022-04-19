include("../../../iron_line_profile.jl")

verticle_line_x = fill(6.4, 100)
verticle_line_y = LinRange(0, 1.05, 100)

mass_vars = [10, 100, 10e6]
spin_vars = [0.0, 0.5, 0.998]
angle_vars = [15.0, 40.0, 85.0]

mass_fixed = 10
spin_fixed = 0.998
angle_fixed = 40.0

x_vals_vals = []
bins_vals = []
labels = []
for spin in spin_vars
    x_vals, bins = iron_line_profile(;
                                    spin=spin,
                                    mass=mass_fixed,
                                    obs_angle=angle_fixed,
                                    tolerance=1e-12,
                                    dtmax=0.5,
                                    size_multiplier=10,
                                    resolution=3000,
                                    output="data",
                                    nbins=200,
                                    fov=10
                                    )
    push!(x_vals_vals, x_vals)
    push!(bins_vals, bins)
    push!(labels, "a = $spin")
end
fixed_labels=permutedims(labels)
plt = plot(
            x_vals_vals, 
            bins_vals,
            xlabel="Energy (keV)",
            ylabel="Flux (Arbitrary Units)",
            legend=:topleft,
            label=fixed_labels, 
            grid=false, 
            framestyle=:box,
            ylims=(0,1.05)
            )
plot!(
    verticle_line_x, 
    verticle_line_y, 
    linestyle=:dash, 
    linecolor=:black, 
    label=false,
    ylims=(0,1.05),
    lw=2
    )
png(plt, "iron_line_comp_spin.png")

x_vals_vals = []
bins_vals = []
labels = []
for angle in angle_vars
    x_vals, bins = iron_line_profile(;
                                    obs_angle=angle,
                                    spin=spin_fixed,
                                    mass=mass_fixed,
                                    tolerance=1e-12,
                                    dtmax=0.5,
                                    size_multiplier=10,
                                    resolution=3000,
                                    output="data",
                                    nbins=200,
                                    fov=10,
                                    )
    push!(x_vals_vals, x_vals)
    push!(bins_vals, bins)
    push!(labels, "i = $angle")
end
fixed_labels=permutedims(labels)
plt = plot(
            x_vals_vals, 
            bins_vals,            
            xlabel="Energy (keV)",
            ylabel="Flux (Arbitrary Units)",
            legend=:topleft,
            label=fixed_labels, 
            framestyle=:box, 
            grid=false,
            ylims=(0,1.05)
            )
plot!(
    verticle_line_x, 
    verticle_line_y, 
    linestyle=:dash, 
    linecolor=:black, 
    label=false,
    ylims=(0,1.05),
    lw=2
    )
png(plt, "iron_line_comp_angle.png")