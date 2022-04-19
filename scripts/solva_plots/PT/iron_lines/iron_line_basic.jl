include("../../../iron_line_profile.jl")

verticle_line_x = fill(6.4, 100)
verticle_line_y = LinRange(0, 1.05, 100)

mass_vars = [10, 100, 10e6]
spin_vars = [0.0, 0.5, 0.998]
angle_vars = [15.0, 40.0, 85.0]

mass_fixed = 10
spin_fixed = 0.998
angle_fixed = 40.0

# for mass in mass_vars
x_vals, bins = iron_line_profile(;
                                mass=mass_fixed,
                                spin=spin_fixed,
                                obs_angle=angle_fixed,
                                fov=10,
                                tolerance=1e-12,
                                dtmax=0.5,
                                size_multiplier=10,
                                output="data",
                                nbins=200
                                )
plt = plot(
            x_vals, 
            bins,
            xlabel="Energy (keV)",
            ylabel="Flux (Arbitrary Units)",
            legend=false,
            grid=false, 
            framestyle=:box, 
            ylims=(0,1.05)
            )

# adding vertical line at 6.4 keV
plot!(
    verticle_line_x, 
    verticle_line_y, 
    linestyle=:dash, 
    linecolor=:black, 
    label=false,
    ylims=(0,1.05),
    lw=2
    )
png(plt, "iron_line_basic.png")