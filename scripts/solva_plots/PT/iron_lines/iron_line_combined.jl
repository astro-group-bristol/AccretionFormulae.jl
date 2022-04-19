include("../../../good_plots/line_and_render.jl")

angle_vars = [15.0, 40.0, 85.0]

mass_fixed = 10
spin_fixed = 0.998
angle_fixed = 40.0

for angle in angle_vars
    plt = combined_plot(;
                        mass=mass_fixed,
                        spin=spin_fixed,
                        obs_angle=angle,
                        tolerance=1e-12,
                        dtmax=0.5,
                        size_multiplier=10,
                        resolution=3000,
                        nbins=200,
                        fov=10
                        )
    png(plt, "iron_line_combined_angle=$angle.png")
end