include("../../arbitrary_f/temperature_render_approx.jl")

mass_vars = [10, 100, 10e6]
spin_vars = [0.0, 0.5, 0.998]
angle_vars = [15.0, 40.0, 85.0]
ϵ_3_vars = [-2, 2]

mass_fixed = 10
spin_fixed = 0.97
angle_fixed = 85.0

for ϵ_3 in ϵ_3_vars
    hmap, cache, title, new_img = temperature_render_approx(;
                                                    mass=mass_fixed,
                                                    spin=spin_fixed,
                                                    obs_angle=angle_fixed,
                                                    tolerance=1e-12,
                                                    dtmax=0.5,
                                                    size_multiplier=4,
                                                    resolution=3000,
                                                    ϵ_3=ϵ_3
                                                    )
    title!(title)
    png(hmap, "render_ϵ_3_approx=$ϵ_3.png")
end

for angle in angle_vars
    hmap, cache, title, new_img = temperature_render_approx(;
                                                    obs_angle=angle,
                                                    spin=spin_fixed,
                                                    mass=mass_fixed,
                                                    tolerance=1e-12,
                                                    dtmax=0.5,
                                                    size_multiplier=4,
                                                    resolution=3000
                                                    )
    title!(title)
    png(hmap, "render_angle_approx=$angle.png")
end