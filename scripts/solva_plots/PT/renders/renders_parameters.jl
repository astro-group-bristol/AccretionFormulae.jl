include("../../../temperature_render.jl")

mass_vars = [10, 100, 10e6]
spin_vars = [0.0, 0.5, 0.998]
angle_vars = [15.0, 40.0, 85.0]

mass_fixed = 10
spin_fixed = 0.998
angle_fixed = 85.0

for mass in mass_vars
    hmap, cache, title, new_img = temperature_render(;
                                                    mass=mass,
                                                    spin=spin_fixed,
                                                    obs_angle=angle_fixed,
                                                    tolerance=1e-12,
                                                    dtmax=0.5,
                                                    size_multiplier=10,
                                                    resolution=3000
                                                    )
    title!(title)
    png(hmap, "render_mass=$mass.png")
end

for spin in spin_vars
    hmap, cache, title, new_img = temperature_render(;
                                                    spin=spin,
                                                    mass=mass_fixed,
                                                    obs_angle=angle_fixed,
                                                    tolerance=1e-12,
                                                    dtmax=0.5,
                                                    size_multiplier=10,
                                                    resolution=3000
                                                    )
    title!(title)
    png(hmap, "render_spin=$spin.png")
end

for angle in angle_vars
    hmap, cache, title, new_img = temperature_render(;
                                                    obs_angle=angle,
                                                    spin=spin_fixed,
                                                    mass=mass_fixed,
                                                    tolerance=1e-12,
                                                    dtmax=0.5,
                                                    size_multiplier=10,
                                                    resolution=3000
                                                    )
    title!(title)
    png(hmap, "render_angle=$angle.png")
end