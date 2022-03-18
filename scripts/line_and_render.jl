include("temperature_render.jl")
include("line_profile_cached.jl")
using Measures
gr()

function combined_plot(
                        mass = 1,
                        spin = 0.998,
                        obs_angle = 85.0,
                        disc_angle = 90.0,
                        tolerance = 1e-8,
                        dtmax = 1000.0,
                        size_multiplier::Int64 = 4,
                        resolution = 400,
                        fov = 10
                        )
    hmap, cache = temperature_render(
                                    mass = mass,
                                    spin = spin,
                                    obs_angle = obs_angle,
                                    disc_angle = disc_angle,
                                    tolerance = tolerance,
                                    dtmax = dtmax,
                                    size_multiplier = size_multiplier,
                                    resolution = resolution,
                                    fov = fov
                                    )
    energy, line_flux = energy_histogram_cache(cache)
    profile = histogram(energy, weights=line_flux, size=(400,200))

    l = @layout[a{0.6w} b{0.4w}]

    combined_plot = plot(hmap, profile, layout=l, size=(1000,500))
end


combined_plot()