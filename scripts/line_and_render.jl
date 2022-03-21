include("temperature_render.jl")
include("line_profile_cached.jl")
using Measures
gr()

function combined_plot(;
                        mass = 1,
                        spin = 0.998,
                        obs_angle = 85.0,
                        disc_angle = 90.0,
                        tolerance = 1e-8,
                        dtmax = 1000.0,
                        size_multiplier::Int64 = 5,
                        resolution = 400,
                        fov = 10
                        )
    hmap, cache, title = temperature_render(
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
    profile = histogram(
                        energy, 
                        weights=line_flux, 
                        size=(400,200),
                        xlim=(3,10),
                        ylim=(0,100),
                        legend=false)

    verticle_line_x = fill(6.4, 100)
    verticle_line_y = LinRange(0, 100, 100)
    plot!(
        verticle_line_x, 
        verticle_line_y, 
        linestyle=:dash, 
        linecolor=:black, 
        label=false, 
        # xlims=(4,7.5), 
        # ylims=(0,maximum(last(line_profiles*1.1))),
        lw=2
    )

    # combined plot
    # subplot just for title
    titleplot = plot(title=title, grid = false, showaxis = false, bottom_margin = -50Plots.px, ticks=false)

    l = @layout[A{0.1h}; B{0.6w} C{0.4w}]
    combined_plot = plot(titleplot, hmap, profile, layout=l, size=(1000,500))
end

# combined_plot()

combined_plot(
                mass=10,
                spin=0.998,
                obs_angle=60.0,
                tolerance=1e-12,
                size_multiplier=5,
                fov=10,
                dtmax=20
)

combined_plot(
                mass=10,
                spin=0.998,
                obs_angle=62.5,
                tolerance=1e-12,
                size_multiplier=5,
                fov=10,
                dtmax=20
)

combined_plot(
                mass=10,
                spin=0.998,
                obs_angle=65,
                tolerance=1e-12,
                size_multiplier=5,
                fov=10,
                dtmax=20
)