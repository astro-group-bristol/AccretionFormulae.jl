include("temperature_render.jl")
include("new_line_profile.jl")
using Measures
gr()

function combined_plot(;
                        mass = 1,
                        spin = 0.998,
                        obs_angle = 85.0,
                        tolerance = 1e-8,
                        dtmax = 1000.0,
                        size_multiplier::Int64 = 1,
                        resolution = 400,
                        fov = 3.0
                        )
    plt, hmap, title, cache = iron_line_profile(
                                    mass = mass,
                                    spin = spin,
                                    obs_angle = obs_angle,
                                    tolerance = tolerance,
                                    dtmax = dtmax,
                                    size_multiplier = size_multiplier,
                                    resolution = resolution,
                                    fov = fov
                                    )

    verticle_line_x = fill(6.4, 100)
    verticle_line_y = LinRange(0, 1, 100)
    plot!(
        plt,
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
    combined_plot = plot(titleplot, hmap, plt, layout=l, size=(1000,500))
end

comb = combined_plot(
                mass=10,
                spin=0.998,
                obs_angle=15.0,
                tolerance=1e-12,
                size_multiplier=1,
                fov=10,
                dtmax=20
)
display(comb)