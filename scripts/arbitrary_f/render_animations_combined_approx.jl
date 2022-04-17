include("temperature_render_approx.jl")
include("line_and_render_approx.jl")

# choosing the range to animate over
# vars = 75.0:1.0:85.0 # changing angle small
# vars = [-179.0:1.0:-1.0;1.0:1.0:179.0] # changing angle full
vars = -2.0:0.2:2.0                 # smaller range of angles
# vars = 1.0:5.0:200 # changing mass
# vars = 0:0.1:0.998 # changing spin

# generating frames
n_frames = length(vars)
anim = @animate for (i, var) in enumerate(vars)
    print("Frame: $(i-1) / $n_frames complete")
    fvar = Float64(var)
    # temperature_render(
    #                     mass=10,
    #                     spin=var, 
    #                     tolerance=1e-8, 
    #                     resolution=1080, 
    #                     size_multiplier=2
    #                     )
    combined_plot_approx(
                    mass=10,
                    spin=0.96,
                    obs_angle=40,
                    tolerance=1e-12,
                    size_multiplier=1,
                    fov=10,
                    dtmax=20,
                    α_22=fvar
    )
end

# saving animation
gif(anim, "combined_gif_approx.gif", fps=5)
