include("temperature_render.jl")

angles = 75.0:1.0:85.0
# angles = [-179.0:5.0:-1.0;1.0:5.0:179.0]
n_frames = length(angles)
anim = @animate for (i, angle) in enumerate(angles)
    print("Frame: $i / $n_frames complete")
    fangle = Float64(angle)
    temperature_render(obs_angle=fangle, tolerance = 1e-4)
end
gif(anim, "hello.gif", fps=10)