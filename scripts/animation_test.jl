include("temperature_render.jl")

angles = 75.0:1.0:85.0
anim = @animate for angle in angles
    fangle = Float64(angle)
    temperature_render(obs_angle=fangle, tolerance = 1e-4)
end
gif(anim, "hello.gif", fps=10)