include("histogram_test.jl")

# vars = 75.0:1.0:85.0                 # angle
vars = [-179.0:5.0:-1.0; 1.0:5.0:179.0] # angle
# vars = 1.0:5.0:200                   # mass
n_frames = length(vars)
anim = @animate for (i, var) in enumerate(vars)
    print("Frame: $i / $n_frames complete")
    fvar = Float64(var)
    energy_histogram(obs_angle = fvar)
end
gif(anim, "histgif.gif", fps = 20)
