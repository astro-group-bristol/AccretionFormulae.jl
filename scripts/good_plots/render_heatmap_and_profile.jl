include("iron_line_heatmap.jl")
include("../temperature_render.jl")
include("../iron_line_profile.jl")
using Plots
gr()

hmp=iron_line_heatmap()
x_data=0:100
y_data=repeat(x_data, 100)
plot!(x_data, y_data)