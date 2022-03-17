# imports
include("histogram_test_andy_method.jl")
using Measures

# initialising plot backend
multiplier = 1
font_size = multiplier*10
line_width = multiplier*1
size = multiplier.*(600,400)
margin = multiplier*10mm
gr(
    legendfontsize=font_size, 
    margin=margin, 
    right_margin=margin*0.5,
    framestyle=:box, 
    grid=false, 
    lw=line_width, 
    xguidefontsize=font_size,
    yguidefontsize=font_size,
    size=size,
    tickfontsize=font_size,
    thickness_scaling=1
    )

"""
Automatically determines the labels for various lines to go into the legend
"""
function label(value, name, units)
    scale = floor(log(10, value)) - 1
    if scale > 4
        scale = Int32(scale)
        value = "\$10^{$scale}\$"
    else
        value = "\$$value\$"
    end
    label = "$name = $value $(units) "
end

# values to change
vars = 5:20:85 # obs angles
# vars = [0.0,0.25,0.50,0.75,0.998] # spins

# values to be plotted
energy_vals = []
line_profiles = []

# generating data for different spins
for var in vars
    energies, lineProfile = energy_histogram(
                                            spin=0.998,
                                            # spin=var,
                                            # obs_angle=40.0,
                                            obs_angle=var,
                                            size_multiplier=3, 
                                            fov=10, 
                                            # tolerance=1e-9, 
                                            # dtmax=50
                                            )
    # scale = maximum(lineProfile)
    # lineProfile = lineProfile ./ scale
    push!(energy_vals, energies)
    push!(line_profiles, lineProfile)
end

# legend labels, colours and styles
labels = permutedims(label.(vars, "OA",""))
styles = [:solid :dash :dot :dashdot :dashdotdot]
colors = [:black :blue :red :green :purple]

# adding a line to represent the lab frame 6.4 keV emission line
verticle_line_x = fill(6.4, 100)
verticle_line_y = LinRange(-1, maximum(last(line_profiles*1.2)), 100)

# plots
iron_line_plot = plot(
            last(energy_vals), 
            line_profiles, 
            linestyle=styles, 
            # xlims=(4,7.5), 
            # ylims=(0,maximum(last(line_profiles*1.1))), 
            label=labels,
            linecolor=colors,
            legend=:bottomleft,
            xlabel="Energy (keV)",
            ylabel="Flux (Arbitrary Units)"
            )

# vertical plot
plot!(
    verticle_line_x, 
    verticle_line_y, 
    linestyle=:dash, 
    linecolor=:black, 
    label=false, 
    # xlims=(4,7.5), 
    # ylims=(0,maximum(last(line_profiles*1.1))),
    lw=2*multiplier
    )

# saving image
png(iron_line_plot, "iron_line_plot.png")
display(iron_line_plot)