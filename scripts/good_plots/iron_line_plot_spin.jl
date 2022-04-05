# imports
include("../iron_line_profile.jl")
using Measures
gr()

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
# vars = 5:40:85 # obs angles
# vars=[5,15]
vars = [0.0,0.25,0.50,0.75,0.998] # spins
# vars = [10, 30, 60, 75]

# legend labels, colours and styles
labels = permutedims(label.(vars, "Spin",""))
styles = [:solid :dash :dot :dashdot :dashdotdot]
colors = [:black :blue :red :green :purple]

# values to be plotted
bins_vals = []

# generating data for different spins
for var in vars
    fvar = Float64(var)
    x_vals, bins = iron_line_profile(tolerance=1e-12,
                                            size_multiplier=3,
                                            fov=12,
                                            dtmax=5,
                                            spin=fvar,
                                            output="data")
    push!(bins_vals, bins)
end

# adding a line to represent the lab frame 6.4 keV emission line
verticle_line_x = fill(6.4, 100)
verticle_line_y = LinRange(0, 1.2, 100)

# plotting
plt = plot(
        x_vals,
        bins_vals,
        framestyle=:box,
        grid=false,
        xlabel="Energy (keV)",
        ylabel="Flux (Arbitrary Units)",
        label=labels,
        linestyles=styles,
        linecolor=colors,
        legend=:topleft
        )
# vertical plot
plot!(
    verticle_line_x, 
    verticle_line_y, 
    linestyle=:dash, 
    linecolor=:black, 
    label=false,
    ylims=(0,1.2),
    lw=2
    )

# saving image
png(plt, "iron_line_plot_spin.png")
display(plt)