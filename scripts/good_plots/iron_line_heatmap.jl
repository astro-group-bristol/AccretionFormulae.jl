# imports
include("../iron_line_profile.jl")
using Measures
using DelimitedFiles
gr()

nbins=600
vars = LinRange(5.0, 85.0, 400)
# vars = 5.0:5.0:85.0

# values to be plotted
bins_vals = Array{Float64}[]
x_vals_vals = []

# x_vals_vals, bins_vals = iron_line_profile.(
#                                             tolerance=1e-8,
#                                             # size_multiplier=3,
#                                             # fov=12,
#                                             # dtmax=5,
#                                             obs_angle=vars,
#                                             output="data",
#                                             normalised=true
#                                             )

# generating data for different spins
for var in vars
    fvar = Float64(var)
    x_vals, bins = iron_line_profile(
                                    tolerance=1e-12,
                                    size_multiplier=8,
                                    fov=12,
                                    dtmax=2,
                                    obs_angle=fvar,
                                    output="data",
                                    normalised=true,
                                    nbins=nbins
                                    )

    push!(bins_vals, bins)
    push!(x_vals_vals, x_vals)
end

new = mapreduce(permutedims, vcat, bins_vals)
hmp = heatmap(
            new,
            xticks=(LinRange(0,nbins,5),string.(LinRange(0,10,5))),
            yticks=(LinRange(0,length(vars),5),string.(round.(vars))),
            xlabel="Energy (keV)",
            ylabel="Observation Angle"
            )
writedlm("heatmap_data.txt", new)
display(hmp)
png(hmp, "iron_heatmap.png")