# imports
include("iron_line_profile_approx.jl")
using Measures
using DelimitedFiles
gr()


function iron_line_heatmap_approx(
                            tolerance=1e-12,
                            size_multiplier=1,
                            fov=12,
                            dtmax=10000,
                            nbins=100,
                            vars=LinRange(5.0, 85.0, 5),
                            ϵ_3=0, 
                            α_13=0, 
                            α_22=0,
                            α_52=0
                            )
    # nbins=600
    # vars = LinRange(5.0, 85.0, 400)
    # vars = 5.0:5.0:85.0

    # values to be plotted
    bins_vals = Array{Float64}[]
    x_vals_vals = []

    # generating data for different spins
    for var in vars
        fvar = Float64(var)
        x_vals, bins = iron_line_profile_approx(
                                        tolerance=tolerance,
                                        size_multiplier=size_multiplier,
                                        fov=fov,
                                        dtmax=dtmax,
                                        obs_angle=fvar,
                                        output="data",
                                        normalised=true,
                                        nbins=nbins,
                                        ϵ_3=ϵ_3, 
                                        α_13=α_13, 
                                        α_22=α_22,
                                        α_52=α_52
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
                ylabel="Observation Angle",
                cb_title="Flux (Arbitrary Units)"
                )
end

hmp = iron_line_heatmap_approx()
writedlm("heatmap_data.txt", new)
display(hmp)
png(hmp, "iron_heatmap.png")