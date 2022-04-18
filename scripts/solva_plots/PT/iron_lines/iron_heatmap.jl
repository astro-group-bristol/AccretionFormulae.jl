include("../../../good_plots/iron_line_heatmap.jl")

vars = LinRange(5.0, 85.0, 500)

hmp = iron_line_heatmap(tolerance=1e-12,
                        dtmax=0.5,
                        size_multiplier=10,
                        nbins=400,
                        vars=vars)
png(hmp, "iron_heatmap.png")