include("../../../good_plots/iron_line_heatmap.jl")

vars = 5.0:85.0

hmp = iron_line_heatmap(
                        tolerance=1e-12,
                        dtmax=5,
                        size_multiplier=4,
                        nbins=200,
                        vars=vars,
                        grid=false
                        )
png(hmp, "iron_heatmap.png")