using Plots
using Measures

function colourbar()
    data = permutedims(LinRange(0,3,1000))
    new = repeat(data,10)

    heatmap(new, 
            xticks=(LinRange(0,1000,4),
            string.(0.0:1:3.0)), 
            yticks=false, 
            size=(400,50), 
            bottom_margin=5mm, 
            left_margin=5mm,
            right_margin=5mm,
            tickfontsize=20,
            cbar=false
            )
end