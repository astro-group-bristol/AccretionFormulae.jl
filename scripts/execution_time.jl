include("temperature_render.jl")
using Dates
using Plots

size_multipliers = 1:20
times = []

for mult in size_multipliers
    start=now()
    temperature_render(;size_multiplier=mult)
    push!(times, now()-start)
end

plt = plot(size_multipliers, times, xlabel="Size Multiplier", ylabel="Code Execution Time")
display(plt)