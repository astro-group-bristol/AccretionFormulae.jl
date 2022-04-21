include("temperature_render.jl")
using BenchmarkTools
using Plots

size_multipliers = 1:20
times = []

for size_multiplier in size_multipliers
    push!(times, @btime temperature_render(;size_multiplier=size_multiplier))
end

plt = plot(size_multiplier, times, xlabel="Size Multiplier", ylabel="Code Execution Time")
display(plt)