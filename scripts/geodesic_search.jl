using Revise 
using GeodesicTracer
using ComputedGeodesicEquations
using StaticArrays
using Plots
using AccretionFormulae
using Dates
gr()

init_radii = 10:1.0:100
vϕ_vals = 3.00:0.01:10.00
total = length(init_radii)*length(vϕ_vals) # total number of combinations
found = 0   # number of combinations found

start = now() # timing start
for init_radius in init_radii

    m = BoyerLindquist(M=1.0, a=1.0)
    R_isco = AccretionFormulae.r_isco(m.a, m.M)
    u = @SVector [0.0, init_radius, deg2rad(90.0), 0.0]

    # bypass map_impact_parameters cus it's not ideal for this
    # still have to brute force search this bit
    # or use some analytic result from a paper
    # vs = [@SVector [0.0, 0.0, 0.0, -vϕ/100] for vϕ in 3.96:0.01:3.97]
    vs = [@SVector [0.0, 0.0, 0.0, -vϕ/100] for vϕ in vϕ_vals]
    us = [u for _ in 1:length(vs)]

    sols = tracegeodesics(
        m, 
        us, vs, 
        (0.0, 300.0), 
        μ=1.0;
        abstol=1e-8,
        reltol=1e-8,
        verbose=true
    )
    # @show(typeof(sols))
    for (i, sol) in enumerate(sols)
        # selecting only radius values, and only looking at the end
        new_sols = copy(selectdim(sol, 1, 6))
        min, max =  extrema(new_sols)
        diff = max - min
        
        # if the difference is small enough, and it doesn't crash into the ISCO
        # it is a valid solution, and is displayed
        if (diff <= 1) & (min > R_isco*1.1)
            found +=1
            plt = plot(sol, vars=(8, 6), projection=:polar, range=(0, 21), legend=false)
            title!("Initial Radius=$init_radius, vϕ=$(vϕ_vals[i])")
            display(plt)
        end
    end
end
duration = now() - start # timing finish
println("Done")
println("Found $found solutions in $duration.")
println("Success rate of $(found/total) .")