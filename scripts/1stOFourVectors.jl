using Revise 
using GeodesicTracer
using CarterBoyerLindquist
using StaticArrays
using Plots
using Random
using AccretionFormulae
using ComputedGeodesicEquations
using Dates
gr()

# setup
init_radii = 5:5
min_α = 3.0
max_α = 10.0
stp_α = 1
min_β = 0.0
max_β = 0.0
stp_β = 0.5

# values of α and β to search
α_vals = min_α:stp_α:max_α
β_vals = min_β:stp_β:max_β

total = length(α_vals)*length(β_vals)*length(init_radii) # number of combinations to test
found = 0   # number of solutions found so far
i = 1       # number of combinations tested
start = now() # timing start
for init_radius in init_radii
    for α in α_vals
        for β in β_vals
            println("Testing: $i / $total")

            # path details for the given parameters
            m = BoyerLindquist(M=1.0, a=1.0) #chosen metric
            R_isco = AccretionFormulae.r_isco(m.a, m.M)
            u = @SVector [0.0, init_radius, deg2rad(90.0), 0.0] #initial position (time, radius, theta, phi)
            v = map_impact_parameters(m, u, α, β) #impact parameters α, β to velocity vector at some position u in metric m
            sols = tracegeodesics(m, u, v, (0.0, 10000.0), μ=1.0;
                                    abstol=1e-14,
                                    reltol=1e-14,
                                    dtmax=1.0
                                    )
            
            # selecting only radius values, and only looking at the end
            new_sols = copy(selectdim(sols, 1, 6))
            new_sols = new_sols[Int32(round(length(sols)*0.75)):end]
            # finding the difference between its extremities
            min, max =  extrema(new_sols)
            diff = max - min
            
            # if the difference is small enough, and it doesn't crash into the ISCO
            # it is a valid solution, and is displayed
            if (diff <= 5) & (min > R_isco*1.1)
                found += 1
                plt = plot(sols, 
                            vars=(4+4,2+4), 
                            proj=:polar, 
                            range=(0.0, 30.0), 
                            legend=false, 
                            titlefontsize=8
                            ) #plotting 4 (phi) and 2 (r), in specified radius range
                title!("Initial Radius=$init_radius, α=$α, β=$β; min=$min, max=$max")
                display(plt)
            end
            i+=1
        end
    end
end
duration = now() - start # timing finish
println("Done")
println("Found $found solutions in $duration.")
println("Success rate of $(found/total) .")