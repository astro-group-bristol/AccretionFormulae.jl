using Revise 
using GeodesicTracer
using CarterBoyerLindquist
using StaticArrays
using Plots
using Random
using Distributions
gr()


for i in 0:10000
    print(i)
    α = rand(Uniform(0,20))
    # @show(α)
    β = 0.0

    m = CarterMethodBL(M=1.0, a=1.0) #chosen metric
    u = @SVector [0.0, 20.0, deg2rad(90.0), 0.0] #initial position (time, radius, theta, phi)
    # v = @SVector [0.0, -1.0, 0.0, sin(deg2rad(0.0))] #initial velocity (time, -1/+1= in/out, theta, phi)
    v = map_impact_parameters(m, u, α, β) #impact parameters α, β to velocity vector at some position u in metric m
    sols = tracegeodesics(m, u, v, (0.0, 100.0); #itime period 0-100
        abstol=1e-9,
        reltol=1e-9
    )
    # @show(unique(sols))

    w = CarterBoyerLindquist.carter_velocity(u, m.E, m.M, m.a, sols.prob.p) #gives Carter 4-velocity
    # print(w)

    new_sols = []
    for val in sols
        new_val = string(val)
        push!(new_sols, new_val)
    end

    if length(unique(new_sols)) != length(sols)
        @show(α)
        @show(β)
        plt = plot(sols, vars=(4,2), proj=:polar, range=(0.0, 30.0), legend=false) #plotting 4 (phi) and 2 (r), in specified radius range
        display(plt)
    end
end

print("Done")