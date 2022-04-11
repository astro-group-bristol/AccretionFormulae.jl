using Revise 
using GeodesicTracer
using ComputedGeodesicEquations
using StaticArrays
using Plots
using AccretionFormulae
using Dates

using GeodesicBase
using CarterBoyerLindquist
gr()

function search(;
                init_radii = 10:0.1:20,
                vϕ_vals = 1.00:0.01:10.00
                )       
    # values to test
    init_radii = init_radii
    vϕ_vals = vϕ_vals

    # dictionaries for the multiple solutions of each result at each radius
    r_dict_vϕ = Dict()
    r_dict_Lz = Dict()
    r_dict_En = Dict()
    for init_radius in init_radii
        r_dict_vϕ[init_radius] = []
        r_dict_Lz[init_radius] = []
        r_dict_En[init_radius] = []
    end

    total = length(init_radii)*length(vϕ_vals) # total number of combinations
    found = 0   # number of combinations found
    r_vals=[]   # found vals (maybe not needed)
    v_vals=[]
    i = 0 # current test

    start = now() # timing start
    for init_radius in init_radii
                    
        println("Testing Radius: $i / $(length(init_radii))")
        m = BoyerLindquist(M=1.0, a=1.0)
        R_isco = AccretionFormulae.r_isco(m.a, m.M)
        u = @SVector [0.0, init_radius, deg2rad(90.0), 0.0]

        # bypass map_impact_parameters cus it's not ideal for this
        # still have to brute force search this bit
        # or use some analytic result from a paper
        # vs = [@SVector [0.0, 0.0, 0.0, -vϕ/100] for vϕ in 3.96:0.01:3.97]
        vs = [@SVector [0.0, 0.0, 0.0, vϕ/100] for vϕ in vϕ_vals]
        us = [u for _ in 1:length(vs)]

        sols = tracegeodesics(
            m, 
            us, vs, 
            (0.0, 8.0*init_radius^(3/2)), 
            μ=1.0;
            abstol=1e-8,
            reltol=1e-8,
            verbose=true
        )

        # @show(typeof(sols))
        for (j, sol) in enumerate(sols)
            # selecting only radius values, and only looking at the end
            new_sols = copy(selectdim(sol, 1, 6))
            min, max =  extrema(new_sols)
            diff = max - min
            
            # if the difference is small enough, and it doesn't crash into the ISCO
            # it is a valid solution, and is displayed
            if (diff <= 0.05) & (min > R_isco*1.1)
                found +=1
                # plt = plot(sol, vars=(8, 6), projection=:polar, range=(0, 21), legend=false)
                # title!("Initial Radius=$init_radius, vϕ=$(vϕ_vals[j])")
                # display(plt)
                vϕ = vϕ_vals[j]
                u = [0, init_radius, π/2, 0]
                v = @SVector [0.0, 0.0, 0.0, -vϕ/100]
                m = BoyerLindquist(M = 1.0, a = 0.998)
                
                En = GeodesicBase.E(m, u, v)
                Lz = GeodesicBase.Lz(m, u, v)

                push!(r_vals, init_radius)
                push!(v_vals, vϕ)
                push!(r_dict_vϕ[init_radius], vϕ)
                push!(r_dict_Lz[init_radius], Lz)
                push!(r_dict_En[init_radius], En)
            end
            
        end
        i += 1
    end
    duration = now() - start # timing finish
    println("Done")
    println("Found $found solutions in $duration.")
    println("Success rate of $(found/total) .")

    # r_vals, v_vals = collect(zip(stable_orbits))
    plt = plot(r_vals, v_vals, xlabel="Initial Radius", ylabel="vϕ")
    display(plt)
return r_vals, v_vals, r_dict_vϕ, r_dict_Lz, r_dict_En
end