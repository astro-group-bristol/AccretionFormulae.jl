using Plots
gr()

pl = plot()

m = BoyerLindquist(M=1.0, a=0.998)
CarterBoyerLindquist.rms(m.M, m.a)
# m = JohannsenAD(M=1.0, a=a)
r_range = R_isco:3.0:30
vϕs = @time find_vϕ_for_orbit_range(m, r_range; upper=0.1)
geodesics = map(i -> geodesic_for(m, r_range[i], vϕs[i]), eachindex(r_range))


for sol in geodesics
    plot!(
        sol, 
        vars=(8, 6), 
        projection=:polar, 
        range=(0.0, last(r_range)), 
        legend=false
    )
end

pl