using GeodesicTracer
using ComputedGeodesicEquations
using GeodesicBase
using CarterBoyerLindquist

using StaticArrays
using Optim

# include("/Users/lucyackland-snow/Desktop/DevEnv1/dev/AccretionFormulae/src/johannsen.jl")
include("hunting-circular-orbits.jl")
include("alt-emissivity-deviations.jl")

using Plots

# import GeodesicBase, GeodesicTracer
# function GeodesicBase.metric(m::GeodesicTracer.AbstractAutoDiffStaticAxisSymmetricParams{T}, u) where {T}
#     rθ = @SVector [u[2], u[3]]
#     comps = GeodesicTracer.metric_components(m, rθ)
#     @SMatrix [
#     comps[1] 0 0 comps[5]
#     0 comps[2] 0 0
#     0 0 comps[3] 0
#     comps[5] 0 0 comps[4] 
#     ]
# end


a=0.0
M=1.0
# R_isco = CarterBoyerLindquist.rms(M, a)
r_range = 4.0:0.1:8.0
# m = BoyerLindquist(M=M, a=a)

e_lz = find_orbit_range(;m=m, r_range = r_range, a = a, upper = 0.1)

E_vals = []

for i in 1:length(r_range)
        E = e_lz[i, 1]
        push!(E_vals, E)
end

plt = plot(r_range, E_vals, xlabel="r (M)", ylabel="E")
display(plt)

r_isco_num = r_range[findmin(E_vals)[2]]
r_isco_true = AccretionFormulae.r_isco(a, M)

@show(r_isco_num)
@show(r_isco_true)