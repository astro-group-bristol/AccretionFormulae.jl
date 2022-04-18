"""
This file defined functions using Optim.jl to automatically discover circular orbits with the `BoyerLindquist`
metric parameters. This method is able to achieve approx 10^-9 numerical accuracy on the true E and Lz of a 
circular orbit.
"""

# using Revise

using GeodesicTracer
using ComputedGeodesicEquations
using GeodesicBase
using CarterBoyerLindquist

using StaticArrays
using Optim

using Plots

# quality of stability function which has a minimum at stable orbits
# this is just a sum of the normalised residuals
Qs(rs) = sqrt(sum((rs ./ rs[1] .- 1.0) .^ 2) / length(rs))

# utility function for tracing a single orbit given some `vϕ`
function trace_single_geodesic(m, u, vϕ)
    v = @SVector [0.0, 0.0, 0.0, vϕ]
    tracegeodesics(
        m,
        u,
        v,
        # we pick a duration that is quite long
        # but you can toy with this
        # the method is pretty good for `r > r_isco` at even just
        # (0.0, 30.0) for the affine time
        (0.0, 300.0),
        μ = 1.0,
        # amazingly, we don't even need _that_ good tolerances
        abstol = 1e-9,
        reltol = 1e-9,
        # important, else it will spam your terminal
        verbose = false,
    )
end

# a little utility method to we don't always have to specify `u` everywhere
function geodesic_for(m, r_init, vϕ)
    u = @SVector [0.0, r_init, deg2rad(90.0), 0.0]
    trace_single_geodesic(m, u, vϕ)
end

# utility function which returns to stability after integrating a geodesic
# for some trial `vϕ`
function estimate_stability(m, u, vϕ)
    sol = trace_single_geodesic(m, u, vϕ)
    rs = selectdim(sol, 1, 6)
    Qs(rs)
end

function find_vϕ_for_orbit(m, r_init; lower_bound = 0.0, upper_bound = 0.1)
    u = @SVector [0.0, r_init, deg2rad(90.0), 0.0]
    # here we use the Optim.jl library
    res = optimize(
        # use an anonymous function to wrap our stability estimator into a univariate function
        vϕ -> estimate_stability(m, u, vϕ),
        lower_bound,
        upper_bound,
        # there are to univariate solvers, the other being `Brent`, but
        # personally I've found more success with `GoldenSection` (https://en.wikipedia.org/wiki/Golden-section_search)
        # even though it is slower
        GoldenSection(),
    )
    # uncomment this to see information about each solve
    # @show res 

    # return the `vϕ` which minimized our function
    Optim.minimizer(res)
end


# this uses a sliding window for the lower and upper bounds on the `vϕ` estimate
# basically, if we overestimate the boundaries, the solution takes too many iterations and
# doesn't converge, and if we miss the estimation, we risk excluding the true solution
#
# the sliding window uses the a priori knowledge that radii closer to the origin will have
# greater `vϕ`, and so adjusts the lower and upper bounds accordingly (though still permits
# slightly lower `vϕ` to avoid possible numerical errors)
#
# it seems to work pretty well, but you may need to adjust the upper bound a little
# the tolerance is good to about a factor of 1000, so if your true velocity is 1 and you set
# the bounds to 1000, you will still find it. setting it higher than this, e.g. 10_000 will
# not converge in the set number of iterations
function find_vϕ_for_orbit_range(m, rs; lower = 0.0, upper = 0.1)
    lower = lower
    upper = upper

    # we solve backwards, from the furthest radii inwards
    # so that we can use a asymptotically flat space guess for the window (which should be metric
    # independent), and roll our way in from there
    map(reverse(rs)) do r
        vϕ = find_vϕ_for_orbit(m, r; lower_bound = lower, upper_bound = upper)
        # continuously adjust based on heuristic
        lower = vϕ * 0.99
        upper = vϕ * 2.0
        vϕ
        # reverse again to make sure the order is correct
    end |> reverse
end


# range is set from approx schwarzschild radius to whatever you like it to be
# maybe it's best just to do from r_isco outwards?
function find_orbit_range(; r_range = 2.0:0.1:10.0, a = -0.4, upper = 0.1)
    m = BoyerLindquist(M = 1.0, a = a)
    vϕs = @time find_vϕ_for_orbit_range(m, r_range; upper = upper)

    # calculate the paths
    geodesics = map(i -> geodesic_for(m, r_range[i], vϕs[i]), eachindex(r_range))

    # assemble a matrix of energies and angular momenta
    # such that `e_lz[i, 1]` gives you the energy of geodesic `i``
    e_lz = [
        f(m, geo.u[1].x[2], geo.u[1].x[1]) for geo in geodesics,
        f in (GeodesicBase.E, GeodesicBase.Lz)
    ]

    # get the time velocity
    vts = [geo.u[1].x[1][1] for geo in geodesics]

    # join time velocities in the last column
    # and angular velocities as the second last column
    hcat(e_lz, vϕs, vts)
end


# test function for a single radius
# incase you want to trial indivial things
function test_single(; r_init = 10.0, a = 0.0, upper = 0.1)

    m = BoyerLindquist(M = 1.0, a = a)
    vϕ = find_vϕ_for_orbit(m, r_init; upper_bound = upper)
    sol = geodesic_for(m, r_init, vϕ)

    plot!(sol, vars = (8, 6), projection = :polar, range = (0.0, r_init), legend = false)
end



# pl = plot()

# m = BoyerLindquist(M=1.0, a=-0.5)
# r_range = 3.0:1.0:13.0
# vϕs = @time find_vϕ_for_orbit_range(m, r_range; upper=0.2)
# geodesics = map(i -> geodesic_for(m, r_range[i], vϕs[i]), eachindex(r_range))

# for sol in geodesics
#     plot!(
#         sol, 
#         vars=(8, 6), 
#         projection=:polar, 
#         range=(0.0, last(r_range)), 
#         legend=false
#     )
# end

# plot!(
#     _ -> GeodesicBase.inner_radius(m),
#     0.0:0.01:2π,
#     c=:black,
#     lw=5
# )
# plot!(
#     _ -> CarterBoyerLindquist.rms(m.M, m.a),
#     0.0:0.01:2π,
#     c=:black,
#     lw=2,
#     ls=:dot
# )

# pl

# savefig("example_circs.svg")
