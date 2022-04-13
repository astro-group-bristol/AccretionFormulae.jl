"""
Uses `hunting-circular-orbits.jl` and the Page & Thorne (1974) eq. (12)
to simulate a flux model.

Works for a < 0.97, then the orbit finder starts to panic close to the ISCO.
"""

# using Revise

using GeodesicTracer
using ComputedGeodesicEquations
using GeodesicBase

import CarterBoyerLindquist # for rms function
import LinearAlgebra

using StaticArrays
using Plots
gr()

# hoist other functions
include("hunting-circular-orbits.jl")

# for Kerr, this is just `u[2]`, but to illustrate
# how you could do this for other metrics
function sqrt_g_tilde(m, u)
    metric = GeodesicBase.metric(m, u)
    metric_no_θ = metric[setdiff(1:4, 3), setdiff(1:4, 3)]

    sqrt(-LinearAlgebra.det(metric_no_θ))
end

# gradient function
#
# we just compute discrete derivatives, since this is
# sufficient when we have small `r`.
# using interpolations can make the resolution better
# and increase the accuracy, but this is a good quick-and-dirty way
# to see if our results are going in the right direction
function ∇r(array, r)
    deriv = diff(array) ./ diff(r)
    # ensure same dimension using a lazy extension
    push!(deriv, deriv[end])
    deriv
end

"""
    f_approx(Lz, ∇rwt, ∇rLz_array, wt_array)

A discretized version of `f` from Page & Thorne (1974) eq. (12):

```math
f(r_n) \\approx - \\left. \\left[
   \\frac{1}{L_z} \\left(\\frac{\\partial \\dot{u}^t}{\\partial r} \\right)
\\right] \\right\\rvert_{r = r_n} 
\\sum_{i=0}^{n} \\left\\{ \\left. \\left[
   \\frac{1}{\\dot{u}^t} \\left( \\frac{\\partial L_z}{\\partial r}\\right)
\\right] \\right\\rvert_{r = r_i}
\\right\\}
```

"""
function f_approx(Lz, ∇rwt, ∇rLz_array, wt_array)
    -(∇rwt / Lz) * sum(∇rLz_array ./ wt_array)
end

# use this to see where things break down / if they break down
function plot_derivative_test(r_range, res)
    lzs = res[:, 2]
    # naive derivatives
    derivs = ∇r(lzs, r_range)

    plot(r_range, lzs, label = "Lz", c = :black, legend = false)
    plot!(r_range, derivs, label = "derivative", c = :black, style = :dot, lw = 1.2)
end


# test plot
function f_approx_test(r_range, res)
    wts = res[:, 4]
    Lzs = res[:, 2]

    ∇rwt_array = ∇r(wts, r_range)
    ∇rLz_array = ∇r(Lzs, r_range)

    fs = map(2:length(r_range)) do i
        @views f_approx(Lzs[i], ∇rwt_array[i], ∇rLz_array[1:i], wts[1:i])
    end

    # we rescale here for a better comparison
    fs = fs ./ maximum(fs)

    # we lose a point, so select all but last `r_range`
    # since our derivative function needs at least 2 points to work
    plot(r_range[2:end], fs, c = :black, xlabel = "r", ylabel = "f", label = "simulated f")
end

# this works up until about 0.97 then, the orbit search
# starts to break down
m = BoyerLindquist(M = 1.0, a = 0.97)
u = @SVector [0.0, 13.0, deg2rad(90), 0.0]

println("Starting approximation...")

# making the steps too big over-estimates the gradient
r_range = CarterBoyerLindquist.rms(m.M, m.a):0.01:30.0
res = find_orbit_range(; r_range = r_range, a = m.a)

# use this to check if the derivatives are okay
# plot_derivative_test(r_range, res)

f_approx_test(r_range, res)
println("Approximation done!")

# compare to "true" value
import AccretionFormulae
f_true = AccretionFormulae.flux.(r_range, m.a, m.M)
# rescale as well
f_true = f_true ./ maximum(f_true)
plot!(r_range, f_true, label = "true f")

# savefig("simulated-f-not-log.svg")
