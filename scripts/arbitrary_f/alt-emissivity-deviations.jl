"""
Uses `hunting-circular-orbits.jl` and the Page & Thorne (1974) eq. (12)
to simulate a flux model.

Works for a < 0.97, then the orbit finder starts to panic close to the ISCO.
"""

# using Revise

using GeodesicTracer
using ComputedGeodesicEquations
using GeodesicBase
using AccretionFormulae

import CarterBoyerLindquist # for rms function
import LinearAlgebra

using StaticArrays
using Plots
using Interpolations
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
"""
Johannsen 2014 eq 28
"""
function sqrt_g_tilde(r; a=0.97, M=0, ϵ_3=0, α_13=0, α_22=0, α_52=0)
    sqrt_g_tilde =   (r/abs(1 + (α_13 * M^3)/(r^3) -  
                        (α_22 * a^2 * M^2)/(r^4) + 
                        (α_13 * a^2 * M^3)/(r^5))) *
                √(((1 + (ϵ_3 * M^3)/(r^3))^3)/
                (1 + (α_52 * M^2)/(r^2)))
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
function f_approx(Lz, ∇rwt, ∇rLz_array, wt_array, r; a=0.97, M=1.0, ϵ_3=0, α_13=0, α_22=0, α_52=0)
    -(∇rwt / Lz) * sum(∇rLz_array ./ wt_array) * r/sqrt_g_tilde(r; a, M, ϵ_3, α_13, α_22, α_52)
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
function f_approx_test(r_range, res; a=0.97, M=1.0, ϵ_3=0, α_13=0, α_22=0, α_52=0)
    wts = res[:, 4]
    Lzs = res[:, 2]

    ∇rwt_array = ∇r(wts, r_range)
    ∇rLz_array = ∇r(Lzs, r_range)

    fs = map(2:length(r_range)) do i
        @views f_approx(Lzs[i], ∇rwt_array[i], ∇rLz_array[1:i], wts[1:i], r_range[i]; a=a, M=M, ϵ_3=ϵ_3, α_13=α_13, α_22=α_22, α_52=α_52)
    end

    # we rescale here for a better comparison
    fs = fs ./ maximum(fs)

    # we lose a point, so select all but last `r_range`
    # since our derivative function needs at least 2 points to work
    # plot(r_range[2:end], fs, c = :black, xlabel = "r", ylabel = "f", label = "simulated f")
end

# this works up until about 0.97 then, the orbit search
# starts to break down

function init(;a=0.97, M=1.0, disc_radius=50.0, ϵ_3=0, α_13=0, α_22=0, α_52=0)
    m = BoyerLindquist(M = M, a = a)
    u = @SVector [0.0, 13.0, deg2rad(90), 0.0]

    # making the steps too big over-estimates the gradient
    r_range = CarterBoyerLindquist.rms(m.M, m.a):0.01:disc_radius
    res = find_orbit_range(; r_range = r_range, a = m.a)

    # creating an interpolated function of f

    fs = f_approx_test(r_range, res; a=a, M=M, ϵ_3=ϵ_3, α_13=α_13, α_22=α_22, α_52=α_52)
    return fs, r_range
end

function f_approx_func(r, fs, r_range; a=0.97, M=1.0, disc_radius=50.0)
    
    r_range_short = collect(r_range)
    pop!(r_range_short)
    interp_linear_f = LinearInterpolation(r_range_short, fs,extrapolation_bc=Line())
    interp_linear_f(r)
end

G = 6.67e-11
c = 3e8
L_☼ = 3.8e26
M_☼ = 1.99e30
σ_SB = 5.67e-8
η = 0.1

function diss_approx(mdot, r, a_star, M, fs, r_range)
    diss_approx = ((c^6) / (G^2)) * mdot * f_approx_func(r, fs, r_range; a=a_star, M=M) / (4 * π * r)
end

function temperature_approx(r, a_star, M, fs, r_range)
    m_dot = AccretionFormulae.mdot(M)
    temperature_approx = (diss_approx(m_dot, r, a_star, M, fs, r_range) / σ_SB)^(1 / 4)
end

function observed_temperature_approx(r, a_star, M, g, fs, r_range)
    T = temperature_approx(r, a_star, M, fs, r_range)
    observed_temperature_approx = g * T
end