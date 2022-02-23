using GeodesicRendering
using AccretionGeometry
using CarterBoyerLindquist
using AccretionFormulae

using Plots
gr()


m = CarterMethodBL(M=1.0, a=1.0)
# observer position
u = [0.0, 1000.0, deg2rad(85.0), 0.0]

# disc has inner radius 10, outer radius 50, perpendicular to the spin axis
d = GeometricThinDisc(1.235, 50.0, deg2rad(90.0))

# wrapper around the redshift function

function redshift(m, sol, max_time; kwargs...)
    u = sol.u[end]
    AccretionFormulae.regular_pdotu_inv(u, sol.prob.p, m)
end

function temperature(m, sol, max_time; kwargs...)
    g = redshift(m, sol, max_time; kwargs...)
    u = sol.u[end]
    M_BH = m.M
    a_star = -m.a
    temperature = AccretionFormulae.temp_obs(u[2], a_star, M_BH, g)
end

# create and compose the ValueFunction
redshift_vf = (
    # calculate redshift
    ValueFunction(temperature)
    # filter only pixels with r > 9, i.e. on the disc
    ∘ FilterValueFunction((m, sol, max_time; kwargs...) -> sol.u[end][2] > 9.0, NaN)
    # filter only pixels that terminated early (i.e., intersected with something)
    ∘ ConstValueFunctions.filter_early_term
)
# do the render

img = @time rendergeodesics(
    m, u, 2000.0,
    d,
    fov_factor=3.0, abstol=1e-9, reltol=1e-9,
    image_width = 175 * 2,
    image_height = 125 * 2,
    vf = redshift_vf
)

# plot
heatmap(reverse(img, dims=1))

