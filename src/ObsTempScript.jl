using GeodesicRendering
using AccretionGeometry
using CarterBoyerLindquist
using AccretionFormulae

using Plots
gr()

m = CarterMethodBL(M=1.0, a=0.0)

# observer position
u = [0.0, 1000.0, deg2rad(85.0), 0.0]

M_BH = m.M
a_star = m.a

R_isco = AccretionFormulae.r_isco(a_star, M_BH)

# disc has inner radius 10, outer radius 50, perpendicular to the spin axis
d = GeometricThinDisc(R_isco, 50.0, deg2rad(90.0))

# wrapper around the redshift function
function redshift(m, sol, max_time; kwargs...)
    u = sol.u[end] 
    AccretionFormulae.regular_pdotu_inv(u, sol.prob.p, m)
end

function temperature(m, sol, max_time; kwargs...)
    g = redshift(m, sol, max_time; kwargs...)

    u = sol.u[end]

    temperature = AccretionFormulae.temp_obs(u[2], a_star, M_BH, g)
    # @show(temperature)
end


# create and compose the ValueFunction
redshift_vf = (
    # calculate redshift
    ValueFunction(temperature)
    # filter only pixels with r > 9, i.e. on the disc
    ∘ FilterValueFunction((m, sol, max_time; kwargs...) -> sol.u[end][2] > R_isco, NaN)
    # filter only pixels that terminated early (i.e., intersected with something)
    ∘ ConstValueFunctions.filter_early_term
)

# do the render
img = @time rendergeodesics(
    m, u, 2000.0, 
    d, 
    fov_factor=3.3, abstol=1e-9, reltol=1e-9,
    image_width = 350,
    image_height = 250,
    vf =redshift_vf
)

# plot
heatmap(reverse(img, dims=1))

