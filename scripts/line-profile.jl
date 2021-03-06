# i include this just incase
using Revise

using GeodesicRendering
using AccretionGeometry
using CarterBoyerLindquist
using AccretionFormulae
using StaticArrays

using Plots
gr()

m = CarterMethodBL(M=1.0, a=1.0)

risco = AccretionFormulae.r_isco(m.a, m.M)
d = GeometricThinDisc(risco, 50.0, deg2rad(90))

u = @SVector [0.0, 1000.0, deg2rad(30), 0.0]

cache = @time prerendergeodesics(
    m,
    u,
    2000.0,
    d,
    fov_factor = 9.0 * 4,
    abstol = 1e-9,
    reltol = 1e-9,
    image_width = 350 * 4,
    image_height = 250 * 4
)

# define a common filter 
filter = ConstPointFunctions.filter_intersected

function flux(m, gp, max_time)
    g = AccretionFormulae.redshift(m, gp, max_time)
    r = gp.u[2]
    g^4 / r^3
end

flux_vf = (
    PointFunction(flux) ∘ filter
)

redshift_vf = (
    PointFunction(AccretionFormulae.redshift) ∘ filter
)

flux_img = GeodesicRendering.apply(flux_vf, cache)
redshift_img = GeodesicRendering.apply(redshift_vf, cache)

# select non-NaN
mask = @. ! isnan(redshift_img)
line_flux = flux_img[mask]
energy = (redshift_img .* 6.4)[mask]

# use `vec` to flatten the matrix into a 1d array
his = histogram(vec(energy), weights=vec(line_flux), legend=false)
yaxis!("Line Flux")
xaxis!("Energy")