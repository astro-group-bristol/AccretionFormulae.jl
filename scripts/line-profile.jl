# i include this just incase
using Revise

using GeodesicRendering
using AccretionGeometry
using CarterBoyerLindquist
using AccretionFormulae
using StaticArrays

using Plots
gr()

function flux(m, gp, max_time)
    g = AccretionFormulae.redshift(m, gp, max_time)
    r = gp.u[2]
    g^4 / r^3
end

function energy_histogram_ferg(;
                            obs_angle=40,
                            spin=0.998,
                            fov=9.0,
                            tolerance=1e-9,
                            size_multiplier=4,
                            dtmax=1000
                        )
    m = CarterMethodBL(M=1.0, a=spin)

    risco = AccretionFormulae.r_isco(m.a, m.M)
    d = GeometricThinDisc(risco, 50.0, deg2rad(90))

    u = @SVector [0.0, 1000.0, deg2rad(obs_angle), 0.0]

    cache = @time prerendergeodesics(
        m,
        u,
        2000.0,
        d,
        fov_factor = fov * size_multiplier,
        abstol = tolerance,
        reltol = tolerance,
        image_width = 350 * size_multiplier,
        image_height = 250 * size_multiplier,
        dtmax = dtmax
    )

    # define a common filter 
    filter = ConstPointFunctions.filter_intersected


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
    return vec(energy), vec(line_flux)

    # his = histogram(vec(energy), weights=vec(line_flux), legend=false)
    # display(his)
    # png(his, "histogram_method.png")
    # yaxis!("Line Flux")
    # xaxis!("Energy")
end


# energy_histogram_ferg()