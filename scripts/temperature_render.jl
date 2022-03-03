using GeodesicRendering
using AccretionGeometry
using CarterBoyerLindquist
using AccretionFormulae
using StaticArrays

using Plots
gr()

function temperature(m, sol, max_time; kwargs...)
    g =  AccretionFormulae.redshift(m, sol, max_time; kwargs...)

    u = sol.u[end]

    M = m.M*1.99e30
    # M = 5.0
    a = m.a
    a_star = a/M
    temperature = AccretionFormulae.observed_temperature(u[2], a_star, M, g)
end

function temperature_render(mass=1.0, spin=0.998;obs_angle=85.0, disc_angle=90.0, 
                            tolerance=1e-8, dtmax=1000.0, 
                            size_multiplier::Int64=1, resolution=400, fov=3.0
                            )
                            
    m = CarterMethodBL(M=mass, a=spin)

    # observer position
    u = [0.0, 1000.0, deg2rad(obs_angle), 0.0]
    R_isco = AccretionFormulae.r_isco(m.a, m.M)

    # new method using the CarterBoyerLindquist RMS gave very different values
    # and wasn't working:
    # R_isco = CarterBoyerLindquist.rms(m.M, m.a)

    # disc has inner radius 10, outer radius 50, perpendicular to the spin axis
    d = GeometricThinDisc(R_isco, 50.0, deg2rad(disc_angle))


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
        fov_factor=fov*size_multiplier, abstol=tolerance, reltol=tolerance,
        image_width = 350*size_multiplier,
        image_height = 250*size_multiplier,
        vf = redshift_vf,
        dtmax = dtmax
    )

    # plot
    scale = 1e4
    new_img = reverse(img, dims=1)
    new_img ./= scale

    heatmap(new_img, aspect_ratio=1.0, size=(resolution*3/2, resolution))
    title!("Temperature $scale, $(m.M)")
end

temperature_render(obs_angle=85.0, size_multiplier=2)