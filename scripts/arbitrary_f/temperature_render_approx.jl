using GeodesicRendering
using AccretionGeometry
using CarterBoyerLindquist
using AccretionFormulae
using StaticArrays
using Printf
include("alt-emissivity-deviations.jl")

using Plots
gr()

function radius(m, gp, max_time; kwargs...)
    gp.u[2]
end

function temperature_render_approx(;
    mass = 1,
    spin = 0.97,
    disc_radius=50.0,
    obs_angle = 85.0,
    disc_angle = 90.0,
    tolerance = 1e-8,
    dtmax = 1000.0,
    size_multiplier::Int64 = 1,
    resolution = 400,
    fov = 3.0,
    η = 0.1,
    η_phys = 0.1,
    edd_ratio = 0.1,
    edd_ratio_phys = 0.1,
    ϵ_3=0, 
    α_13=0, 
    α_22=0,
    α_52=0
)
    println("Approximating f")
    fs, r_range = init(;a=spin, M=1.0, disc_radius=disc_radius, ϵ_3=ϵ_3, α_13=α_13, α_22=α_22, α_52=α_52)
    println("Approximation Finished")
    m = CarterMethodBL(M = 1.0, a = spin)
    M = m.M

    function temperature(m, gp, max_time; kwargs...)
        g = AccretionFormulae.redshift(m, gp, max_time; kwargs...)
    
        M = m.M
        a_star = m.a
        temperature = observed_temperature_approx(gp.u[2], a_star, M, g, fs, r_range
         )
    
    end

    # observer position
    u = [0.0, 1000.0, deg2rad(obs_angle), 0.0]
    R_isco = AccretionFormulae.r_isco(m.a, m.M)

    # disc
    d = GeometricThinDisc(R_isco + 1, disc_radius, deg2rad(disc_angle))


    # cache the render
    cache = @time prerendergeodesics(
        m,
        u,
        2000.0,
        d,
        fov_factor = fov * size_multiplier,
        abstol = tolerance,
        reltol = tolerance,
        image_width = 350 * size_multiplier,
        image_height = 350 * size_multiplier,
        dtmax = dtmax,
    )

    # create and compose the PointFunctions
    temperature_vf = (
        PointFunction(temperature) ∘
        FilterPointFunction((m, gp, max_time; kwargs...) -> gp.u[2] > R_isco, NaN) ∘
        ConstPointFunctions.filter_early_term
    )

    radius_vf = (
        PointFunction(radius) ∘
        FilterPointFunction((m, gp, max_time; kwargs...) -> gp.u[2] > R_isco, NaN) ∘
        ConstPointFunctions.filter_early_term
    )

    radius_img = GeodesicRendering.apply(radius_vf, cache)
    temperature_img = GeodesicRendering.apply(temperature_vf, cache)

    # correcting for physical mass
    M_phys = mass * 1.99e30
    r_isco_phys = AccretionFormulae.r_isco(M_phys, 0.998)
    r_g_phys = M_phys

    numerators =
        AccretionFormulae.mass_scale_fraction.(
            M_phys,
            η_phys,
            edd_ratio_phys,
            r_isco_phys,
            r_g_phys,
            radius_img * M_phys,
        )
    denominators =
        AccretionFormulae.mass_scale_fraction.(M, η, edd_ratio, R_isco, M, radius_img)

    fractions = numerators ./ denominators
    temperature_img .*= fractions

    # # autoscale
    # scale = maximum(filter(!isnan,temperature_img))
    # scale = floor(log(10, scale))
    # scale = 10^scale

    # scaling image
    scale = 1e7
    scalestr = @sprintf "%.E" scale
    exponent =  Int32(floor(log(10, scale)))
    new_img = reverse(temperature_img, dims = 1)
    new_img ./= scale

    hmap = heatmap(
        new_img,
        aspect_ratio = 1.0,
        size = (resolution * 3 / 2, resolution),
        clim = (0, 3),
    )
    # contour(new_img, aspect_ratio=1.0, size=(resolution*3/2, resolution), clim=(0,3))
    # title!("Temperature Scale = $scalestr, Mass = $mass M_☼, Obs Angle = $obs_angle")
    # title = "Temperature Scale = $scalestr K, Mass = $mass M_☼, Obs Angle = $obs_angle"
    title = "Temperature Scale = \$10^{$exponent}\$ K, Mass = $mass \$\\mathrm{M}_{\\odot}\$, Obs Angle = $obs_angle\$^{\\circ}\$"
    # "")
    return hmap, cache, title, fs, r_range
end

# α_13_vals = LinRange(0,10,3)

# for α_13 in α_13_vals
# hmap, cache, title = temperature_render_approx(obs_angle = 85, mass = 10, resolution=1080, spin=0.97, α_22=α_13)
# # # hmap, cache, title = temperature_render(
# # #                                         mass=10,
# # #                                         spin=0.998,
# # #                                         obs_angle=0.01,
# # #                                         tolerance=1e-12,
# # #                                         size_multiplier=6,
# # #                                         dtmax=1,
# # #                                         resolution=2000
# # #                                         )
# # title!(title)
# display(hmap)
# end