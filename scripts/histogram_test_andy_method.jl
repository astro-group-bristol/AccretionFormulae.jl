using GeodesicRendering
using AccretionGeometry
using CarterBoyerLindquist
using AccretionFormulae
using StaticArrays
using Printf

using Plots
gr()

function temperature(m, gp, max_time; kwargs...)
    g = AccretionFormulae.redshift(m, gp, max_time; kwargs...)

    M = m.M
    a_star = m.a
    temperature = AccretionFormulae.observed_temperature(gp.u[2], a_star, M, g)

end

function radius(m, gp, max_time; kwargs...)
    gp.u[2]
end

function temperature_render(;
    mass = 1,
    spin = 0.998,
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
)

    m = CarterMethodBL(M = 1.0, a = spin)
    M = m.M

    # observer position
    u = [0.0, 1000.0, deg2rad(obs_angle), 0.0]
    R_isco = AccretionFormulae.r_isco(m.a, m.M)

    # disc
    d = GeometricThinDisc(R_isco + 1, 50.0, deg2rad(disc_angle))


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

    redshift_vf = (
        PointFunction(AccretionFormulae.redshift) ∘
        FilterPointFunction((m, gp, max_time; kwargs...) -> gp.u[2] > R_isco, NaN) ∘
        ConstPointFunctions.filter_early_term
    )

    radius_img = GeodesicRendering.apply(radius_vf, cache)
    temperature_img = GeodesicRendering.apply(temperature_vf, cache)
    redshift_img = GeodesicRendering.apply(redshift_vf, cache)


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
    # scale = 1e7
    # scalestr = @sprintf "%.E" scale
    # new_img = reverse(temperature_img, dims=1)
    # new_img ./= scale

    # hm = heatmap(new_img, aspect_ratio=1.0, size=(resolution*3/2, resolution), 
    # clim=(0,3)
    # )
    # # contour(new_img, aspect_ratio=1.0, size=(resolution*3/2, resolution), clim=(0,3))
    # title!("Temperature Scale = $scalestr, Mass = $mass M_☼")
    # # display(hm)
    return temperature_img, redshift_img, radius_img, M_phys
end


function energy_histogram(;obs_angle=30.0, spin=0.998, size_multiplier=1, fov=6.0, tolerance=1e-8, dtmax=1000)
    temperature_img, redshift_img, radius_img, M_phys= temperature_render(obs_angle=obs_angle,
                                                                                spin=spin, 
                                                                                size_multiplier=size_multiplier,
                                                                                fov=fov, 
                                                                                tolerance=tolerance,
                                                                                dtmax=dtmax)
    m_dot = AccretionFormulae.mdot(M_phys)
    # Define an energy array and a line profile array
    n_energies = 100
    # Minimum energy in keV
    e_min = 1.0
    # Maximum energy in keV
    e_max = 10.0
    energies = range(e_min, stop = e_max, length = n_energies)
    lineProfile = zeros(n_energies)
    
    # Define inner and outer radius of disk (could do this direclty using GeometricThinDisc)
    r_in = 5.0
    r_out = 10.0

    for (i, x) in enumerate(temperature_img)
        for (j, y) in enumerate(x)
            g = redshift_img[i][j]
            if !isnan(g)
                en = 6.4 * g
                index = trunc(Int, 1 + (en - e_min) * (n_energies - 1) / (e_max - e_min))
                if (index >= 1) && (index <= n_energies)
                    # Example line profile where the emissivity is proportional to radius^-3 (an arbitrary choice)
                    if (radius_img[i][j] > r_in) && (radius_img[i][j] < r_out)
                        lineProfile[index] += AccretionFormulae.diss(m_dot, radius_img[i][j], spin, 1)* g^4
                    end
                end
            end
        end
    end

    # plot line profile
    profile =plot(
        energies,
        lineProfile,
        label = string(obs_angle),
        xlabel = "Energy",
        ylabel = "Line flux (arbitrary units)",
    )
    display(profile)

    # scale = maximum(filter(!isnan,lineProfile))
    # scale = floor(log(10, scale))
    # scale = 10^scale
    # lineProfile /= scale
    # profile = plot(energies, lineProfile, label=string(obs_angle), xlabel="Energy", ylabel="Line flux (arbitrary units)")
    # display(profile)
    # plot r image
    # heatmap(radius_img, aspect_ratio=1.0, clim=(0,50))
    return energies, lineProfile

end

# energy_histogram()
