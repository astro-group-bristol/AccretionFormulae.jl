using GeodesicRendering
using AccretionGeometry
using CarterBoyerLindquist
using AccretionFormulae
using StaticArrays
using Printf

using Plots
gr()

function temperature(m, sol, max_time; kwargs...)
    g =  AccretionFormulae.redshift(m, sol, max_time; kwargs...)

    u = sol.u[end]
    M = m.M
    a_star = m.a
    temperature = AccretionFormulae.observed_temperature(u[2], a_star, M, g)
    
end

function radius(m, sol, max_time; kwargs...)
    u = sol.u[end]
    radius = u[2]
end

function temperature_render(;mass=1, spin=0.998, obs_angle=85.0, disc_angle=90.0, 
                            tolerance=1e-8, dtmax=1000.0, 
                            size_multiplier::Int64=1, resolution=400, fov=3.0,
                            η=0.1, η_phys=0.1, edd_ratio=0.1, edd_ratio_phys=0.1
                            )
                            
    m = CarterMethodBL(M=1.0, a=spin)
    M = m.M

    # observer position
    u = [0.0, 1000.0, deg2rad(obs_angle), 0.0]
    R_isco = AccretionFormulae.r_isco(m.a, m.M)

    # disc
    d = GeometricThinDisc(R_isco+1, 50.0, deg2rad(disc_angle))

    # create and compose the ValueFunctions
    temperature_vf = (
        ValueFunction(temperature)
        ∘ FilterValueFunction((m, sol, max_time; kwargs...) -> sol.u[end][2] > R_isco, NaN)
        ∘ ConstValueFunctions.filter_early_term
    )

    radius_vf = (
        ValueFunction(radius)
        ∘ FilterValueFunction((m, sol, max_time; kwargs...) -> sol.u[end][2] > R_isco, NaN)
        ∘ ConstValueFunctions.filter_early_term
    )

    redshift_vf = (
        ValueFunction(AccretionFormulae.redshift)
        ∘ FilterValueFunction((m, sol, max_time; kwargs...) -> sol.u[end][2] > R_isco, NaN)
        ∘ ConstValueFunctions.filter_early_term   
    )

    # do the render
    temperature_img = @time rendergeodesics(
        m, u, 2000.0, 
        d, 
        fov_factor=fov*size_multiplier, abstol=tolerance, reltol=tolerance,
        image_width = 350*size_multiplier,
        image_height = 350*size_multiplier,
        vf = temperature_vf,
        dtmax = dtmax
    )

    radius_img = @time rendergeodesics(
        m, u, 2000.0, 
        d, 
        fov_factor=fov*size_multiplier, abstol=tolerance, reltol=tolerance,
        image_width = 350*size_multiplier,
        image_height = 350*size_multiplier,
        vf = radius_vf,
        dtmax = dtmax
    )

    redshift_img = @time rendergeodesics(
        m, u, 2000.0, 
        d, 
        fov_factor=fov*size_multiplier, abstol=tolerance, reltol=tolerance,
        image_width = 350*size_multiplier,
        image_height = 350*size_multiplier,
        vf = redshift_vf,
        dtmax = dtmax
    )
    # correcting for physical mass
    M_phys = mass*1.99e30
    r_isco_phys = AccretionFormulae.r_isco(M_phys, 0.998)
    r_g_phys = M_phys

    numerators =   AccretionFormulae.mass_scale_fraction.(M_phys, η_phys, edd_ratio_phys, r_isco_phys, r_g_phys, radius_img*M_phys)
    denominators = AccretionFormulae.mass_scale_fraction.(M, η, edd_ratio, R_isco, M, radius_img)

    fractions = numerators./denominators
    temperature_img .*= fractions

    # # autoscale
    # scale = maximum(filter(!isnan,temperature_img))
    # scale = floor(log(10, scale))
    # scale = 10^scale
    
    # scaling image
    scale = 1e7
    scalestr = @sprintf "%.E" scale
    new_img = reverse(temperature_img, dims=1)
    new_img ./= scale

    # hm = heatmap(new_img, aspect_ratio=1.0, size=(resolution*3/2, resolution), 
    # clim=(0,3)
    # )
    # # contour(new_img, aspect_ratio=1.0, size=(resolution*3/2, resolution), clim=(0,3))
    # title!("Temperature Scale = $scalestr, Mass = $mass M_☼")
    # # display(hm)
    return temperature_img, redshift_img, radius_img
end


function energy_histogram(;obs_angle=30.0)
    temperature_img, redshift_img, radius_img = temperature_render(obs_angle=obs_angle, size_multiplier=6, fov=6.0, tolerance=1.0E-9)

    # Define an energy array and a line profile array
    n_energies = 100
    # Minimum energy in keV
    e_min = 1.0
    # Maximum energy in keV
    e_max = 10.0
    energies = range(e_min, stop=e_max, length=n_energies)
    lineProfile = zeros(n_energies)

    # Define inner and outer radius of disk (could do this direclty using GeometricThinDisc)
    r_in = 5.0
    r_out = 10.0

    for (i, x) in enumerate(temperature_img)
        for (j, y) in enumerate(x)
            g = redshift_img[i][j]
            if !isnan(g)
                en = 6.4 * g
                index = trunc(Int, 1 + (en-e_min)*(n_energies-1)/(e_max-e_min))
                if (index >= 1) && (index <= n_energies)
                    # Example line profile where the emissivity is proportional to radius^-3 (an arbitrary choice)
                    if (radius_img[i][j] > r_in) && (radius_img[i][j] < r_out)
                        lineProfile[index] += radius_img[i][j]^-3 * g^4
                    end
                end
            end
        end
    end

    # plot line profile
    plot(energies, lineProfile, label=string(obs_angle), xlabel="Energy", ylabel="Line flux (arbitrary units)")

    # plot r image
    # heatmap(radius_img, aspect_ratio=1.0, clim=(0,50))

end

energy_histogram()