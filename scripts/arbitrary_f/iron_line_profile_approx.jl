using Revise
using GeodesicRendering
using AccretionGeometry
using CarterBoyerLindquist
using AccretionFormulae
using StaticArrays
using Plots

include("temperature_render_approx.jl")
gr()

"""
Calculating the index for each flux
"""
function en(pixel, e_max, n_bins)
    if typeof(pixel) == Tuple{Float64, Float64}
        flux = pixel[1]
        g = pixel[2]
        en = 6.4*g
        index = n_bins*(en/e_max)

        # rounding down to the nearest bin
        index = Int32(round(index-0.5))
        return (flux, index)
    else
        return NaN
    end
end


function iron_line_profile_approx(;
                            mass = 1,
                            spin = 0.998,
                            obs_angle = 85.0,
                            tolerance = 1e-8,
                            dtmax = 1000.0,
                            size_multiplier::Int64 = 1,
                            resolution = 400,
                            fov = 3.0,
                            output = "plot",
                            normalised = true,
                            nbins = 100,
                            ϵ_3=0, 
                            α_13=0, 
                            α_22=0,
                            α_52=0
                            )

                            @show(ϵ_3)
                            @show(α_13)
                            @show(α_22)
                            @show(α_52)

    hmap, cache, title, fs, r_range = temperature_render_approx(;
                                            mass = mass,
                                            spin = spin,
                                            obs_angle = obs_angle,
                                            tolerance = tolerance,
                                            dtmax = dtmax,
                                            size_multiplier = size_multiplier,
                                            resolution = resolution,
                                            fov = fov,
                                            ϵ_3=ϵ_3, 
                                            α_13=α_13, 
                                            α_22=α_22,
                                            α_52=α_52
                                            )

    function flux(m, gp, max_time)
        g = AccretionFormulae.redshift(m, gp, max_time)
        r = gp.u[2]
        M = m.M
        a_star = m.a
        mdot = AccretionFormulae.mdot(M)
        flux = g^4 * diss_approx(mdot, r, a_star, M, fs, r_range)
        return (flux, g)
    end

    # define a common filter 
    filter = ConstPointFunctions.filter_intersected

    # define the value function
    flux_redshift_vf = (
        PointFunction(flux) ∘ filter
    )

    # generate the images from the cache
    img = GeodesicRendering.apply(flux_redshift_vf, cache)

    # parameters
    e_min = 0
    e_max = 10
    n_bins = nbins

    bins = zeros(n_bins)
    x_vals = LinRange(0, e_max, n_bins)

    new = en.(img, e_max, n_bins)

    for tuple in new
        if typeof(tuple) == Tuple{Float64, Int32}
            bins[tuple[2]] += tuple[1]
        end
    end

    # normalising to peak at 1
    if normalised
        bins ./= maximum(bins)
    end
    plt = plot(x_vals, bins, xlims=(e_min, e_max), grid=false, framestyle=:box)
    xaxis!("Energy (keV)")
    yaxis!("Flux (Arbitrary Units)")

    if output == "plot"
        return plt, hmap, title, cache
    elseif output == "data"
        return x_vals, bins
    end
end

# plt, hmap, title, cache = iron_line_profile_approx(tolerance=1e-12,
#                                             # size_multiplier=3,
#                                             # fov=12,
#                                             # dtmax=5,
#                                             obs_angle=30,
#                                             ϵ_3=0)
# display(plt)