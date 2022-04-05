using Revise
using GeodesicRendering
using AccretionGeometry
using CarterBoyerLindquist
using AccretionFormulae
using StaticArrays
using Plots

include("temperature_render.jl")
gr()

function flux(m, gp, max_time)
    g = AccretionFormulae.redshift(m, gp, max_time)
    r = gp.u[2]
    M = m.M
    a_star = m.a
    mdot = AccretionFormulae.mdot(M)
    flux = g^4 * AccretionFormulae.diss(mdot, r, a_star, M)
    return (flux, g)
end

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


function iron_line_profile(;
                            mass = 1,
                            spin = 0.998,
                            obs_angle = 85.0,
                            tolerance = 1e-8,
                            dtmax = 1000.0,
                            size_multiplier::Int64 = 1,
                            resolution = 400,
                            fov = 3.0,
                            output = "plot"
                            )
    hmap, cache, title = temperature_render(;
                                            mass = mass,
                                            spin = spin,
                                            obs_angle = obs_angle,
                                            tolerance = tolerance,
                                            dtmax = dtmax,
                                            size_multiplier = size_multiplier,
                                            resolution = resolution,
                                            fov = fov
                                            )

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
    n_bins = 100

    bins = zeros(n_bins)
    x_vals = LinRange(0, e_max, n_bins)

    new = en.(img, e_max, n_bins)

    for tuple in new
        if typeof(tuple) == Tuple{Float64, Int32}
            bins[tuple[2]] += tuple[1]
        end
    end

    # normalising to peak at 1
    bins ./= maximum(bins)
    plt = plot(x_vals, bins, xlims=(e_min, e_max), grid=false, framestyle=:box)
    xaxis!("Energy (keV)")
    yaxis!("Flux (Arbitrary Units)")

    if output == "plot"
        return plt, hmap, title, cache
    elseif output == "data"
        return x_vals, bins
    end
end
# plt, hmap, title, cache = iron_line_profile(tolerance=1e-12,
#                                             # size_multiplier=3,
#                                             # fov=12,
#                                             # dtmax=5,
#                                             obs_angle=30)
# display(plt)