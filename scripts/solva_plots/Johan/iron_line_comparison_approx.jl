include("../../iron_line_profile.jl")
include("../../arbitrary_f/iron_line_profile_approx.jl")

verticle_line_x = fill(6.4, 100)
verticle_line_y = LinRange(0, 1.05, 100)

x_vals_exact, bins_exact = iron_line_profile(                                            
                                            spin=0.8,                                            
                                            obs_angle=30,
                                            fov=12,
                                            # tolerance=1e-12,
                                            # size_multiplier=5,
                                            # dtmax=0.5,
                                            output="data",
                                            nbins=200,
                                            normalised=true
                                            )

x_vals_approx, bins_approx = iron_line_profile_approx(
                                                    spin=0.8,                                            
                                                    obs_angle=30,
                                                    fov=12,
                                                    # tolerance=1e-12,
                                                    # size_multiplier=5,
                                                    # dtmax=0.5,
                                                    output="data",
                                                    nbins=200,
                                                    normalised=true
                                                    )

plt = plot(
            x_vals_exact, 
            bins_exact, 
            xlabel="Energy (keV)",
            ylabel="Flux (Arbitrary Units)",
            label="Analytic f", 
            grid=false, 
            framestyle=:box, 
            ylims=(0,1.05)
            )
plot!(
        x_vals_approx, 
        bins_approx, 
        label="Numerical f", 
        legend=:topleft, 
        linestyle =:dash
        )
# adding vertical line at 6.4 keV
plot!(
    verticle_line_x, 
    verticle_line_y, 
    linestyle=:dash, 
    linecolor=:black, 
    label=false,
    ylims=(0,1.05),
    lw=2
    )
png(plt, "iron_line_approx_ana_comp.png")

#changing ϵ_3
vars = [-2,0,2]

labels = []
x_vals_vals = []
bins_vals = []

for var in vars
    x_vals_approx, bins_approx = iron_line_profile_approx(
                                                        ϵ_3=var,
                                                        spin=0.8,                                            
                                                        obs_angle=30,
                                                        fov=12,
                                                        # tolerance=1e-12,
                                                        # size_multiplier=5,
                                                        # dtmax=0.5,
                                                        output="data",
                                                        nbins=200,
                                                        normalised=false
                                                        )
    push!(bins_vals, bins_approx)
    push!(x_vals_vals, x_vals_approx)
    push!(labels, "ϵ₃ = $var")
end
fixed_labels=permutedims(labels)
plt = plot(
            x_vals_vals, 
            bins_vals, 
            xlabel="Energy (keV)",
            ylabel="Flux (Arbitrary Units)",
            label=fixed_labels, 
            grid=false, 
            framestyle=:box
            )
png(plt, "iron_line_eps_comp.png")

#changing α_22
vars = [-5,0,5]

labels = []
x_vals_vals = []
bins_vals = []

for var in vars
    x_vals_approx, bins_approx = iron_line_profile_approx(
                                                        α_22=var,
                                                        spin=0.8,                                            
                                                        obs_angle=30,
                                                        fov=12,
                                                        # tolerance=1e-12,
                                                        # size_multiplier=5,
                                                        # dtmax=0.5,
                                                        output="data",
                                                        nbins=200,
                                                        normalised=false
                                                        )
    push!(bins_vals, bins_approx)
    push!(x_vals_vals, x_vals_approx)
    push!(labels, "α₂₂ = $var")
end
fixed_labels=permutedims(labels)
plt = plot(
            x_vals_vals, 
            bins_vals, 
            xlabel="Energy (keV)",
            ylabel="Flux (Arbitrary Units)",
            label=fixed_labels, 
            grid=false, 
            framestyle=:box
            )
png(plt, "iron_line_a22_comp.png")

#changing α_13
vars = [-2,0,2]

labels = []
x_vals_vals = []
bins_vals = []

for var in vars
    x_vals_approx, bins_approx = iron_line_profile_approx(
                                                        α_13=var,
                                                        spin=0.8,                                            
                                                        obs_angle=30,
                                                        fov=12,
                                                        # tolerance=1e-12,
                                                        # size_multiplier=5,
                                                        # dtmax=0.5,
                                                        output="data",
                                                        nbins=200,
                                                        normalised=false
                                                        )
    push!(bins_vals, bins_approx)
    push!(x_vals_vals, x_vals_approx)
    push!(labels, "α₁₃ = $var")
end
fixed_labels=permutedims(labels)
plt = plot(
            x_vals_vals, 
            bins_vals, 
            xlabel="Energy (keV)",
            ylabel="Flux (Arbitrary Units)",
            label=fixed_labels, 
            grid=false, 
            framestyle=:box
            )
png(plt, "iron_line_a13_comp.png")

#changing α_52
vars = [-5,0,5]

labels = []
x_vals_vals = []
bins_vals = []

for var in vars
    x_vals_approx, bins_approx = iron_line_profile_approx(
                                                        α_52=var,
                                                        spin=0.8,                                            
                                                        obs_angle=30,
                                                        fov=12,
                                                        tolerance=1e-12,
                                                        size_multiplier=5,
                                                        dtmax=0.5,
                                                        output="data",
                                                        nbins=200,
                                                        normalised=false
                                                        )
    push!(bins_vals, bins_approx)
    push!(x_vals_vals, x_vals_approx)
    push!(labels, "α₅₂ = $var")
end
fixed_labels=permutedims(labels)
plt = plot(
            x_vals_vals, 
            bins_vals, 
            xlabel="Energy (keV)",
            ylabel="Flux (Arbitrary Units)",
            label=fixed_labels, 
            grid=false, 
            framestyle=:box
            )
png(plt, "iron_line_a52_comp.png")