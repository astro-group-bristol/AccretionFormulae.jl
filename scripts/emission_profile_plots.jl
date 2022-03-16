# imports
using AccretionFormulae
using Plots
using Measures

# initialising plot backend
multiplier = 3
font_size = multiplier*10
line_width = multiplier*2
margin = multiplier*10mm
gr(
    legendfontsize=font_size, 
    margin=margin, 
    right_margin=margin*0.5,
    framestyle=:box, 
    grid=false, 
    lw=line_width, 
    xlabel="R (M)",
    xguidefontsize=font_size,
    ylabel="Temperature (10⁶ K)",
    yguidefontsize=font_size,
    tickfontsize=font_size,
    thickness_scaling=1
    )

"""
Automatically determines the labels for various lines to go into the legend
"""
function label(value, name, units)
    scale = floor(log(10, value)) - 1
    if scale > 4
        scale = Int32(scale)
        value = "\$10^{$scale}\$"
    else
        value = "\$$value\$"
    end
    label = "$name = $value $(units) "
end

# constants
M_☼ = 1.99e30
G = 6.67e-11
c = 3e8
M = 10*M_☼
a_star = 0.0

# masses and spins to be changed
masses = [10, 100, 200, 10e6]
spins = [0.0, 0.5, 0.8, 0.95, 0.998]

# legend labels, colours and styles
mass_labels = permutedims(label.(masses, "M", "\$\\mathrm{M}_{\\odot}\$"))
spin_labels = permutedims(label.(spins, "Spin", ""))
colors = [:black :blue :red :green :purple]
styles = [:solid :dash :dot :dashdot :dashdotdot]

# changing the mass
temperatures_mass = []
for mass in masses
    mass *= M_☼
    R_isco = AccretionFormulae.r_isco(a_star, mass)
    r_vals = LinRange(R_isco, 50*mass, 1000)
    temperature_vals = AccretionFormulae.temperature.(r_vals, a_star, mass)
    temperature_vals ./= 1e6
    push!(temperatures_mass, temperature_vals)
end

# changing the radius units to gravitational radii
R_isco_mass = AccretionFormulae.r_isco(a_star, M)
r_vals_mass = LinRange(R_isco_mass, 50*M, 1000)
r_vals_mass /= M

# changing the spin
temperatures_spin = []
r_vals_spins = []
for spin in spins
    R_isco = AccretionFormulae.r_isco(spin, M)
    r_vals = LinRange(R_isco, 30*M, 1000)
    temperature_vals = AccretionFormulae.temperature.(r_vals, spin, M)
    temperature_vals ./= 1e6
    # every line has a ISCO, so need to create r_vals in gravitational radii individually
    r_vals = r_vals ./ M
    push!(r_vals_spins, r_vals)
    push!(temperatures_spin, temperature_vals)
end

# plots
changing_mass_plot = plot(
    r_vals_mass, 
    temperatures_mass,
    label=mass_labels,
    linecolor = :black,
    linestyle = styles
    )

changing_spin_plot = plot(
    r_vals_spins, 
    temperatures_spin,
    label=spin_labels,
    linecolor = :black,
    linestyle = styles
    )

combined_plot = plot(changing_mass_plot, changing_spin_plot, layout=(1,2), size=multiplier.*(1200,400))

# saving image
png(combined_plot, "emission_profiles.png")
display(combined_plot)