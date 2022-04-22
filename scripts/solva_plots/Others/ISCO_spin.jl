using AccretionFormulae
using Plots

function r_hor(a, M)
    (M +  √(M^2 - (a*M)^2))/M
end

M = 10*1.99e30
spin_vals_pos = LinRange(0, 1, 1000)
spin_vals_neg = LinRange(-1, 0 , 1000)
ISCO_pos_vals = AccretionFormulae.r_isco.(spin_vals_pos, M) ./ M
ISCO_neg_vals = reverse(AccretionFormulae.r_isco.(spin_vals_neg, M)) ./ M
EH_pos_vals = r_hor.(spin_vals_pos, M)

plt = plot(
            spin_vals_pos, 
            ISCO_pos_vals, 
            legend=:topleft, 
            label="\$r_{ISCO} \\textrm{Prograde}\$", 
            xlabel="a (M)", 
            ylabel="R (M)", 
            grid=false, 
            framestyle=:box
            )
plot!(spin_vals_pos, ISCO_neg_vals, label="\$r_{ISCO} \\textrm{Retrograde}\$")
plot!(spin_vals_pos, EH_pos_vals, label="\$r_{s}\$")
png(plt, "ISCO_EH_spin_comp.png")