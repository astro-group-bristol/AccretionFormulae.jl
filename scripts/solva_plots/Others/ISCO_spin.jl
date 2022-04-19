using AccretionFormulae

function r_hor(a, M)
    M +  √(M^2 - (a*M)^2)
end

M = 10*1.99e30
spin_vals_pos = LinRange(0, 1, 1000)
spin_vals_neg = LinRange(-1, 0 , 1000)
ISCO_pos_vals = AccretionFormulae.r_isco.(spin_vals_pos, M)
ISCO_neg_vals = reverse(AccretionFormulae.r_isco.(spin_vals_neg, M))
EH_pos_vals = r_hor.(spin_vals_pos, M)

plt = plot(
            spin_vals_pos, 
            ISCO_pos_vals, 
            legend=:topleft, 
            label="Prograde ISCO", 
            xlabel="Spin", 
            ylabel="Radius (m)", 
            grid=false, 
            framestyle=:box
            )
plot!(spin_vals_pos, ISCO_neg_vals, label="Retrograde ISCO")
plot!(spin_vals_pos, EH_pos_vals, label="Event Horizon Radius")
png(plt, "ISCO_EH_spin_comp.png")