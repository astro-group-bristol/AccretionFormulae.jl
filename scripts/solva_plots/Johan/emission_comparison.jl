using AccretionFormulae
include("../../arbitrary_f/alt-emissivity-deviations.jl")

a = 0.97
M = 1.0

# epsilon
# creating arbitrary function
fs, r_range = init()
fs_epsneg, r_range_epsneg = init(ϵ_3=-2)
fs_epspos, r_range_epspos = init(ϵ_3=2)

r_vals = collect(r_range)
pop!(r_vals)

f_vals_analytic = AccretionFormulae.flux.(r_range, a, M)
f_vals_analytic ./= maximum(f_vals_analytic)
f_vals_approx = []
f_vals_approx_epsneg = []
f_vals_approx_epspos = []
for r in r_range
    push!(f_vals_approx, f_approx_func(r, fs, r_range, a=a, M=M))
    push!(f_vals_approx_epsneg, f_approx_func(r, fs_epsneg, r_range_epsneg, a=a, M=M))
    push!(f_vals_approx_epspos, f_approx_func(r, fs_epspos, r_range_epspos, a=a, M=M))
end

plt = plot(r_range, f_vals_analytic, label="Analytic Solution", grid=false, framestyle=:box, ylims=(0,1.05), xlabel="Radius (M)", ylabel="f (Arbitrary Units)")
plot!(r_range, f_vals_approx, label="ϵ₃=0", linestyle=:dash)
plot!(r_range, f_vals_approx_epsneg, label="ϵ₃=-2")
plot!(r_range, f_vals_approx_epspos, label="ϵ₃=2")

png(plt, "emission_comparison_ana_arb_eps.png")

# alpha 22
# creating arbitrary function
fs, r_range = init()
fs_epsneg, r_range_epsneg = init(α_22 =-5)
fs_epspos, r_range_epspos = init(α_22 =5)

r_vals = collect(r_range)
pop!(r_vals)

f_vals_approx = []
f_vals_approx_epsneg = []
f_vals_approx_epspos = []
for r in r_range
    push!(f_vals_approx, f_approx_func(r, fs, r_range, a=a, M=M))
    push!(f_vals_approx_epsneg, f_approx_func(r, fs_epsneg, r_range_epsneg, a=a, M=M))
    push!(f_vals_approx_epspos, f_approx_func(r, fs_epspos, r_range_epspos, a=a, M=M))
end

plt = plot(r_range, f_vals_analytic, label="Analytic Solution", grid=false, framestyle=:box, ylims=(0,1.05), xlabel="Radius (M)", ylabel="f (Arbitrary Units)")
plot!(r_range, f_vals_approx, label="α₂₂=0", linestyle=:dash)
plot!(r_range, f_vals_approx_epsneg, label="α₂₂=-5")
plot!(r_range, f_vals_approx_epspos, label="α₂₂=5")

png(plt, "emission_comparison_ana_arb_alpha22.png")

# alpha 13
# creating arbitrary function
fs, r_range = init()
fs_epsneg, r_range_epsneg = init(α_13 =-2)
fs_epspos, r_range_epspos = init(α_13 =2)

r_vals = collect(r_range)
pop!(r_vals)

f_vals_approx = []
f_vals_approx_epsneg = []
f_vals_approx_epspos = []
for r in r_range
    push!(f_vals_approx, f_approx_func(r, fs, r_range, a=a, M=M))
    push!(f_vals_approx_epsneg, f_approx_func(r, fs_epsneg, r_range_epsneg, a=a, M=M))
    push!(f_vals_approx_epspos, f_approx_func(r, fs_epspos, r_range_epspos, a=a, M=M))
end

plt = plot(r_range, f_vals_analytic, label="Analytic Solution", grid=false, framestyle=:box, ylims=(0,1.05), xlabel="Radius (M)", ylabel="f (Arbitrary Units)")
plot!(r_range, f_vals_approx, label="α₁₃ =0", linestyle=:dash)
plot!(r_range, f_vals_approx_epsneg, label="α₁₃=-2")
plot!(r_range, f_vals_approx_epspos, label="α₁₃=2")

png(plt, "emission_comparison_ana_arb_alpha13.png")

# alpha 52
# creating arbitrary function
fs, r_range = init()
fs_epsneg, r_range_epsneg = init(α_52 =-2)
fs_epspos, r_range_epspos = init(α_52 =2)

r_vals = collect(r_range)
pop!(r_vals)

f_vals_approx = []
f_vals_approx_epsneg = []
f_vals_approx_epspos = []
for r in r_range
    push!(f_vals_approx, f_approx_func(r, fs, r_range, a=a, M=M))
    push!(f_vals_approx_epsneg, f_approx_func(r, fs_epsneg, r_range_epsneg, a=a, M=M))
    push!(f_vals_approx_epspos, f_approx_func(r, fs_epspos, r_range_epspos, a=a, M=M))
end

plt = plot(r_range, f_vals_analytic, label="Analytic Solution", grid=false, framestyle=:box, ylims=(0,1.05), xlabel="Radius (M)", ylabel="f (Arbitrary Units)")
plot!(r_range, f_vals_approx, label="α₅₂ =0", linestyle=:dash)
plot!(r_range, f_vals_approx_epsneg, label="α₅₂ =-2")
plot!(r_range, f_vals_approx_epspos, label="α₅₂ =2")

png(plt, "emission_comparison_ana_arb_alpha52.png")