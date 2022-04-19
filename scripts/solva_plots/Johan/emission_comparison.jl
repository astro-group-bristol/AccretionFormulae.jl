using AccretionFormulae
include("../../arbitrary_f/alt-emissivity-deviations.jl")

a = 0.97
M = 1.0

# creating arbitrary function
fs, r_range = init()
fs_epsneg, r_range_epsneg = init(ϵ_3=2)
fs_epspos, r_range_epspos = init(ϵ_3=-2)

r_vals = r_range
r_range = collect(r_range)
@show(r_range)
pop!(r_range)

f_vals_analytic = AccretionFormulae.temperature.(r_range, a, M)
f_vals_approx = temperature_approx.(r_vals, a, M, fs, r_range)
f_vals_approx_epsneg = temperature_approx.(r_vals, a, M, fs_epsneg, r_range_epsneg)
f_vals_approx_epspos = temperature_approx.(r_vals, a, M, fs_epspos, r_range_epspos)

plt = plot(r_range, f_vals_analytic, label="Analytic Solution")
plot!(r_range, f_vals_approx)
plot!(r_range, f_vals_approx_epsneg)
plot!(r_range, f_vals_approx_epspos)
