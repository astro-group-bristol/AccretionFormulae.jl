"""
    eⱽ(M, r, a, θ)
Modified from Cunningham et al. (1975) eq. (A2a):
```math
e^\\nu = \\sqrt{\\frac{\\Delta \\Sigma}{A}}.
```
"""
eⱽ(M, r, a, θ) = √(Σ(r, a, θ) * Δ(M, r, a) / A(M, r, a, θ))


"""
    eᶲ(M, r, a, θ) 
Modified from Cunningham et al. (1975) eq. (A2b):
```math
e^\\Phi = \\sin \\theta \\sqrt{\\frac{A}{\\Sigma}}.
```
"""
eᶲ(M, r, a, θ) = sin(θ) * √(A(M, r, a, θ) / Σ(r, a, θ))


"""
    ω(M, r, a, θ)
From Cunningham et al. (1975) eq. (A2c):
```math
\\omega = \\frac{2 a M r}{A}.
```
"""
ω(M, r, a, θ) = 2 * a * M * r / A(M, r, a, θ)


"""
    Ωₑ(M, r, a)
Coordinate angular velocity of an accreting gas. 
Taken from Cunningham et al. (1975) eq. (A7b):
```math
\\Omega_e = \\frac{\\sqrt{M}}{a \\sqrt{M} + r_e^{3/2}}.
```
# Notes
Fanton et al. (1997) use
```math
\\Omega_e = \\frac{\\sqrt{M}}{a \\sqrt{M} \\pm r_e^{3/2}},
```
where the sign is dependent on co- or contra-rotation. This function may be extended in the future to support this definition.
"""
Ωₑ(M, r, a) = √M / (r^1.5 + a * √M)


"""
    Vₑ(M, r, a, θ)  
Velocity of an accreting gas in a locally non-rotating reference frame (see Bardeen et al. 1973).
Taken from Cunningham et al. (1975) eq. (A7b):
```math
V_e = (\\Omega_e - \\omega) e^{\\Phi - \\nu}.
```
"""
Vₑ(M, r, a, θ) = (Ωₑ(M, r, a) - ω(M, r, a, θ)) * eᶲ(M, r, a, θ) / eⱽ(M, r, a, θ)


"""
    Lₑ(M, rms, a) 
Angular momentum of an accreting gas within ``r_ms``.
Taken from Cunningham et al. (1975) eq. (A11b):
```math
L_e = \\sqrt{M} \\frac{
        r_{\\text{ms}}^2 - 2 a \\sqrt{M r_{\\text{ms}}} + a^2
    }{
        r_{\\text{ms}}^{3/2} - 2 M \\sqrt{r_{\\text{ms}}} + a \\sqrt{M}
    }.
```
"""
Lₑ(M, rms, a) = √M * (rms^2 - 2 * a * √(M * rms) + a^2) / (rms^1.5 - 2 * M * √rms + a * √M)


"""
    H(M, rms, r, a)
Taken from Cunningham et al. (1975) eq. (A12e):
```math
H = \\frac{2 M r_e - a \\lambda_e}{\\Delta},
```
where we distinguing ``r_e`` as the position of the accreting gas. 
"""
H(M, rms, r, a) = (2 * M * r - a * Lₑ(M, rms, a)) / Δ(M, r, a)


"""
    γₑ(M, rms) 
Taken from Cunningham et al. (1975) eq. (A11c):
```math
\\gamma_e = \\sqrt{1 - \\frac{
        2M
    }{
        3 r_{\\text{ms}} 
    }}.
```
"""
γₑ(M, rms) = √(1 - (2 * M) / (3 * rms))


"""
    uʳ(M, rms, r)
Taken from Cunningham et al. (1975) eq. (A12b):
```math
u^r = - \\sqrt{\\frac{
        2M
    }{
        3 r_{\\text{ms}} 
    }} \\left(
        \\frac{ r_{\\text{ms}} }{r_e} - 1
    \\right)^{3/2}.
```
"""
uʳ(M, rms, r) = -√((2 * M) / (3 * rms)) * (rms / r - 1)^1.5


"""
    uᶲ(M, rms, r, a)
Taken from Cunningham et al. (1975) eq. (A12c):
```math
u^\\phi = \\frac{\\gamma_e}{r_e^2} \\left( 
        L_e + aH 
    \\right).
```
"""
uᶲ(M, rms, r, a) = γₑ(M, rms) / r^2 * (Lₑ(M, rms, a) + a * H(M, rms, r, a))


"""
    uᵗ(M, rms, r, a) 
Taken from Cunningham et al. (1975) eq. (A12b):
```math
u^t = \\gamma_e \\left(
        1 + \\frac{2 M (1 + H)}{r_e}
    \\right).
```
"""
uᵗ(M, rms, r, a) = γₑ(M, rms) * (1 + 2 * M * (1 + H(M, rms, r, a)) / r)


"""
Experimental API:
"""
function regular_pdotu_inv(L, M, r, a, θ)
    (eⱽ(M, r, a, θ) * √(1 - Vₑ(M, r, a, θ)^2)) / (1 - L * Ωₑ(M, r, a))
end
@inline function regular_pdotu_inv(u, p, m)
    @inbounds regular_pdotu_inv(p.L, m.M, u[2], m.a, u[3])
end

# Ignore plunging region for now 
# until we have an RMS function

#function plunging_p_dot_u(E, L, M, Q, rms, r, a, sign_r)
#    inv(
#        uᵗ(M, rms, r, a) - uᶲ(M, rms, r, a) * L -
#        sign_r * uʳ(M, rms, r) * Σδr_δλ(E, L, M, Q, r, a) / Δ(M, r, a)
#    )
#end
#function plunging_p_dot_u(u, p::CarterMethodBL{T}, sign_r) where {T}
#    plunging_p_dot_u(p.E, p.L, p.M, p.Q, p.rms, u[2], p.a, sign_r)
#end
#function plunging_p_dot_u(u, v, p::AbstractMetricParams{T}) where {T}
#    let r = u[2], a = p.a, M = p.M, rms = rms(p.M, p.a)
#        inv(uᵗ(M, rms, r, a) * v[1] - uᶲ(M, rms, r, a) * v[4] - uʳ(M, rms, r) * v[2])
#    end
#end
#@inline function redshift_function(val, λ, u, p::CarterMethodBL{T}, d) where {T}
#    @inbounds if u[2] > rms(p.M, p.a)
#        return regular_pdotu_inv(u, p)
#    else
#        return plunging_p_dot_u(u, p, λ < p.λr_change ? -1 : 1)
#    end
#end

@inline function redshift_function(m::CarterMethodBL{T}, u, p) where {T}
    regular_pdotu_inv(u, p, m)
end

@inline function redshift_function(m::AbstractMetricParams{T}, u, v) where {T}
    metric = GeodesicBase.metric(m, u)
    energy = GeodesicBase.E(metric, v)
    angmom = GeodesicBase.Lz(metric, v)

    # calculate momentum vector; only need 1 and 4 components
    # since currently u 2 and 3 are both 0
    # TODO: create a momentum function
    norm = metric[1,1]*metric[4,4] - metric[1,4]^2
    pt = - metric[4,4]*energy + metric[1,4]*angmom
    pϕ = metric[1,4]*energy + metric[1,1]*angmom

    # make momentum vector
    p = @SVector [
        pt/norm, 0, 0, pϕ/norm
    ]

    # TODO: this only works for Kerr
    disc_norm = (AccretionFormulae.eⱽ(m.M, u[2], m.a, u[3]) * √(1 - AccretionFormulae.Vₑ(m.M, u[2], m.a, u[3])^2))

    u_disc = @SVector [
        1/disc_norm, 0, 0, AccretionFormulae.Ωₑ(m.M, u[2], m.a)/disc_norm
    ]

    # use Tullio to do the einsum
    @tullio g := metric[i,j] * (u_disc[i]) * p[j]
    1/g
end

# value functions exports

function _redshift_guard(m::CarterMethodBL{T}, sol, max_time; kwargs...) where {T}
    redshift_function(m, sol.u[end], sol.prob.p)
end
function _redshift_guard(m::AbstractMetricParams{T}, sol, max_time; kwargs...) where {T}
    redshift_function(m, sol.u[end].x[2], sol.u[end].x[1])
end

const redshift = ValueFunction(_redshift_guard)