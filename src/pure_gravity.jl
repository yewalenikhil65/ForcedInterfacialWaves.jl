"""
    Pure-gravity two-fluid IVP solution (α = 0), Figure 6.

    Two independent evaluations:
    1. Analytical CPV — T₀ + T₁⁻…T₄⁻ (left) / T₁⁺…T₄⁺ (right)
    2. Direct numerical CPV — combined integrand G(k; x, t) split around k = β
"""

# ═══════════════════════════════════════════════════════════════════════════════
# Fresnel integrals — using FresnelIntegrals.jl
# Convention:  C(x) = ∫₀ˣ cos(πt²/2) dt,  S(x) = ∫₀ˣ sin(πt²/2) dt
# Matches MATLAB's  fresnelc / fresnels exactly.
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fresnel_C(x)

Fresnel cosine integral (matches MATLAB's `fresnelc`):

```math
C(x) = \\int_0^x \\cos\\!\\left(\\frac{\\pi t^2}{2}\\right)dt
```

Wraps `FresnelIntegrals.fresnelc` — O(1) cost, no quadrature.
"""
@inline fresnel_C(x::Real) = fresnelc(x)

"""
    fresnel_S(x)

Fresnel sine integral (matches MATLAB's `fresnels`):

```math
S(x) = \\int_0^x \\sin\\!\\left(\\frac{\\pi t^2}{2}\\right)dt
```

Wraps `FresnelIntegrals.fresnels` — O(1) cost, no quadrature.
"""
@inline fresnel_S(x::Real) = fresnels(x)

# ═══════════════════════════════════════════════════════════════════════════════
# Unified T₀–T₄ functions (accept any x, t)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    T₀(x, p)

Time-independent (steady) contribution. Symmetric in x:

```math
T_0(x) = \\frac{1}{\\pi(1+\\rho_r)}\\left[-\\pi\\sin(\\beta|x|)
+ \\int_0^\\infty \\frac{y\\,e^{-|x|y}}{\\beta^2+y^2}\\,dy\\right]
```
"""
function T₀(x::Real, p::PureGravityParams)
    xa = abs(x)
    β² = p.beta^2
    term_a = -π * sin(p.beta * xa)
    integrand = y -> exp(-y * xa) * y / (β² + y^2)
    term_b, _ = quadgk(integrand, 0.0, Inf;
                       atol=p.atol_transformed, rtol=p.rtol_transformed)
    return (term_a + term_b) / (π * (1.0 + p.rho_r))
end

"""
    T₁(x, t, p)

Closed-form transient term. Sign determined by position relative to wavefront:

```math
T_1(x,t) = \\mp\\frac{\\sin(\\beta x)}{1+\\rho_r}
```
Upper sign (−) for x < t, lower sign (+) for x > t.
"""
@inline function T₁(x::Real, t::Real, p::PureGravityParams)
    s = x < t ? -1.0 : 1.0
    return s * sin(p.beta * x) / (1.0 + p.rho_r)
end

"""
    T₂(x, t, p)

Exponentially-damped oscillatory integral. Valid for any x ≠ t:

```math
T_2 = -\\frac{4}{\\pi(1+\\rho_r)\\beta}\\int_0^{K_{\\max}} v^2\\,
\\frac{\\exp(\\mp 2v^2 a \\pm vt\\sqrt\\beta)}{\\beta + (2v-\\sqrt\\beta)^2}
\\left[\\sqrt\\beta\\cos(vt\\sqrt\\beta) \\pm (2v-\\sqrt\\beta)\\sin(vt\\sqrt\\beta)\\right]dv
```
"""
function T₂(x::Real, t::Real, p::PureGravityParams)
    a   = t - x
    sβ  = p.sqrt_beta
    β   = p.beta
    tsβ = t * sβ
    s   = a > 0 ? 1.0 : -1.0   # +1 for left (x<t), -1 for right (x>t)

    integrand = v -> begin
        arg = v * tsβ
        num = sβ * cos(arg) + s * (2.0v - sβ) * sin(arg)
        den = β + (2.0v - sβ)^2
        return v^2 * exp(-s * 2.0 * v^2 * a + s * v * tsβ) * num / den
    end

    I, _ = quadgk(integrand, 0.0, p.k_max_analytical;
                  atol=p.atol_transformed, rtol=p.rtol_transformed)
    return -4.0 / (π * (1.0 + p.rho_r) * β) * I
end

"""
    T₃(x, t, p)

Fresnel-integral contribution. No quadrature — uses [`fresnel_C`](@ref), [`fresnel_S`](@ref):

```math
T_3 = \\frac{\\beta^{-1/2}}{\\pi(1+\\rho_r)}\\left(1+\\frac{t}{2a}\\right)\\sqrt{\\frac{\\pi}{2|a|}}
\\left[\\cos\\!\\left(\\frac{b^2}{|a|}\\right)\\left(\\frac{1}{2}\\mp C(X)\\right)
+ \\sin\\!\\left(\\frac{b^2}{|a|}\\right)\\left(\\frac{1}{2}\\mp S(X)\\right)\\right]
```
where ``a = t-x``, ``b = t\\sqrt\\beta/2``, ``X = b\\sqrt{2/(\\pi|a|)}``.
"""
function T₃(x::Real, t::Real, p::PureGravityParams)
    a    = t - x
    absa = abs(a)
    sβ   = p.sqrt_beta
    b    = 0.5 * t * sβ
    X    = b * sqrt(2.0 / (π * absa))
    s    = a > 0 ? 1.0 : -1.0   # +1 for left, -1 for right

    coeff = 1.0 / (π * (1.0 + p.rho_r) * sβ)
    pref  = (1.0 + t / (2.0 * a)) * sqrt(π / (2.0 * absa))
    b2a   = b^2 / absa

    return coeff * pref * (cos(b2a) * (0.5 - s * fresnel_C(X)) +
                           sin(b2a) * (0.5 - s * fresnel_S(X)))
end

"""
    T₄(x, t, p)

Oscillatory integral. Same formula for both regions — sign enters through ``a = t-x``:

```math
T_4 = -\\frac{1}{\\pi(1+\\rho_r)}\\int_0^{K_{\\max}}
\\frac{\\cos(av^2 + vt\\sqrt\\beta)}{v + \\sqrt\\beta}\\,dv
```
"""
function T₄(x::Real, t::Real, p::PureGravityParams)
    a   = t - x
    sβ  = p.sqrt_beta
    tsβ = t * sβ

    integrand = v -> cos(v^2 * a + v * tsβ) / (v + sβ)
    I, _ = quadgk(integrand, 0.0, p.k_max_analytical;
                  atol=p.atol_transformed, rtol=p.rtol_transformed, order=15)
    return -1.0 / (π * (1.0 + p.rho_r)) * I
end

# ═══════════════════════════════════════════════════════════════════════════════
# Backward-compatible wrappers (deprecated)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    gravity_T0(x, p)

Deprecated: use [`T₀`](@ref) instead.
"""
gravity_T0(x::Real, p::PureGravityParams) = T₀(x, p)

gravity_T1_left(x::Real, p::PureGravityParams) = -sin(p.beta * x) / (1.0 + p.rho_r)
gravity_T1_right(x::Real, p::PureGravityParams) = sin(p.beta * x) / (1.0 + p.rho_r)

# For T2–T4, the old API passed `a` explicitly — just forward to unified version
gravity_T2_left(x::Real, t::Real, a::Real, p::PureGravityParams) = T₂(x, t, p)
gravity_T3_left(x::Real, t::Real, a::Real, p::PureGravityParams) = T₃(x, t, p)
gravity_T4_left(x::Real, t::Real, a::Real, p::PureGravityParams) = T₄(x, t, p)
gravity_T2_right(x::Real, t::Real, a::Real, p::PureGravityParams) = T₂(x, t, p)
gravity_T3_right(x::Real, t::Real, a::Real, p::PureGravityParams) = T₃(x, t, p)
gravity_T4_right(x::Real, t::Real, a::Real, p::PureGravityParams) = T₄(x, t, p)

"""
    gravity_analytical_left(x, t, p) → (η_total, η_steady, η_transient)

Deprecated: use `solve(ForcedGravityProblem(p, x, t))` instead.
"""
function gravity_analytical_left(x::Real, t::Real, p::PureGravityParams)
    t0 = T₀(x, p)
    η_s  = p.F0 * t0
    η_tr = p.F0 * (T₁(x, t, p) + T₂(x, t, p) + T₃(x, t, p) + T₄(x, t, p))
    return (η_s + η_tr, η_s, η_tr)
end

"""
    gravity_analytical_right(x, t, p) → (η_total, η_steady, η_transient)

Deprecated: use `solve(ForcedGravityProblem(p, x, t))` instead.
"""
function gravity_analytical_right(x::Real, t::Real, p::PureGravityParams)
    t0 = T₀(x, p)
    η_s  = p.F0 * t0
    η_tr = p.F0 * (T₁(x, t, p) + T₂(x, t, p) + T₃(x, t, p) + T₄(x, t, p))
    return (η_s + η_tr, η_s, η_tr)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Direct numerical CPV
# ═══════════════════════════════════════════════════════════════════════════════

"""
    gravity_combined_integrand(k, x, t, p)

Three-term integrand ``G(k; x, t)`` for direct numerical CPV evaluation.
Terms 1 and 2 are individually singular at ``k = \\beta``; they cancel when combined:

```math
G = \\frac{\\cos(kx)}{\\pi(1+\\rho_r)(k-\\beta)}
  - \\frac{k\\cos[k(t-x)-t\\sqrt{\\beta k}]}{2\\pi(1-\\rho_r)(k-\\sqrt{\\beta k})}
  - \\frac{k\\cos[k(t-x)+t\\sqrt{\\beta k}]}{2\\pi(1-\\rho_r)(k+\\sqrt{\\beta k})}
```
"""
@inline function gravity_combined_integrand(k::Real, x::Real, t::Real, p::PureGravityParams)
    @fastmath begin
        sqrt_bk = sqrt(p.beta * k)
        a_phase = k * (t - x)
        t_sqrt_bk = t * sqrt_bk

        inv_pi_1pr = 1.0 / (π * (1.0 + p.rho_r))
        inv_2pi_1mr = 1.0 / (2.0π * (1.0 - p.rho_r))

        term1 = cos(k * x) * inv_pi_1pr / (k - p.beta)
        term2 = -k * cos(a_phase - t_sqrt_bk) * inv_2pi_1mr / (k - sqrt_bk)
        term3 = -k * cos(a_phase + t_sqrt_bk) * inv_2pi_1mr / (k + sqrt_bk)

        return term1 + term2 + term3
    end
end

"""
    gravity_numerical_cpv(x, t, p)

Direct numerical Cauchy principal value evaluation, split symmetrically around ``k = \\beta``:

```math
\\eta_{\\mathrm{CPV}}(x,t) = F_0\\left[
\\int_0^{\\beta-\\varepsilon} G(k;x,t)\\,dk +
\\int_{\\beta+\\varepsilon}^{K_{\\max}} G(k;x,t)\\,dk
\\right]
```

where ``G`` is [`gravity_combined_integrand`](@ref). Uses `order=15` for the oscillatory interval above the pole.
"""
function gravity_numerical_cpv(x::Real, t::Real, p::PureGravityParams)
    integrand = k -> gravity_combined_integrand(k, x, t, p)
    ε = p.epsilon_cpv

    # Below pole: short interval, default order fine
    I_lo, _ = quadgk(integrand, 0.0, p.beta - ε;
                     atol=p.atol_cpv, rtol=p.rtol_cpv)

    # Above pole: long oscillatory interval, higher order helps
    I_hi, _ = quadgk(integrand, p.beta + ε, p.k_max_cpv;
                     atol=p.atol_cpv, rtol=p.rtol_cpv, order=15)

    return p.F0 * (I_lo + I_hi)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Profile helpers
# ═══════════════════════════════════════════════════════════════════════════════

"""
    compute_gravity_profile(x_grid, t, p) → (eta_analytical, eta_steady, eta_transient, eta_cpv)

Batch evaluation over spatial grid. Points within |x−t| ≤ front_band are NaN.
Pre-allocates a segment buffer for repeated quadgk calls to reduce allocations.
"""
function compute_gravity_profile(x_grid::AbstractVector{<:Real}, t::Real, p::PureGravityParams)
    N = length(x_grid)
    eta_analytical = fill(NaN, N)
    eta_steady     = Vector{Float64}(undef, N)
    eta_transient  = fill(NaN, N)
    eta_cpv        = fill(NaN, N)

    # Pre-allocate segment buffer for reuse across quadgk calls
    segbuf = alloc_segbuf()

    Threads.@threads for ix in 1:N
        xv = x_grid[ix]
        a  = t - xv

        eta_steady[ix] = p.F0 * gravity_T0(xv, p)

        if a > p.front_band
            η_tot, _, η_tr = gravity_analytical_left(xv, t, p)
            eta_analytical[ix] = η_tot
            eta_transient[ix]  = η_tr
            eta_cpv[ix] = gravity_numerical_cpv(xv, t, p)
        elseif a < -p.front_band
            η_tot, _, η_tr = gravity_analytical_right(xv, t, p)
            eta_analytical[ix] = η_tot
            eta_transient[ix]  = η_tr
            eta_cpv[ix] = gravity_numerical_cpv(xv, t, p)
        end
    end

    return eta_analytical, eta_steady, eta_transient, eta_cpv
end

"""
    make_gravity_xgrid(p; Nx=2001)

Nondimensional spatial grid for pure-gravity case, excluding x ≈ 0.
"""
function make_gravity_xgrid(p::PureGravityParams; Nx::Int=2001)
    xg = collect(range(-p.L / 2.0, p.L / 2.0; length=Nx))
    filter!(x -> abs(x) > 1.0e-6, xg)
    return xg
end
