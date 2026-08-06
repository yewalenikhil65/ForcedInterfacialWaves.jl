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
# T₀ : time-independent (steady) contribution
# ═══════════════════════════════════════════════════════════════════════════════

"""
    gravity_T0(x, p)

Time-independent (steady) contribution:

```math
T_0(x) = \\frac{1}{\\pi(1+\\rho_r)}\\left[-\\pi\\sin(\\beta|x|)
+ \\int_0^\\infty e^{-y|x|}\\frac{y}{\\beta^2+y^2}\\,dy\\right]
```
"""
function gravity_T0(x::Real, p::PureGravityParams)
    xa = abs(x)
    β² = p.beta^2

    term_a = -π * sin(p.beta * xa)

    # This integrand decays exponentially — default order=7 is fine
    integrand = y -> exp(-y * xa) * y / (β² + y^2)
    term_b, _ = quadgk(integrand, 0.0, Inf;
                       atol=p.atol_transformed, rtol=p.rtol_transformed)

    return (term_a + term_b) / (π * (1.0 + p.rho_r))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Left region  x < t − δ_x  (behind front),  a = t − x > 0
# ═══════════════════════════════════════════════════════════════════════════════

"""
    gravity_T1_left(x, p)

First transient contribution (left region, closed-form):

```math
T_1^-(x) = -\\frac{\\sin(\\beta x)}{1+\\rho_r}
```
"""
@inline function gravity_T1_left(x::Real, p::PureGravityParams)
    return -sin(p.beta * x) / (1.0 + p.rho_r)
end

"""
    gravity_T2_left(x, t, a, p)

Exponentially-damped oscillatory integral (left region, ``a = t-x > 0``):

```math
T_2^- = -\\frac{4}{\\pi(1+\\rho_r)\\beta}\\int_0^{K_{\\max}} v^2\\,
e^{-2v^2 a + vt\\sqrt\\beta}\\,
\\frac{\\sqrt\\beta\\cos(vt\\sqrt\\beta) + (2v-\\sqrt\\beta)\\sin(vt\\sqrt\\beta)}
{\\beta + (2v-\\sqrt\\beta)^2}\\,dv
```
"""
function gravity_T2_left(x::Real, t::Real, a::Real, p::PureGravityParams)
    sβ  = p.sqrt_beta
    β   = p.beta
    tsβ = t * sβ

    integrand = v -> begin
        arg = v * tsβ
        num = sβ * cos(arg) + (2.0v - sβ) * sin(arg)
        den = β + (2.0v - sβ)^2
        return v^2 * exp(-2.0 * v^2 * a + v * tsβ) * num / den
    end

    # Exponential decay — converges quickly; order=7 sufficient
    I, _ = quadgk(integrand, 0.0, p.k_max_analytical;
                  atol=p.atol_transformed, rtol=p.rtol_transformed)

    return -4.0 / (π * (1.0 + p.rho_r) * β) * I
end

"""
    gravity_T3_left(x, t, a, p)

Fresnel-integral contribution (left region). No quadrature — uses [`fresnel_C`](@ref), [`fresnel_S`](@ref):

```math
T_3^- = \\frac{1}{\\pi(1+\\rho_r)\\sqrt\\beta}
\\left(1+\\frac{t}{2a}\\right)\\sqrt{\\frac{\\pi}{2a}}\\,
\\left[\\cos\\!\\left(\\frac{b^2}{a}\\right)\\!\\left(\\tfrac12 - C(X)\\right)
+ \\sin\\!\\left(\\frac{b^2}{a}\\right)\\!\\left(\\tfrac12 - S(X)\\right)\\right]
```

where ``b = t\\sqrt\\beta/2``, ``X = b\\sqrt{2/(\\pi a)}``.
"""
function gravity_T3_left(x::Real, t::Real, a::Real, p::PureGravityParams)
    sβ = p.sqrt_beta
    b  = 0.5 * t * sβ
    X  = b * sqrt(2.0 / (π * a))

    coeff = 1.0 / (π * (1.0 + p.rho_r) * sβ)
    pref  = (1.0 + t / (2.0 * a)) * sqrt(π / (2.0 * a))
    b2a   = b^2 / a

    return coeff * pref * (cos(b2a) * (0.5 - fresnel_C(X)) +
                           sin(b2a) * (0.5 - fresnel_S(X)))
end

"""
    gravity_T4_left(x, t, a, p)

Oscillatory integral (left region):

```math
T_4^- = -\\frac{1}{\\pi(1+\\rho_r)}\\int_0^{K_{\\max}}
\\frac{\\cos(v^2 a + vt\\sqrt\\beta)}{v + \\sqrt\\beta}\\,dv
```
"""
function gravity_T4_left(x::Real, t::Real, a::Real, p::PureGravityParams)
    sβ  = p.sqrt_beta
    tsβ = t * sβ

    # Oscillatory integrand — use higher order for better polynomial interpolation
    integrand = v -> cos(v^2 * a + v * tsβ) / (v + sβ)
    I, _ = quadgk(integrand, 0.0, p.k_max_analytical;
                  atol=p.atol_transformed, rtol=p.rtol_transformed, order=15)

    return -1.0 / (π * (1.0 + p.rho_r)) * I
end

"""
    gravity_analytical_left(x, t, p) → (η_total, η_steady, η_transient)

Complete analytical solution for ``x < t - \\delta_x`` (behind the wavefront):

```math
\\eta^-(x,t) = F_0\\left[T_0(x) + T_1^- + T_2^- + T_3^- + T_4^-\\right]
```

Returns a tuple `(η_total, η_steady, η_transient)` where `η_steady = F₀ T₀` and
`η_transient = F₀(T₁⁻ + T₂⁻ + T₃⁻ + T₄⁻)`.
"""
function gravity_analytical_left(x::Real, t::Real, p::PureGravityParams)
    a = t - x
    T0 = gravity_T0(x, p)
    T1 = gravity_T1_left(x, p)
    T2 = gravity_T2_left(x, t, a, p)
    T3 = gravity_T3_left(x, t, a, p)
    T4 = gravity_T4_left(x, t, a, p)

    η_s  = p.F0 * T0
    η_tr = p.F0 * (T1 + T2 + T3 + T4)
    return (η_s + η_tr, η_s, η_tr)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Right region  x > t + δ_x  (ahead of front),  a = t − x < 0
# ═══════════════════════════════════════════════════════════════════════════════

"""
    gravity_T1_right(x, p)

First transient contribution (right region, closed-form):

```math
T_1^+(x) = +\\frac{\\sin(\\beta x)}{1+\\rho_r}
```
"""
@inline function gravity_T1_right(x::Real, p::PureGravityParams)
    return sin(p.beta * x) / (1.0 + p.rho_r)
end

"""
    gravity_T2_right(x, t, a, p)

Exponentially-damped oscillatory integral (right region, ``a = t-x < 0``):

```math
T_2^+ = -\\frac{4}{\\pi(1+\\rho_r)\\beta}\\int_0^{K_{\\max}} v^2\\,
e^{2v^2 a - vt\\sqrt\\beta}\\,
\\frac{\\sqrt\\beta\\cos(vt\\sqrt\\beta) - (2v-\\sqrt\\beta)\\sin(vt\\sqrt\\beta)}
{\\beta + (2v-\\sqrt\\beta)^2}\\,dv
```

Note: ``a < 0`` so ``\\exp(2v^2 a)`` decays.
"""
function gravity_T2_right(x::Real, t::Real, a::Real, p::PureGravityParams)
    sβ  = p.sqrt_beta
    β   = p.beta
    tsβ = t * sβ

    integrand = v -> begin
        arg = v * tsβ
        num = sβ * cos(arg) - (2.0v - sβ) * sin(arg)
        den = β + (2.0v - sβ)^2
        return v^2 * exp(2.0 * v^2 * a - v * tsβ) * num / den
    end

    I, _ = quadgk(integrand, 0.0, p.k_max_analytical;
                  atol=p.atol_transformed, rtol=p.rtol_transformed)

    return -4.0 / (π * (1.0 + p.rho_r) * β) * I
end

"""
    gravity_T3_right(x, t, a, p)

Fresnel-integral contribution (right region, ``a < 0``). Uses ``|a|`` and flipped signs:

```math
T_3^+ = \\frac{1}{\\pi(1+\\rho_r)\\sqrt\\beta}
\\left(1+\\frac{t}{2a}\\right)\\sqrt{\\frac{\\pi}{2|a|}}\\,
\\left[\\cos\\!\\left(\\frac{b^2}{|a|}\\right)\\!\\left(\\tfrac12 + C(X)\\right)
+ \\sin\\!\\left(\\frac{b^2}{|a|}\\right)\\!\\left(\\tfrac12 + S(X)\\right)\\right]
```

where ``b = t\\sqrt\\beta/2``, ``X = b\\sqrt{2/(\\pi|a|)}``.
"""
function gravity_T3_right(x::Real, t::Real, a::Real, p::PureGravityParams)
    sβ   = p.sqrt_beta
    absa = abs(a)
    b    = 0.5 * t * sβ
    X    = b * sqrt(2.0 / (π * absa))

    coeff = 1.0 / (π * (1.0 + p.rho_r) * sβ)
    pref  = (1.0 + t / (2.0 * a)) * sqrt(π / (2.0 * absa))
    b2a   = b^2 / absa

    return coeff * pref * (cos(b2a) * (0.5 + fresnel_C(X)) +
                           sin(b2a) * (0.5 + fresnel_S(X)))
end

"""
    gravity_T4_right(x, t, a, p)

Oscillatory integral (right region). Same formula as ``T_4^-``; sign enters through ``a < 0``:

```math
T_4^+ = -\\frac{1}{\\pi(1+\\rho_r)}\\int_0^{K_{\\max}}
\\frac{\\cos(v^2 a + vt\\sqrt\\beta)}{v + \\sqrt\\beta}\\,dv
```
"""
function gravity_T4_right(x::Real, t::Real, a::Real, p::PureGravityParams)
    sβ  = p.sqrt_beta
    tsβ = t * sβ

    integrand = v -> cos(v^2 * a + v * tsβ) / (v + sβ)
    I, _ = quadgk(integrand, 0.0, p.k_max_analytical;
                  atol=p.atol_transformed, rtol=p.rtol_transformed, order=15)

    return -1.0 / (π * (1.0 + p.rho_r)) * I
end

"""
    gravity_analytical_right(x, t, p) → (η_total, η_steady, η_transient)

Complete analytical solution for ``x > t + \\delta_x`` (ahead of wavefront):

```math
\\eta^+(x,t) = F_0\\left[T_0(x) + T_1^+ + T_2^+ + T_3^+ + T_4^+\\right]
```
"""
function gravity_analytical_right(x::Real, t::Real, p::PureGravityParams)
    a = t - x
    T0 = gravity_T0(x, p)
    T1 = gravity_T1_right(x, p)
    T2 = gravity_T2_right(x, t, a, p)
    T3 = gravity_T3_right(x, t, a, p)
    T4 = gravity_T4_right(x, t, a, p)

    η_s  = p.F0 * T0
    η_tr = p.F0 * (T1 + T2 + T3 + T4)
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
