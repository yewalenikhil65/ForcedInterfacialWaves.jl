"""
    validate_integrals.jl

    Open validation of individual integrals:
    - Pure gravity: T₀, T₁, T₂, T₃, T₄ individually for left/right regions
    - Capillary-gravity: partial integrals I₁, I₂, I₃ and G(x)
    - Cross-check: analytical terms sum to CPV for α=0
    - Cross-check: IVP converges to steady for α>0

    This script prints all intermediate values for reproducibility
    and comparison against MATLAB.

    Usage:
        julia --project=. scripts/validate_integrals.jl
"""

push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using WavesDelta
using Printf

println("="^72)
println("  WavesDelta.jl — INTEGRAL-LEVEL VALIDATION")
println("  Each integral validated independently for MATLAB equivalence")
println("="^72)

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1: Pure-Gravity (α = 0) — Individual T-terms
# ═══════════════════════════════════════════════════════════════════════════════

println("\n", "─"^72)
println("  PURE GRAVITY (α = 0)  —  Individual analytical terms")
println("─"^72)

pg = compute_gravity_parameters()
@printf("\n  Parameters:\n")
@printf("    β          = %.15e\n", pg.beta)
@printf("    √β         = %.15e\n", pg.sqrt_beta)
@printf("    F₀         = %.15e\n", pg.F0)
@printf("    K_max      = %.1f\n", pg.k_max_analytical)
@printf("    ε_cpv      = %.1e\n", pg.epsilon_cpv)
@printf("    front_band = %.1f\n", pg.front_band)

# Test at multiple (x, t) combinations
test_cases_left = [
    (x=-5.0, t_dim=1.0),
    (x=-2.0, t_dim=1.0),
    (x= 0.5, t_dim=1.0),
    (x=-3.0, t_dim=0.60),
]

println("\n  ── LEFT REGION (x < t − δ) ──")
println("  Each T-term is an integral; sum multiplied by F₀ gives η_transient.")
println()

for case in test_cases_left
    x = case.x
    t_dim = case.t_dim
    t = t_dim / pg.t_c
    a = t - x

    if a <= pg.front_band
        continue
    end

    @printf("  x = %6.2f,  t_dim = %.2f s,  t = %.6f,  a = t−x = %.6f\n", x, t_dim, t, a)
    println("  " * "·"^60)

    T0 = gravity_T0(x, pg)
    T1 = gravity_T1_left(x, pg)
    T2 = gravity_T2_left(x, t, a, pg)
    T3 = gravity_T3_left(x, t, a, pg)
    T4 = gravity_T4_left(x, t, a, pg)

    @printf("    T₀(x)       = %+.15e\n", T0)
    @printf("    T₁⁻(x)      = %+.15e   [closed-form: −sin(βx)/(1+ρr)]\n", T1)
    @printf("    T₂⁻(x,t)    = %+.15e   [quadgk integral, 0 to K_max]\n", T2)
    @printf("    T₃⁻(x,t)    = %+.15e   [fresnelc/fresnels, no quadrature]\n", T3)
    @printf("    T₄⁻(x,t)    = %+.15e   [quadgk integral, 0 to K_max]\n", T4)

    η_s  = pg.F0 * T0
    η_tr = pg.F0 * (T1 + T2 + T3 + T4)
    η_tot = η_s + η_tr

    @printf("    η_steady     = %+.15e  (= F₀·T₀)\n", η_s)
    @printf("    η_transient  = %+.15e  (= F₀·(T₁+T₂+T₃+T₄))\n", η_tr)
    @printf("    η_analytical = %+.15e  (= η_s + η_tr)\n", η_tot)

    # Cross-check with numerical CPV
    η_cpv = gravity_numerical_cpv(x, t, pg)
    @printf("    η_CPV        = %+.15e  [direct numerical CPV]\n", η_cpv)
    @printf("    |analytical − CPV| = %.4e\n\n", abs(η_tot - η_cpv))
end

# Right region
test_cases_right = [
    (x_offset=3.0, t_dim=1.0),
    (x_offset=5.0, t_dim=1.0),
    (x_offset=8.0, t_dim=1.0),
]

println("  ── RIGHT REGION (x > t + δ) ──\n")

for case in test_cases_right
    t_dim = case.t_dim
    t = t_dim / pg.t_c
    x = t + case.x_offset
    a = t - x  # a < 0

    @printf("  x = %6.2f,  t_dim = %.2f s,  t = %.6f,  a = t−x = %.6f\n", x, t_dim, t, a)
    println("  " * "·"^60)

    T0 = gravity_T0(x, pg)
    T1 = gravity_T1_right(x, pg)
    T2 = gravity_T2_right(x, t, a, pg)
    T3 = gravity_T3_right(x, t, a, pg)
    T4 = gravity_T4_right(x, t, a, pg)

    @printf("    T₀(x)       = %+.15e\n", T0)
    @printf("    T₁⁺(x)      = %+.15e   [closed-form: +sin(βx)/(1+ρr)]\n", T1)
    @printf("    T₂⁺(x,t)    = %+.15e   [quadgk integral]\n", T2)
    @printf("    T₃⁺(x,t)    = %+.15e   [fresnelc/fresnels]\n", T3)
    @printf("    T₄⁺(x,t)    = %+.15e   [quadgk integral]\n", T4)

    η_s  = pg.F0 * T0
    η_tr = pg.F0 * (T1 + T2 + T3 + T4)
    η_tot = η_s + η_tr

    @printf("    η_steady     = %+.15e\n", η_s)
    @printf("    η_transient  = %+.15e\n", η_tr)
    @printf("    η_analytical = %+.15e\n", η_tot)

    η_cpv = gravity_numerical_cpv(x, t, pg)
    @printf("    η_CPV        = %+.15e\n", η_cpv)
    @printf("    |analytical − CPV| = %.4e\n\n", abs(η_tot - η_cpv))
end

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2: Capillary-Gravity (α > 0) — Partial integrals
# ═══════════════════════════════════════════════════════════════════════════════

println("\n", "─"^72)
println("  CAPILLARY-GRAVITY (α > 0)  —  Partial integrals I₁, I₂, I₃")
println("─"^72)

p = compute_cg_parameters()
@printf("\n  Parameters:\n")
@printf("    α          = %.15e\n", p.alpha)
@printf("    ρ_r        = %.15e\n", p.rho_r)
@printf("    k_s        = %.15e\n", p.k_s)
@printf("    k_l        = %.15e\n", p.k_l)
@printf("    F₀         = %.15e\n", p.F0)
@printf("    ε_pv       = %.1e\n", p.epsilon_pv)
@printf("    k_max      = %s\n", p.k_max == Inf ? "∞" : string(p.k_max))
@printf("    k_max_stdy = %s\n", p.k_max_steady == Inf ? "∞" : string(p.k_max_steady))

println("\n  ── Integrand values at sample k ──")
println("  (Verifies Julia integrand functions match MATLAB pointwise)\n")

x_test, t_test = 3.0, 110.0
@printf("  At x = %.1f, t = %.1f:\n", x_test, t_test)
for k in [0.5, 1.0, p.k_s - 0.01, p.k_s + 0.01, 3.0, p.k_l - 0.01, p.k_l + 0.01, 10.0]
    Is = cg_integrand_steady(k, x_test, p)
    I3 = cg_integrand_I3(k, x_test, t_test, p)
    I4 = cg_integrand_I4(k, x_test, t_test, p)
    It = cg_combined_integrand(k, x_test, t_test, p)
    @printf("    k=%8.4f  I_s=%+.8e  I₃=%+.8e  I₄=%+.8e  𝕀=%+.8e\n", k, Is, I3, I4, It)
end

println("\n  ── Partial integrals over sub-intervals ──")
println("  Matches MATLAB's integral_1, integral_2, integral_3\n")

cg_test_points = [
    (x=-3.0, t_dim=3.0),
    (x= 1.0, t_dim=3.0),
    (x= 3.0, t_dim=3.0),
    (x= 5.0, t_dim=3.0),
    (x= 3.0, t_dim=0.15),
]

for case in cg_test_points
    x = case.x
    t_dim = case.t_dim
    t = t_dim / p.t_c

    I1, I2, I3 = cg_partial_integrals(x, t, p)
    Itotal = I1 + I2 + I3
    η_ivp = -p.F0 / (2π) * Itotal

    @printf("  x=%6.2f, t_dim=%.2f s (t=%.4f):\n", x, t_dim, t)
    @printf("    I₁ = ∫₀^{k_s−ε}        = %+.15e\n", I1)
    @printf("    I₂ = ∫_{k_s+ε}^{k_l−ε} = %+.15e\n", I2)
    @printf("    I₃ = ∫_{k_l+ε}^∞       = %+.15e\n", I3)
    @printf("    I₁+I₂+I₃               = %+.15e\n", Itotal)
    @printf("    η = −F₀/(2π)·(sum)      = %+.15e\n", η_ivp)

    # Steady comparison at large time
    if t_dim >= 3.0
        η_s = steady_surface_elevation(x, p)
        @printf("    η_steady                = %+.15e\n", η_s)
        @printf("    |η_IVP − η_steady|      = %.4e\n", abs(η_ivp - η_s))
    end

    # G(x) integral
    Gx = cg_Gx_integral(x, p)
    @printf("    G(x)                    = %+.15e\n\n", Gx)
end

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3: Fresnel validation (no quadrature — direct closed-form)
# ═══════════════════════════════════════════════════════════════════════════════

println("─"^72)
println("  FRESNEL INTEGRALS  —  Used in T₃ (closed-form, faster than quadgk)")
println("─"^72)
println()
@printf("  fresnelc(0)   = %.15e  (exact: 0)\n", fresnel_C(0.0))
@printf("  fresnels(0)   = %.15e  (exact: 0)\n", fresnel_S(0.0))
@printf("  fresnelc(1)   = %.15e  (NIST: 0.779893400376823)\n", fresnel_C(1.0))
@printf("  fresnels(1)   = %.15e  (NIST: 0.438259147390355)\n", fresnel_S(1.0))
@printf("  fresnelc(2)   = %.15e  (NIST: 0.488253406075341)\n", fresnel_C(2.0))
@printf("  fresnels(2)   = %.15e  (NIST: 0.343415678363698)\n", fresnel_S(2.0))
@printf("  fresnelc(5)   = %.15e  (→ 0.5)\n", fresnel_C(5.0))
@printf("  fresnels(5)   = %.15e  (→ 0.5)\n", fresnel_S(5.0))

println("\n  Note: fresnelc/fresnels from FresnelIntegrals.jl uses the erf")
println("  representation — O(1) cost per evaluation, no adaptive quadrature.")
println("  This is identical to MATLAB's built-in fresnelc/fresnels.")

# ═══════════════════════════════════════════════════════════════════════════════

println("\n", "="^72)
println("  VALIDATION COMPLETE")
println("  To reproduce in MATLAB: evaluate the same integrands/integrals at")
println("  the same (x, t, k) points and compare to 12+ digits.")
println("="^72)
