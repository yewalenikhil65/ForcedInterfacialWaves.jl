# validation_table.jl — Generate reference values at specific (x,t) points
# for MATLAB-Julia comparison.
#
# Run: julia --project=. test/validation_table.jl

using ForcedInterfacialWaves
using Printf

function generate_gravity_table()
    p = compute_gravity_parameters()
    
    println("=" ^ 100)
    println("PURE GRAVITY WAVE VALIDATION TABLE")
    println("=" ^ 100)
    println()
    println("Parameters:")
    @printf("  U = %.10f cm/s\n", p.U)
    @printf("  g = %.10f cm/s²\n", p.g)
    @printf("  ρ_l = %.10f g/cm³\n", p.rho_l)
    @printf("  ρ_u = %.10f g/cm³\n", p.rho_u)
    @printf("  ρ_r = %.16e\n", p.rho_r)
    @printf("  β = %.16e\n", p.beta)
    @printf("  √β = %.16e\n", p.sqrt_beta)
    @printf("  F0 = %.16e\n", p.F0)
    @printf("  λ_grav = %.16e\n", p.gravity_wavelength)
    @printf("  L = %.16e\n", p.L)
    println()

    # T0 values
    println("-" ^ 100)
    println("T0(x) — Steady contribution")
    println("-" ^ 100)
    @printf("  %-12s  %-24s\n", "x", "T0(x)")
    for x in [-5.0, -2.0, -1.0, 0.5, 1.0, 2.0, 5.0]
        val = gravity_T0(x, p)
        @printf("  %-12.4f  %+.16e\n", x, val)
    end
    println()

    # Analytical solution at multiple (x, t) points
    t_values = [2.0, 5.0, 10.0]
    
    println("-" ^ 100)
    println("Analytical solution: η_total, η_steady, η_transient")
    println("-" ^ 100)
    @printf("  %-8s  %-8s  %-24s  %-24s  %-24s\n", "x", "t", "η_total", "η_steady", "η_transient")

    for t in t_values
        # Left region points
        x_left = [-5.0, -2.0, -1.0, 0.5]
        for x in x_left
            a = t - x
            if a > p.front_band
                η_tot, η_st, η_tr = gravity_analytical_left(x, t, p)
                @printf("  %-8.2f  %-8.2f  %+.16e  %+.16e  %+.16e\n", x, t, η_tot, η_st, η_tr)
            end
        end

        # Right region points
        x_right = [t + 2.0, t + 4.0, t + 8.0]
        for x in x_right
            a = t - x
            if a < -p.front_band
                η_tot, η_st, η_tr = gravity_analytical_right(x, t, p)
                @printf("  %-8.2f  %-8.2f  %+.16e  %+.16e  %+.16e\n", x, t, η_tot, η_st, η_tr)
            end
        end
        println()
    end

    # CPV solution comparison
    println("-" ^ 100)
    println("Analytical vs CPV comparison")
    println("-" ^ 100)
    @printf("  %-8s  %-8s  %-24s  %-24s  %-12s\n", "x", "t", "η_analytical", "η_CPV", "|diff|")

    comparison_points = [
        (-5.0, 2.0), (-2.0, 2.0), (0.5, 2.0),
        (4.0, 2.0), (6.0, 2.0), (10.0, 2.0),
        (-3.0, 5.0), (1.0, 5.0), (7.0, 5.0), (13.0, 5.0),
    ]

    for (x, t) in comparison_points
        a = t - x
        if a > p.front_band
            η_an, _, _ = gravity_analytical_left(x, t, p)
        elseif a < -p.front_band
            η_an, _, _ = gravity_analytical_right(x, t, p)
        else
            continue
        end
        η_cpv = gravity_numerical_cpv(x, t, p)
        diff = abs(η_an - η_cpv)
        @printf("  %-8.2f  %-8.2f  %+.16e  %+.16e  %.4e\n", x, t, η_an, η_cpv, diff)
    end
    println()

    # Individual T terms at reference point
    println("-" ^ 100)
    println("Individual T terms at (x=-2, t=2)")
    println("-" ^ 100)
    x, t = -2.0, 2.0
    a = t - x
    @printf("  T0      = %+.16e\n", gravity_T0(x, p))
    @printf("  T1_left = %+.16e\n", gravity_T1_left(x, p))
    @printf("  T2_left = %+.16e\n", gravity_T2_left(x, t, a, p))
    @printf("  T3_left = %+.16e\n", gravity_T3_left(x, t, a, p))
    @printf("  T4_left = %+.16e\n", gravity_T4_left(x, t, a, p))
    println()

    println("-" ^ 100)
    println("Individual T terms at (x=4, t=2)")
    println("-" ^ 100)
    x, t = 4.0, 2.0
    a = t - x
    @printf("  T0       = %+.16e\n", gravity_T0(x, p))
    @printf("  T1_right = %+.16e\n", gravity_T1_right(x, p))
    @printf("  T2_right = %+.16e\n", gravity_T2_right(x, t, a, p))
    @printf("  T3_right = %+.16e\n", gravity_T3_right(x, t, a, p))
    @printf("  T4_right = %+.16e\n", gravity_T4_right(x, t, a, p))
    println()
end

function generate_cg_table()
    p = compute_cg_parameters()

    println("=" ^ 100)
    println("CAPILLARY-GRAVITY WAVE VALIDATION TABLE")
    println("=" ^ 100)
    println()
    println("Parameters:")
    @printf("  U = %.10f cm/s\n", p.U)
    @printf("  g = %.10f cm/s²\n", p.g)
    @printf("  T = %.10f dyn/cm\n", p.T)
    @printf("  ρ_l = %.10f g/cm³\n", p.rho_l)
    @printf("  ρ_u = %.10f g/cm³\n", p.rho_u)
    @printf("  α = %.16e\n", p.alpha)
    @printf("  ρ_r = %.16e\n", p.rho_r)
    @printf("  β = %.16e\n", p.beta)
    @printf("  γ_ρ = %.16e\n", p.gamma_rho)
    @printf("  k_l = %.16e\n", p.k_l)
    @printf("  k_s = %.16e\n", p.k_s)
    @printf("  F0 = %.16e\n", p.F0)
    @printf("  l_c = %.16e\n", p.l_c)
    println()

    # Dispersion relation values
    println("-" ^ 100)
    println("Dispersion Ω(k)")
    println("-" ^ 100)
    @printf("  %-12s  %-24s\n", "k", "Ω(k)")
    for k in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, p.k_s, p.k_l]
        @printf("  %-12.6f  %+.16e\n", k, dispersion_omega(k, p))
    end
    println()

    # Steady solution
    println("-" ^ 100)
    println("Steady capillary-gravity solution η_steady(x)")
    println("-" ^ 100)
    @printf("  %-12s  %-24s\n", "x", "η_steady")
    for x in [-10.0, -5.0, -2.0, -1.0, 1.0, 2.0, 5.0, 10.0]
        val = steady_surface_elevation(x, p)
        @printf("  %-12.4f  %+.16e\n", x, val)
    end
    println()

    # IVP solution at various times
    println("-" ^ 100)
    println("IVP vs Steady comparison at x = 2.0")
    println("-" ^ 100)
    @printf("  %-8s  %-24s  %-24s  %-12s\n", "t", "η_IVP", "η_steady", "|diff|")
    x_ref = 2.0
    η_st = steady_surface_elevation(x_ref, p)
    for t in [5.0, 10.0, 20.0, 50.0, 100.0, 200.0]
        η_ivp = ivp_surface_elevation(x_ref, t, p)
        diff = abs(η_ivp - η_st)
        @printf("  %-8.1f  %+.16e  %+.16e  %.4e\n", t, η_ivp, η_st, diff)
    end
    println()

    # Full spatial profiles at t = 50
    println("-" ^ 100)
    println("Full IVP profile at t = 50")
    println("-" ^ 100)
    @printf("  %-12s  %-24s  %-24s  %-12s\n", "x", "η_IVP(x,50)", "η_steady(x)", "|diff|")
    for x in [-10.0, -5.0, -2.0, -1.0, 1.0, 2.0, 5.0, 10.0]
        η_ivp = ivp_surface_elevation(x, 50.0, p)
        η_st_x = steady_surface_elevation(x, p)
        diff = abs(η_ivp - η_st_x)
        @printf("  %-12.4f  %+.16e  %+.16e  %.4e\n", x, η_ivp, η_st_x, diff)
    end
    println()
end

function generate_fresnel_table()
    println("=" ^ 100)
    println("FRESNEL INTEGRALS REFERENCE VALUES")
    println("=" ^ 100)
    println()
    @printf("  %-12s  %-24s  %-24s\n", "x", "C(x)", "S(x)")
    for x in [0.0, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0]
        @printf("  %-12.4f  %+.16e  %+.16e\n", x, fresnel_C(x), fresnel_S(x))
    end
    println()
end

# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
println("\n" * "═" ^ 100)
println("   FORCEDINTERFACIALWAVES.JL VALIDATION TABLE — FOR MATLAB-JULIA COMPARISON")
println("═" ^ 100 * "\n")

generate_fresnel_table()
generate_gravity_table()
generate_cg_table()

println("\n" * "═" ^ 100)
println("   END OF VALIDATION TABLE")
println("═" ^ 100)
