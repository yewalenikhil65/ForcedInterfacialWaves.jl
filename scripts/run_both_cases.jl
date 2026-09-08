"""
    Main driver script for ForcedInterfacialWaves.jl
    Computes profiles for both capillary-gravity (Figure 10) and pure-gravity (Figure 6) cases.

    Usage:
        julia --project=. scripts/run_both_cases.jl
"""

using Printf

push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using ForcedInterfacialWaves

# ═══════════════════════════════════════════════════════════════════════════════
# Capillary-Gravity Case (Figure 10)
# ═══════════════════════════════════════════════════════════════════════════════

function run_capillary_gravity()
    println("="^70)
    println("  CAPILLARY-GRAVITY CASE (Figure 10, α > 0)")
    println("="^70)

    p = compute_cg_parameters()

    @printf("  Parameters:\n")
    @printf("    α     = %.6e\n", p.alpha)
    @printf("    ρ_r   = %.6e\n", p.rho_r)
    @printf("    β     = %.15f\n", p.beta)
    @printf("    γ_ρ   = %.15f\n", p.gamma_rho)
    @printf("    k_s   = %.15f\n", p.k_s)
    @printf("    k_l   = %.15f\n", p.k_l)
    @printf("    F0    = %.6e\n", p.F0)
    @printf("    l_c   = %.6f cm\n", p.l_c)
    @printf("    t_c   = %.6f s\n", p.t_c)

    # Dimensional wavelengths
    k_l_dim = p.k_l / p.l_c
    k_s_dim = p.k_s / p.l_c
    lambda_c = 2π / k_l_dim
    lambda_g = 2π / k_s_dim
    @printf("    λ_capillary = %.4f cm\n", lambda_c)
    @printf("    λ_gravity   = %.4f cm\n", lambda_g)
    println()

    # Time indices to compute (subset for demonstration)
    time_indices = [1, 15, 60, 300]

    # Spatial grid (reduced for speed)
    Nx = 201
    x_grid = make_cg_xgrid(p; Nx=Nx)

    for ti in time_indices
        t_dim = ti / 100.0
        t = t_dim / p.t_c

        @printf("  Computing t_dim = %.2f s (t = %.4f, time_index = %d)...\n", t_dim, t, ti)

        η_ivp, η_steady = compute_cg_profile(x_grid, t, p)

        # Print a few sample values
        mid = div(length(x_grid), 2)
        sample_indices = [1, div(mid, 2), mid, mid + div(mid, 2), length(x_grid)]

        @printf("    %10s  %16s  %16s\n", "x", "1e3*η_IVP", "1e3*η_steady")
        for ix in sample_indices
            @printf("    %10.4f  %16.8e  %16.8e\n",
                    x_grid[ix], 1e3 * η_ivp[ix], 1e3 * η_steady[ix])
        end
        println()
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Pure-Gravity Case (Figure 6)
# ═══════════════════════════════════════════════════════════════════════════════

function run_pure_gravity()
    println("="^70)
    println("  PURE-GRAVITY CASE (Figure 6, α = 0)")
    println("="^70)

    p = compute_gravity_parameters()

    @printf("  Parameters:\n")
    @printf("    β         = %.15f\n", p.beta)
    @printf("    √β        = %.15f\n", p.sqrt_beta)
    @printf("    F0        = %.6e\n", p.F0)
    @printf("    λ_gravity = %.6f\n", p.gravity_wavelength)
    @printf("    L         = %.6f\n", p.L)
    @printf("    l_c       = %.6f cm\n", p.l_c)
    @printf("    t_c       = %.6f s\n", p.t_c)
    println()

    # Time indices
    time_indices = [1, 15, 60, 100]

    # Spatial grid (reduced for speed)
    Nx = 201
    x_grid = make_gravity_xgrid(p; Nx=Nx)

    for ti in time_indices
        t_dim = ti / 100.0
        t = t_dim / p.t_c

        @printf("  Computing t_dim = %.2f s (t = %.4f, time_index = %d)...\n", t_dim, t, ti)

        η_analytical, η_steady, η_transient, η_cpv = compute_gravity_profile(x_grid, t, p)

        # Compute agreement between analytical and numerical CPV
        valid = .!isnan.(η_analytical) .& .!isnan.(η_cpv)
        if any(valid)
            abs_errors = abs.(η_analytical[valid] .- η_cpv[valid])
            max_abs_err = maximum(abs_errors)
            mean_abs_err = sum(abs_errors) / length(abs_errors)
            @printf("    Max |η_A - η_CPV| = %.4e\n", max_abs_err)
            @printf("    Mean|η_A - η_CPV| = %.4e\n", mean_abs_err)
        end

        # Print a few sample values
        left_idx = findall(valid .& (x_grid .< t))
        right_idx = findall(valid .& (x_grid .> t))

        @printf("    Sample (left of front):\n")
        @printf("    %10s  %16s  %16s  %16s\n", "x", "1e3*η_anal", "1e3*η_CPV", "1e3*η_steady")
        for ix in left_idx[1:min(3, length(left_idx))]
            @printf("    %10.4f  %16.8e  %16.8e  %16.8e\n",
                    x_grid[ix], 1e3 * η_analytical[ix], 1e3 * η_cpv[ix], 1e3 * η_steady[ix])
        end

        if !isempty(right_idx)
            @printf("    Sample (right of front):\n")
            for ix in right_idx[1:min(3, length(right_idx))]
                @printf("    %10.4f  %16.8e  %16.8e  %16.8e\n",
                        x_grid[ix], 1e3 * η_analytical[ix], 1e3 * η_cpv[ix], 1e3 * η_steady[ix])
            end
        end
        println()
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Run both cases
# ═══════════════════════════════════════════════════════════════════════════════

function main()
    println("\n" * "╔" * "═"^68 * "╗")
    println("║" * " "^10 * "ForcedInterfacialWaves.jl — Interfacial Waves IVP" * " "^20 * "║")
    println("╚" * "═"^68 * "╝\n")

    t_start = time()

    run_capillary_gravity()
    println()
    run_pure_gravity()

    elapsed = time() - t_start
    @printf("\n  Total elapsed time: %.2f s\n", elapsed)
    println("\n  Done.")
end

main()
