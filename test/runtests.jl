"""
    runtests.jl — Validation test suite for ForcedInterfacialWaves.jl

    Validates Julia implementation against MATLAB reference formulas
    from jfm_matlab_codes.tex.
"""

push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using ForcedInterfacialWaves
using Test
using Printf

@testset "ForcedInterfacialWaves.jl" begin

    # ═══════════════════════════════════════════════════════════════════════════
    @testset "Fresnel integrals" begin
        # Known values (NIST / MATLAB verified)
        @test fresnel_C(0.0) ≈ 0.0 atol=1e-15
        @test fresnel_S(0.0) ≈ 0.0 atol=1e-15
        @test fresnel_C(1.0) ≈ 0.7798934003768228 atol=1e-12
        @test fresnel_S(1.0) ≈ 0.4382591473903547 atol=1e-12
        @test fresnel_C(1000.0) ≈ 0.5 atol=1e-3  # asymptotic
        @test fresnel_S(1000.0) ≈ 0.5 atol=1e-3
        # Odd function
        @test fresnel_C(-1.0) ≈ -fresnel_C(1.0) atol=1e-14
        @test fresnel_S(-1.0) ≈ -fresnel_S(1.0) atol=1e-14
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testset "Capillary-gravity parameters" begin
        p = compute_cg_parameters()

        l_c = 26.7046^2 / 981.0
        @test p.l_c ≈ l_c atol=1e-12
        @test p.t_c ≈ 26.7046 / 981.0 atol=1e-12

        alpha_ref = 72.0 / (1.0 * 26.7046^2 * l_c)
        @test p.alpha ≈ alpha_ref atol=1e-14
        @test p.rho_r ≈ 0.001 atol=1e-16
        @test p.beta ≈ (1.0 - 0.001) / (1.0 + 0.001) atol=1e-14
        @test p.gamma_rho ≈ 1.0 / (1.0 + 0.001) atol=1e-14

        # Roots satisfy quadratic
        disc = (1.001)^2 - 4.0 * p.alpha * 0.999
        @test p.k_l ≈ (1.001 + sqrt(disc)) / (2.0 * p.alpha) atol=1e-12
        @test p.k_s ≈ (1.001 - sqrt(disc)) / (2.0 * p.alpha) atol=1e-12
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testset "CG dispersion relation consistency" begin
        p = compute_cg_parameters()
        # At roots, the steady condition should hold:
        # χ(k)² = β k + γ_ρ α k³
        for k in [p.k_s, p.k_l]
            χ² = dispersion_chi(k, p)^2
            @test χ² ≈ p.beta * k + p.gamma_rho * p.alpha * k^3 atol=1e-12
        end
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testset "CG IVP converges to steady at large t" begin
        p = compute_cg_parameters()
        t = 3.0 / p.t_c  # t_dim = 3.0 s (time_index = 300)

        for x in [-3.0, 1.0, 3.0, 5.0]
            η_ivp = ivp_surface_elevation(x, t, p)
            η_std = steady_surface_elevation(x, p)
            @test isfinite(η_ivp)
            @test isfinite(η_std)
            @test abs(η_ivp - η_std) < 5e-4  # should be small at t=3s
        end
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testset "Pure-gravity parameters" begin
        pg = compute_gravity_parameters()
        @test pg.beta ≈ 0.999 / 1.001 atol=1e-14
        @test pg.sqrt_beta ≈ sqrt(0.999 / 1.001) atol=1e-14
        @test pg.gravity_wavelength ≈ 2π / pg.beta atol=1e-12
        @test pg.L ≈ 4.0 * pg.gravity_wavelength atol=1e-12
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testset "PG T0 symmetry" begin
        pg = compute_gravity_parameters()
        # T0 depends on |x|, so T0(x) == T0(-x)
        @test gravity_T0(3.0, pg) ≈ gravity_T0(-3.0, pg) atol=1e-12
        @test gravity_T0(1.0, pg) ≈ gravity_T0(-1.0, pg) atol=1e-12
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testset "PG analytical vs CPV agreement (left region)" begin
        pg = compute_gravity_parameters()
        t = 1.0 / pg.t_c

        for x in [-5.0, -2.0, 0.5]
            a = t - x
            @assert a > pg.front_band
            η_ana, _, _ = gravity_analytical_left(x, t, pg)
            η_cpv = gravity_numerical_cpv(x, t, pg)
            diff = abs(η_ana - η_cpv)
            @test diff < 1e-4
            @printf("    x=%6.2f  |diff| = %.3e\n", x, diff)
        end
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testset "PG analytical vs CPV agreement (right region)" begin
        pg = compute_gravity_parameters()
        t = 1.0 / pg.t_c

        for x in [t + 3.0, t + 8.0]
            a = t - x
            @assert a < -pg.front_band
            η_ana, _, _ = gravity_analytical_right(x, t, pg)
            η_cpv = gravity_numerical_cpv(x, t, pg)
            diff = abs(η_ana - η_cpv)
            @test diff < 2e-4
            @printf("    x=%6.2f  |diff| = %.3e\n", x, diff)
        end
    end

    # ═══════════════════════════════════════════════════════════════════════════
    @testset "PG individual terms are finite" begin
        pg = compute_gravity_parameters()
        t = 1.0 / pg.t_c
        x = -2.0
        a = t - x

        @test isfinite(gravity_T0(x, pg))
        @test isfinite(gravity_T1_left(x, pg))
        @test isfinite(gravity_T2_left(x, t, a, pg))
        @test isfinite(gravity_T3_left(x, t, a, pg))
        @test isfinite(gravity_T4_left(x, t, a, pg))

        xr = t + 5.0
        ar = t - xr
        @test isfinite(gravity_T1_right(xr, pg))
        @test isfinite(gravity_T2_right(xr, t, ar, pg))
        @test isfinite(gravity_T3_right(xr, t, ar, pg))
        @test isfinite(gravity_T4_right(xr, t, ar, pg))
    end

end  # top-level testset


    # ═══════════════════════════════════════════════════════════════════════════
    @testset "CG optimized profile equivalence" begin
        # Finite tail keeps this test fast; production defaults remain k_max=Inf.
        p = compute_cg_parameters(k_max=50.0, k_max_steady=50.0)
        t = 0.15 / p.t_c
        x_grid = collect(range(-5.0, 10.0; length=33))

        values = zeros(length(x_grid))
        for k in (0.5, 3.0, 10.0)
            cg_profile_integrand!(values, k, x_grid, t, p)
            scalar_values = [cg_combined_integrand(k, x, t, p) for x in x_grid]
            @test maximum(abs.(values .- scalar_values)) < 5e-12
        end

        eta_vector = compute_cg_ivp_profile(x_grid, t, p;
                                             method=:threaded_vector)
        eta_scalar = compute_cg_ivp_profile(x_grid, t, p;
                                             method=:threaded_scalar)
        @test maximum(abs.(eta_vector .- eta_scalar)) < 1e-9

        eta_steady_vector = compute_cg_steady_profile(x_grid, p)
        eta_steady_scalar = [steady_surface_elevation(x, p) for x in x_grid]
        @test maximum(abs.(eta_steady_vector .- eta_steady_scalar)) < 2e-9

        eta_only, steady_skipped = compute_cg_profile(
            x_grid, t, p; compute_steady=false)
        @test steady_skipped === nothing
        @test eta_only ≈ eta_vector atol=1e-12 rtol=1e-10
    end
println("\n✓ All tests passed.")
