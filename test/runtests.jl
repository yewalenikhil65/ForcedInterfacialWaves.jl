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
    @testset "CG spatial grid" begin
        p = compute_cg_parameters()
        x_grid = make_cg_xgrid(p; Nx=401, xlim=(-10.0, 10.0))

        @test first(x_grid) ≈ -10.0
        @test last(x_grid) ≈ 10.0
        @test all(abs.(x_grid) .> 1.0e-12)
        @test length(x_grid) == 400  # x ≈ 0 is intentionally excluded
        @test_throws ArgumentError make_cg_xgrid(p; xlim=(10.0, -10.0))
        @test_throws ArgumentError make_cg_xgrid(p; xlim=(-10.0,))
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

        eta_symmetric_vector = compute_cg_steady_profile(
            x_grid, p; rayleigh_dissipation=false)
        eta_symmetric_scalar = [steady_surface_elevation(
            x, p; rayleigh_dissipation=false) for x in x_grid]
        @test maximum(abs.(eta_symmetric_vector .- eta_symmetric_scalar)) < 2e-9

        eta_only, steady_skipped = compute_cg_profile(
            x_grid, t, p; compute_steady=false)
        @test steady_skipped === nothing
        @test eta_only ≈ eta_vector atol=1e-12 rtol=1e-10

        eta_profile, eta_profile_s = compute_cg_profile(x_grid, t, p)
        @test eta_profile ≈ eta_vector atol=1e-12 rtol=1e-10
        @test eta_profile_s ≈ eta_symmetric_vector

        _, eta_profile_classical = compute_cg_profile(
            x_grid, t, p; rayleigh_dissipation=true)
        @test eta_profile_classical ≈ eta_steady_vector
    end
    @testset "CG I3/I4 component APIs" begin
        p = compute_cg_parameters(k_max=30.0, k_max_steady=30.0)
        t = 0.15 / p.t_c
        x_grid = [-2.0, 1.0, 3.0]

        @test 𝕀₃(3.0, t, p) ≈ cg_I3(3.0, t, p)
        @test 𝕀₄(3.0, t, p) ≈ cg_I4(3.0, t, p)

        i3_scalar = [cg_I3(x, t, p) for x in x_grid]
        i4_scalar = [cg_I4(x, t, p) for x in x_grid]
        i3_vector = compute_cg_I3_profile(x_grid, t, p;
                                           method=:vector)
        i4_vector = compute_cg_I4_profile(x_grid, t, p;
                                           method=:vector)
        @test i3_vector ≈ i3_scalar atol=2e-8 rtol=2e-8
        @test i4_vector ≈ i4_scalar atol=2e-8 rtol=2e-8
        @test compute_cg_I3_profile(x_grid, t, p;
                                    method=:threaded_vector) ≈ i3_scalar atol=2e-8 rtol=2e-8
        @test compute_cg_I4_profile(x_grid, t, p;
                                    method=:threaded_vector) ≈ i4_scalar atol=2e-8 rtol=2e-8
    end
    @testset "solve method dispatches" begin
        pg = compute_gravity_parameters(k_max_analytical=20.0, k_max_cpv=20.0)
        x = -2.0
        t = 1.0 / pg.t_c
        prob_pg = ForcedGravityProblem(pg, x, t)

        default_pg = solve(prob_pg)
        explicit_pg = solve(prob_pg; method=IVP())
        legacy_pg = solve(prob_pg; method=:IVP)
        @test explicit_pg.η ≈ default_pg.η
        @test legacy_pg.η ≈ default_pg.η
        @test explicit_pg.η_steady ≈ default_pg.η_steady
        @test length((explicit_pg.η_s_local, explicit_pg.η_s_farfield)) == 2

        steady_pg_no_rayleigh = solve(prob_pg;
                                      method=steady(rayleigh_dissipation=false))
        steady_pg_rayleigh = solve(prob_pg;
                                   method=steady(rayleigh_dissipation=true))
        @test steady_pg_no_rayleigh.η ≈ pg.F0 * T₀(x, pg)
        @test steady_pg_no_rayleigh.η == steady_pg_no_rayleigh.η_steady
        @test steady_pg_no_rayleigh.η_transient == 0.0
        @test steady_pg_no_rayleigh.η_s_local + steady_pg_no_rayleigh.η_s_farfield ≈ steady_pg_no_rayleigh.η
        @test steady_pg_rayleigh.η_s_local ≈ steady_pg_no_rayleigh.η_s_local
        @test steady_pg_rayleigh.η_s_farfield != steady_pg_no_rayleigh.η_s_farfield
        @test solve(prob_pg; method=:steady).η ≈ steady_pg_rayleigh.η

        x_grid_pg = [-2.0, 2.0]
        steady_pg_profile = solve(ForcedGravityProblem(pg, x_grid_pg);
                                  method=steady(rayleigh_dissipation=false))
        @test steady_pg_profile.t === nothing
        @test steady_pg_profile.η ≈ pg.F0 .* T₀.(x_grid_pg, Ref(pg))
        @test steady_pg_profile.η == steady_pg_profile.η_steady
        @test all(iszero, steady_pg_profile.η_transient)
        @test steady_pg_profile.η_s_local .+ steady_pg_profile.η_s_farfield ≈ steady_pg_profile.η
        @test_throws ArgumentError solve(ForcedGravityProblem(pg, x_grid_pg))

        p = compute_cg_parameters(k_max=20.0, k_max_steady=20.0)
        prob_cg = ForcedGCProblem(p, 3.0, 0.15 / p.t_c)
        default_cg = solve(prob_cg)
        explicit_cg = solve(prob_cg; method=IVP())
        symmetric_cg = solve(prob_cg; method=IVP(asym_cancel=false))
        asymmetric_cg = solve(prob_cg; method=IVP(; asym_cancel=true))
        @test IVP().asym_cancel == false
        @test IVP(; asym_cancel=true).asym_cancel == true
        @test explicit_cg.η ≈ default_cg.η
        @test symmetric_cg.η ≈ default_cg.η
        @test asymmetric_cg.η ≈ default_cg.η

        steady_cg_no_rayleigh = solve(prob_cg;
                                      method=steady(rayleigh_dissipation=false))
        steady_cg_rayleigh = solve(prob_cg;
                                   method=steady(rayleigh_dissipation=true))
        @test explicit_cg.η_steady ≈ default_cg.η_steady
        @test symmetric_cg.η_steady ≈ steady_cg_no_rayleigh.η
        @test asymmetric_cg.η_steady ≈ steady_cg_rayleigh.η
        @test symmetric_cg.η_transient ≈ symmetric_cg.η - symmetric_cg.η_steady
        @test asymmetric_cg.η_transient ≈ asymmetric_cg.η - asymmetric_cg.η_steady
        @test explicit_cg.η_s_local + explicit_cg.η_s_farfield ≈ explicit_cg.η_steady
        @test explicit_cg.η_transient ≈ explicit_cg.η - explicit_cg.η_steady
        @test steady_cg_no_rayleigh.η ≈ steady_surface_elevation(
            3.0, p; rayleigh_dissipation=false)
        @test steady_cg_rayleigh.η ≈ steady_surface_elevation(
            3.0, p; rayleigh_dissipation=true)
        @test steady_cg_no_rayleigh.η_s_local ≈ steady_cg_rayleigh.η_s_local
        @test steady_cg_no_rayleigh.η_s_farfield != steady_cg_rayleigh.η_s_farfield
        @test solve(prob_cg; method=:steady).η ≈ steady_cg_rayleigh.η

        x_grid_cg = [-2.0, 3.0]
        steady_cg_profile = solve(ForcedGCProblem(p, x_grid_cg);
                                  method=steady(rayleigh_dissipation=true))
        @test steady_cg_profile.t === nothing
        @test steady_cg_profile.η ≈ compute_cg_steady_profile(x_grid_cg, p)
        @test steady_cg_profile.η == steady_cg_profile.η_steady
        @test all(iszero, steady_cg_profile.η_transient)
        @test steady_cg_profile.η_s_local .+ steady_cg_profile.η_s_farfield ≈ steady_cg_profile.η

        ivp_cg_profile = solve(ForcedGCProblem(p, x_grid_cg, 0.15 / p.t_c))
        ivp_cg_profile_symmetric = solve(
            ForcedGCProblem(p, x_grid_cg, 0.15 / p.t_c);
            method=IVP(asym_cancel=false))
        ivp_cg_profile_asymmetric = solve(
            ForcedGCProblem(p, x_grid_cg, 0.15 / p.t_c);
            method=IVP(asym_cancel=true))
        eta_s_profile = compute_cg_steady_profile(
            x_grid_cg, p; rayleigh_dissipation=false)
        eta_classical_profile = compute_cg_steady_profile(
            x_grid_cg, p; rayleigh_dissipation=true)
        @test ivp_cg_profile.η_steady ≈ eta_s_profile
        @test ivp_cg_profile_symmetric.η_steady ≈ eta_s_profile
        @test ivp_cg_profile_asymmetric.η_steady ≈ eta_classical_profile
        @test ivp_cg_profile.η_transient ≈
              ivp_cg_profile.η .- ivp_cg_profile.η_steady
        @test ivp_cg_profile_asymmetric.η_transient ≈
              ivp_cg_profile_asymmetric.η .- ivp_cg_profile_asymmetric.η_steady
        @test ivp_cg_profile.η_s_local .+ ivp_cg_profile.η_s_farfield ≈
              ivp_cg_profile.η_steady
        @test ivp_cg_profile_asymmetric.η_s_local .+
              ivp_cg_profile_asymmetric.η_s_farfield ≈
              ivp_cg_profile_asymmetric.η_steady

        @test_throws ArgumentError solve(ForcedGCProblem(p, x_grid_cg))

        @test_throws ArgumentError solve(prob_pg; method=:unknown)
        @test_throws ArgumentError solve(prob_cg; method=:unknown)
    end

println("\n✓ All tests passed.")
