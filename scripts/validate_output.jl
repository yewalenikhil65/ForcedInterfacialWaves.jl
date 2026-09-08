## validate_output.jl — Run in Julia, paste output back
#  Produces numerical values for all key integrals for comparison with MATLAB.
#  Usage: julia --project=. scripts/validate_output.jl

push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using ForcedInterfacialWaves, Printf

## Parameters
p = compute_cg_parameters()
println("=== CG PARAMETERS ===")
@printf("alpha   = %.15e\n", p.alpha)
@printf("k_s     = %.15e\n", p.k_s)
@printf("k_l     = %.15e\n", p.k_l)
@printf("F0      = %.15e\n", p.F0)

## chi(2)
@printf("\nchi(2)  = %.15e\n", dispersion_chi(2.0, p))

## Combined integrand at k=2, x=3, t=110
@printf("\nI(2; x=3, t=110) = %.15e\n", cg_combined_integrand(2.0, 3.0, 110.0, p))

## Partial integrals
I1, I2, I3 = cg_partial_integrals(3.0, 110.0, p)
η_ivp = ivp_surface_elevation(3.0, 110.0, p)
println("\n=== CG PARTIAL INTEGRALS (x=3, t=110) ===")
@printf("I1      = %.15e\n", I1)
@printf("I2      = %.15e\n", I2)
@printf("I3      = %.15e\n", I3)
@printf("eta_IVP = %.15e\n", η_ivp)

## Steady
Gx = cg_Gx_integral(3.0, p)
η_s = steady_surface_elevation(3.0, p)
@printf("G(3)        = %.15e\n", Gx)
@printf("eta_steady  = %.15e\n", η_s)

## ========== PURE GRAVITY ==========
pg = compute_gravity_parameters()
println("\n=== PG PARAMETERS ===")
@printf("beta       = %.15e\n", pg.beta)
@printf("sqrt_beta  = %.15e\n", pg.sqrt_beta)

## T terms at x=-2, t_dim=1.0s
t_pg = 1.0 / pg.t_c
x_pg = -2.0
a_pg = t_pg - x_pg
@printf("\n=== PG T-TERMS (x=-2, t_dim=1s, t=%.6f, a=%.6f) ===\n", t_pg, a_pg)

@printf("T0  = %.15e\n", gravity_T0(x_pg, pg))
@printf("T1- = %.15e\n", gravity_T1_left(x_pg, pg))
@printf("T2- = %.15e\n", gravity_T2_left(x_pg, t_pg, a_pg, pg))
@printf("T3- = %.15e\n", gravity_T3_left(x_pg, t_pg, a_pg, pg))
@printf("T4- = %.15e\n", gravity_T4_left(x_pg, t_pg, a_pg, pg))

η_analytical, _, _ = gravity_analytical_left(x_pg, t_pg, pg)
@printf("eta_analytical = %.15e\n", η_analytical)

## CPV
η_cpv = gravity_numerical_cpv(x_pg, t_pg, pg)
@printf("eta_CPV        = %.15e\n", η_cpv)
@printf("|diff|          = %.4e\n", abs(η_analytical - η_cpv))

println("\nDone.")
