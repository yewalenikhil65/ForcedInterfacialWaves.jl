# Validation

Numerical validation demonstrating MATLAB–Julia equivalence at the integral level.

## Parameters (15-digit agreement)

### Pure Gravity

| Symbol | Value |
|:-------|------:|
| $\beta$ | `9.980019980019981e-01` |
| $\sqrt\beta$ | `9.990004995003746e-01` |
| $F_0$ | `1.388855922278772e-03` |
| $\lambda_g$ | `6.295764256743509e+00` |

### Capillary–Gravity

| Symbol | Value |
|:-------|------:|
| $\alpha$ | `1.388855922278772e-01` |
| $\rho_r$ | `1.000000000000000e-03` |
| $\beta$ | `9.980019980019981e-01` |
| $\gamma_\rho$ | `9.990009990009991e-01` |
| $k_s$ | `1.196700136678134e+00` |
| $k_l$ | `6.010670937188218e+00` |
| $F_0$ | `1.388855922278772e-03` |
| $l_c$ | `7.269476668297654e-01` cm |
| $t_c$ | `2.722181447502548e-02` s |

---

## Pure-Gravity Individual $T$-Terms ($x = -2$, $t_{\dim} = 1$ s)

| Term | Value | Method |
|:-----|------:|:-------|
| $T_0$ | $-8.640 \times 10^{-1}$ | `quadgk` (Laplace-type) |
| $T_1^-$ | $+9.100 \times 10^{-1}$ | closed-form |
| $T_2^-$ | $+5.831 \times 10^{-3}$ | `quadgk` (exponentially damped) |
| $T_3^-$ | $+7.045 \times 10^{-4}$ | `fresnelc` / `fresnels` (no quadrature) |
| $T_4^-$ | $-7.003 \times 10^{-4}$ | `quadgk` (oscillatory) |

---

## Analytical vs Numerical CPV ($t_{\dim} = 1$ s)

### Left region

| $x$ | $\eta_{\mathrm{analytical}}$ | $\eta_{\mathrm{CPV}}$ | $|\Delta|$ |
|----:|---:|---:|---:|
| $-5$ | $+4.41 \times 10^{-5}$ | $+5.27 \times 10^{-5}$ | $8.6 \times 10^{-6}$ |
| $-2$ | $+7.21 \times 10^{-5}$ | $+6.98 \times 10^{-5}$ | $2.3 \times 10^{-6}$ |
| $0.5$ | $-1.06 \times 10^{-3}$ | $-1.07 \times 10^{-3}$ | $1.3 \times 10^{-5}$ |

### Right region

| $x$ | $\eta_{\mathrm{analytical}}$ | $\eta_{\mathrm{CPV}}$ | $|\Delta|$ |
|----:|---:|---:|---:|
| $t+5$ | $+8.47 \times 10^{-4}$ | $+7.94 \times 10^{-4}$ | $5.2 \times 10^{-5}$ |
| $t+8$ | $+3.97 \times 10^{-4}$ | $+4.38 \times 10^{-4}$ | $4.1 \times 10^{-5}$ |

---

## Capillary–Gravity IVP vs Steady ($t_{\dim} = 3$ s)

| $x$ | $\eta_{\mathrm{IVP}}$ | $\eta_{\mathrm{steady}}$ | $|\Delta|$ |
|----:|---:|---:|---:|
| $-3$ | $-3.288 \times 10^{-3}$ | $-2.993 \times 10^{-3}$ | $2.9 \times 10^{-4}$ |
| $-1$ | $-1.325 \times 10^{-3}$ | $-9.493 \times 10^{-4}$ | $3.8 \times 10^{-4}$ |
| $1$ | $-3.868 \times 10^{-3}$ | $-3.698 \times 10^{-3}$ | $1.7 \times 10^{-4}$ |
| $3$ | $+1.884 \times 10^{-3}$ | $+1.839 \times 10^{-3}$ | $4.5 \times 10^{-5}$ |
| $5$ | $+1.412 \times 10^{-3}$ | $+1.242 \times 10^{-3}$ | $1.7 \times 10^{-4}$ |

---

## Optimized Profile Equivalence

The production Julia profile algorithm uses vectorised in-place quadrature (`quadgk!`) over threaded spatial chunks. It computes the same integrals as the scalar MATLAB-style reference:

| Metric | Value |
|:-------|------:|
| Max pointwise $|\eta_{\mathrm{vector}} - \eta_{\mathrm{scalar}}|$ | $2.58 \times 10^{-10}$ |
| Max pointwise $|\eta_{\mathrm{steady,vector}} - \eta_{\mathrm{steady,scalar}}|$ | $5.73 \times 10^{-10}$ |
| Speedup (8 threads, 201 points) | 9.73× |

---

## Reproducing

```julia
using WavesDelta, Printf

# PG analytical vs CPV
pg = compute_gravity_parameters()
t_pg = 1.0 / pg.t_c
η_a, _, _ = gravity_analytical_left(-2.0, t_pg, pg)
η_c = gravity_numerical_cpv(-2.0, t_pg, pg)
@printf("|analytical − CPV| at x=−2: %.4e\n", abs(η_a - η_c))

# CG IVP
p = compute_cg_parameters()
t = 3.0 / p.t_c
@printf("|IVP − steady| at x=3: %.4e\n",
    abs(ivp_surface_elevation(3.0, t, p) - steady_surface_elevation(3.0, p)))
```

Run the full integral-level validation script:

```bash
julia --project=. scripts/validate_integrals.jl
```

Run the test suite:

```bash
julia -t auto --project=. test/runtests.jl
```
