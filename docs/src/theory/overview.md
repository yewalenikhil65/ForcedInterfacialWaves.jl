# Theory

This section contains the mathematical formulation and numerical implementation of the forced interfacial wave IVP, following the notation of [Kadari et al. (2026)](https://arxiv.org/abs/2605.12254).

## Problem description and nondimensionalisation

```@raw html
<figure style="text-align:center;">
  <img src="../assets/Fig3.png" alt="Pressure forcing at a two-fluid interface" style="max-width:100%;">
</figure>
```

A localised pressure forcing $p_e = F_0\,\delta(x)$ acts at the interface of two inviscid, incompressible, irrotational fluids of infinite depth. Both streams move at uniform speed $U$ rightwards. In the absence of forcing, the interface is flat at $z = 0$.

After nondimensionalisation using

```math
l_c = \frac{U^2}{g}, \quad t_c = \frac{U}{g}, \quad p_c = \rho_l U^2,
```

the key dimensionless groups are

```math
\alpha = \frac{g T}{\rho_l U^4}, \quad
\rho_r = \frac{\rho_u}{\rho_l}, \quad
\beta = \frac{1-\rho_r}{1+\rho_r}.
```

The parameter $\alpha$ measures the relative importance of surface tension. Two limiting cases are studied:

| Regime | Condition | Dispersion | Singularity structure |
|:-------|:---------:|:-----------|:----------------------|
| Pure gravity | $\alpha = 0$ | $\chi(k) = \sqrt{\beta k}$ | Single CPV pole at $k = \beta$ |
| Capillary–gravity | $\alpha > 0$ | $\chi(k) = \sqrt{\beta k + \frac{\alpha}{1+\rho_r} k^3}$ | Removable poles at $k_s$, $k_l$ (cancel in combined integrand) |

## Structure of the solution

In both regimes, the interfacial displacement decomposes as

```math
\eta(x,t) = \eta_s(x) + \eta_{\mathrm{tr}}(x,t),
```

where $\eta_s$ is a nominally time-independent (steady) contribution and $\eta_{\mathrm{tr}}$ contains the transient wave dynamics. A central result of the manuscript is that the *transient* terms also contribute to the final steady state through asymmetric cancellations as $t \to \infty$.

## Pages in this section

- **[Steady-State](steady_state.md)** — Manuscript §3.4: the general Fourier integral solution (eqns. 3.7–3.8), the steady far-field and local decomposition (eqn. 3.11), the Rayleigh dissipation form (eqn. 3.12), and numerical evaluation with Julia/MATLAB code.

- **[Pure Gravity](pure_gravity.md)** — Manuscript §4.1: setting $\alpha = 0$, the analytical $T_0$–$T_4^\pm$ decomposition, Fresnel-integral representations, independent CPV verification, and the spatial IVP profile at $t = 183.68$ (manuscript Fig. 6). Includes Julia/MATLAB implementations.

- **[Capillary–Gravity](capillary_gravity.md)** — Manuscript §4.2: the finite-capillarity case with the combined integrand technique, steady-state asymmetric cancellation proof, the $\mathbb{I}_4$ transient decay (Fig. 7), the full IVP profile (Fig. 8), and comparison with nonlinear simulations (Fig. 10). Includes Julia/MATLAB implementations.

Each sub-page is self-contained: it introduces the relevant equations, presents the derivation, and provides working code.

### Manuscript figure index

| Figure | Content | Theory page | Notebook cell |
|:------:|:--------|:------------|:-------------|
| 5a | Steady CG: far-field + local (no Rayleigh) | [Steady-State](steady_state.md) | 9 |
| 5b | Steady CG: far-field + local (Rayleigh dissipation) | [Steady-State](steady_state.md) | 10 |
| 6 | Pure-gravity IVP at $t = 183.68$ | [Pure Gravity](pure_gravity.md) | 13 |
| 7 | $\mathbb{I}_4(x,t)$ transient decay at $t = 0.34$ | [Capillary–Gravity](capillary_gravity.md) | 17 |
| 8 | CG full IVP at $t = 367.35$ | [Capillary–Gravity](capillary_gravity.md) | 19 |
| 10 | CG IVP vs nonlinear simulation at $t_{\dim} = 25$ s | [Capillary–Gravity](capillary_gravity.md) | 22 |

The companion [`forced_interfacial_waves_usage.ipynb`](https://github.com/yewalenikhil65/ForcedInterfacialWaves.jl/blob/main/notebooks/forced_interfacial_waves_usage.ipynb) notebook provides a runnable Julia walkthrough covering all the same cases with identical time instances.

## Computing parameters

Both regimes start from the same dimensional inputs (CGS units) and reduce to a small set of
nondimensional numbers. Julia exposes this through [`compute_gravity_parameters`](@ref) and
[`compute_cg_parameters`](@ref), which return mutable parameter structs consumed directly by
`solve`. The MATLAB equivalent computes the same quantities as plain variables.

```julia
using ForcedInterfacialWaves

# Pure gravity (α = 0)
pg = compute_gravity_parameters()
println("β        = ", pg.β)
println("F₀       = ", pg.F₀)
println("l_c [cm] = ", pg.l_c)
println("t_c [s]  = ", pg.t_c)

# Capillary-gravity (α > 0)
p = compute_cg_parameters()
println("α        = ", p.α)
println("kₛ       = ", p.kₛ)
println("kₗ       = ", p.kₗ)
println("F₀       = ", p.F₀)
```

```text
β        = 0.9980019980019981
F₀       = 0.0013888559222787725
l_c [cm] = 0.7269476668297654
t_c [s]  = 0.027221814475025485
α        = 0.13888559222787725
kₛ       = 1.1967001366781338
kₗ       = 6.010670937188218
F₀       = 0.0013888559222787723
```

```matlab
% Dimensional physical parameters (CGS units)
U     = 26.7046;   % Base flow speed [cm/s]
g     = 981.0;     % Gravitational acceleration [cm/s^2]
T     = 72.0;      % Surface tension [dyn/cm]
rho_l = 1.0;       % Lower-fluid density [g/cm^3]
rho_u = 0.001;     % Upper-fluid density [g/cm^3]

% Characteristic scales
l_c = U^2 / g;
t_c = U / g;
F_c = rho_l * U^2 * l_c;

% Pure-gravity parameters (alpha = 0)
rho_r = rho_u / rho_l;
beta  = (1.0 - rho_r) / (1.0 + rho_r);
F0_pg = 0.01 * T / F_c;

fprintf('beta      = %.6f\n', beta);
fprintf('F0        = %.6e\n', F0_pg);
fprintf('l_c [cm]  = %.6f\n', l_c);
fprintf('t_c [s]   = %.6f\n', t_c);

% Capillary-gravity parameters (alpha > 0)
alpha = T / (rho_l * U^2 * l_c);
disc  = (1.0 + rho_r)^2 - 4.0 * alpha * (1.0 - rho_r);
k_l   = ((1.0 + rho_r) + sqrt(disc)) / (2.0 * alpha);
k_s   = ((1.0 + rho_r) - sqrt(disc)) / (2.0 * alpha);
F0_cg = 0.01 * T / F_c;

fprintf('alpha     = %.6f\n', alpha);
fprintf('k_s       = %.6f\n', k_s);
fprintf('k_l       = %.6f\n', k_l);
fprintf('F0        = %.6e\n', F0_cg);
```

## Parameter names and mutation

Both parameter structs are **mutable** and expose the mathematical fields listed in the
[API reference](../api.md) under ASCII and Unicode names (for example, `p.alpha ≡ p.α`
and `p.k_l ≡ p.kₗ`), including for destructuring:

```julia
(; α, kₗ, kₛ) = p
α, kₗ, kₛ
```

```text
(0.13888559222787725, 6.010670937188218, 1.1967001366781338)
```

Mutating a field does **not** recompute quantities derived from it. For example, `p.α`
influences `k_l`, `k_s`, and `F0`; changing `p.α` in place leaves those derived fields
stale. Build a fresh parameter set with `compute_cg_parameters` or
`compute_gravity_parameters` keywords when changing a physical input:

```julia
p_new = compute_cg_parameters(T=80.0)   # consistently re-derives k_l, k_s, F0, ...
p_new.k_l, p_new.k_s
```

```text
(5.254642763643847, 1.2319912028358704)
```

Direct mutation is intended for quadrature controls without downstream dependents, such
as `atol`, `rtol`, `k_max`, and `epsilon_pv`.

Both parameter structs also carry the quadrature tolerances used throughout the theory pages:

| Parameter | Julia | MATLAB | Value |
|:----------|:------|:-------|------:|
| CPV exclusion radius $\varepsilon$ | `p.ε` / `pg.ε` | `epsilon_pv` | $10^{-6}$ |
| Absolute tolerance | `p.atol` / `pg.atol_cpv` | `AbsTol` | $10^{-10}$ |
| Relative tolerance | `p.rtol` / `pg.rtol_cpv` | `RelTol` | $10^{-8}$ |
| $K_{\max}$ (capillary–gravity) | `p.k_max` | `k_max` | $\infty$ |
| $K_{\max}$ (pure-gravity analytical) | `pg.k_max_analytical` | `KMAX_analytical` | $100$ |
| $K_{\max}$ (pure-gravity CPV) | `pg.k_max_cpv` | `KMAX_cpv` | $100$ |
| Front-band half-width | `pg.front_band` | `front_band` | $1.0$ |
