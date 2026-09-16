# Theory

This section contains the mathematical formulation and numerical implementation of the forced interfacial wave IVP, following the notation of [Kadari et al. (2026)](https://arxiv.org/abs/2605.12254).

## Problem description and nondimensionalisation

```@raw html
<figure style="text-align:center;">
  <img src="../../assets/Fig3.png" alt="Pressure forcing at a two-fluid interface" style="max-width:100%;">
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

The parameter $\alpha$ measures the relative importance of surface tension. 

## Pages in this section

- **[Steady-State](steady_state.md)** — Manuscript §3.4: the general Fourier integral solution (eqns. 3.7–3.8), the steady far-field and local decomposition (eqn. 3.11), the Rayleigh dissipation form (eqn. 3.12), and numerical evaluation with Julia/MATLAB code.

- **[Pure Gravity](pure_gravity.md)** — Manuscript §4.1: setting $\alpha = 0$, the analytical $T_0$–$T_4^\pm$ decomposition, Fresnel-integral representations, independent CPV verification, and the spatial IVP profile at $t = 183.68$ (manuscript Fig. 6). Includes Julia/MATLAB implementations.

- **[Capillary–Gravity](capillary_gravity.md)** — Manuscript §4.2: the finite-capillarity case with the combined integrand technique, steady-state asymmetric cancellation proof, the $\mathbb{I}_4$ transient decay (Fig. 7), the full IVP profile (Fig. 8), and comparison with nonlinear simulations (Fig. 10). Includes Julia/MATLAB implementations.

Each sub-page is self-contained: it introduces the relevant equations, presents the derivation, and provides working code.

## Specifying the parameters

We begin by specifying the dimensional physical parameters $U,\ g,\ T,\ \rho_l,\ \rho_u$, from
which we form the characteristic scales $l_c,\ t_c,\ F_c$ and then the nondimensional groups:
the pure-gravity parameters ($\rho_r,\ \beta,\ F_0$) and the capillary–gravity parameters
($\alpha,\ k_l,\ k_s,\ F_0$, with $\Delta$ the discriminant of the root equation), as follows
in Julia/MATLAB (click below to see the respective code).

```julia
# Dimensional physical parameters (CGS units)
U  = 26.7046   # base flow speed            [cm/s]
g  = 981.0     # gravitational acceleration [cm/s²]
T  = 72.0      # surface tension            [dyn/cm]
ρₗ = 1.0       # lower-fluid density        [g/cm³]
ρᵤ = 0.001     # upper-fluid density        [g/cm³]

# Characteristic scales
l_c = U^2 / g
t_c = U / g
F_c = ρₗ * U^2 * l_c

# Pure-gravity parameters (α = 0)
ρᵣ    = ρᵤ / ρₗ
β     = (1 - ρᵣ) / (1 + ρᵣ)
F₀_pg = 0.01 * T / F_c

println("β        = ", β)
println("F₀       = ", F₀_pg)
println("l_c [cm] = ", l_c)
println("t_c [s]  = ", t_c)

# Capillary-gravity parameters (α > 0)
α     = T / (ρₗ * U^2 * l_c)
Δ     = (1 + ρᵣ)^2 - 4α * (1 - ρᵣ)
kₗ    = ((1 + ρᵣ) + √Δ) / (2α)
kₛ    = ((1 + ρᵣ) - √Δ) / (2α)
F₀_cg = 0.01 * T / F_c

println("α        = ", α)
println("kₛ       = ", kₛ)
println("kₗ       = ", kₗ)
println("F₀       = ", F₀_cg)
```

```text
β        = 0.9980019980019981
F₀       = 0.0013888559222787723
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

