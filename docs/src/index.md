# ForcedInterfacialWaves.jl

*Julia implementation of the initial value problem (IVP) for pressure-forced interfacial waves in a two-fluid system.*

## Citation

> Kadari, V.K., Yewale, N., Farsoiya, P.K., Mayya, Y.S. & Dasgupta, R. (2026).
> Interfacial waves from pressure forcing: revisiting classical theories from an IVP perspective.
> *arXiv preprint* [arXiv:2605.12254](https://arxiv.org/abs/2605.12254).

```bibtex
@article{kadari2026interfacial,
  title={Interfacial waves from pressure forcing: revisiting classical
         theories from an {IVP} perspective},
  author={Kadari, Vinod Kumar and Yewale, Nikhil and Farsoiya, Palas Kumar
          and Mayya, Y. S. and Dasgupta, Ratul},
  journal={arXiv preprint arXiv:2605.12254},
  year={2026}
}
```

## Installation

`ForcedInterfacialWaves.jl` is not yet registered in the Julia General registry. Install it
directly from GitHub (this mutates your active environment, so it is shown as a plain,
non-executed snippet rather than a live example):

```julia
using Pkg
Pkg.add(url="https://github.com/yewalenikhil65/ForcedInterfacialWaves.jl.git")
using ForcedInterfacialWaves
```

## Problem

```@raw html
<figure style="text-align:center;">
  <img src="assets/Fig3.png" alt="Figure 10 comparison" style="max-width:100%;">
</figure>
```

As shown above (Fig. $3$ in the manuscript), a localised pressure $\tilde{p}_e = \tilde{F}_0\delta(\tilde{x})$ force is applied at the interface between two inviscid, incompressible fluids of infinite depth, both streams moving at uniform speed $U$ rightwards. In the absence of forcing, the interface is flat and remains at $\tilde{z}=0$. This package computes the resulting nondimensional interfacial displacement $\eta(x,t)$ for two physical regimes:

| Case | Surface tension | Poles | Manuscript figure |
|:-----|:--------------:|:------|:------:|
| Pure gravity ($\alpha = 0$) | absent | CPV at $k = \beta$ | 6 |
| Capillary–gravity ($\alpha > 0$) | active | Removable at $k_s$, $k_l$ (combined integrand cancels) | 10 |

**Pure gravity (α = 0):** The solution is decomposed analytically into $T_0$ (steady) and $T_1$–$T_4$ (transient) terms — implemented as [`T₀`](@ref)–[`T₄`](@ref), which branch automatically on the side of the wavefront. $T_3$ uses the Fresnel cosine and sine integrals (evaluated in closed form via `FresnelIntegrals.jl`, no quadrature). An independent direct numerical Cauchy principal value (CPV) evaluation verifies the analytical decomposition.

**Capillary–gravity (α > 0):** The combined integrand $\mathbb{I}(k; x, t)$ — which sums the steady part $\eta_s$ and the two transient integrands $\mathbb{I}_3$, $\mathbb{I}_4$ so that their singularities at the gravity root $k_s$ and capillary root $k_l$ cancel — is integrated over $[0,\infty)$ split around both poles. The classical steady solution via residues and a $G(x)$ integral is also computed for large-time comparison.

See the [Theory](theory/overview.md) section for the full mathematical derivation.

## Quick start

The primary interface is [`solve`](@ref) applied to a [`ForcedGravityProblem`](@ref) (pure
gravity) or [`ForcedGCProblem`](@ref) (capillary–gravity). Both accept a scalar `x` or a
vector spatial grid, and return a [`WaveSolution`](@ref) with fields `η`, `η_steady`,
`η_s_local`, `η_s_farfield`, `η_transient`, `x`, `t`.

```@example quickstart
using ForcedInterfacialWaves

# ─── Parameters ───
pg = compute_gravity_parameters()      # PureGravityParams (α = 0)
p  = compute_cg_parameters()           # CapillaryGravityParams (α > 0)
```

### Steady state

```@example quickstart
x = make_cg_xgrid(p; Nx=2001, xlim=(-15.0, 15.0))

# Without Rayleigh dissipation (symmetric, eqn. 3.11)
sol_s = solve(ForcedGCProblem(p, x); method=steady(rayleigh_dissipation=false))
sol_s.η_s_local, sol_s.η_s_farfield   # decomposed components

# With Rayleigh dissipation (asymmetric, eqn. 3.12)
sol_r = solve(ForcedGCProblem(p, x); method=steady(rayleigh_dissipation=true))
extrema(sol_r.η)
```

### Pure-gravity IVP

```@example quickstart
t_pg = 183.68                          # nondimensional time (Fig. 6)
x_pg = make_gravity_xgrid(pg; Nx=2001)

prof_pg = solve(ForcedGravityProblem(pg, x_pg, t_pg))
prof_pg.η, prof_pg.η_steady, prof_pg.η_transient  # full, steady, transient

# Individual T₀–T₄ terms at a point
T₀(-2.0, pg), T₁(-2.0, t_pg, pg), T₂(-2.0, t_pg, pg), T₃(-2.0, t_pg, pg), T₄(-2.0, t_pg, pg)
```

### Capillary–gravity IVP

```@example quickstart
t_cg = 367.35                          # nondimensional time (Fig. 8)

sol_cg = solve(ForcedGCProblem(p, x, t_cg); method=IVP())
sol_cg.η, sol_cg.η_steady, sol_cg.η_transient

# Asymmetric long-time reference decomposition
sol_asym = solve(ForcedGCProblem(p, x, t_cg); method=IVP(asym_cancel=true))
sol_asym.η_steady, sol_asym.η_transient
```

### Parameters: Unicode/ASCII aliases and mutability caveat

Both parameter structs are **mutable** and expose the mathematical fields listed in the
[API reference](api.md) under ASCII and Unicode names (e.g. `p.alpha ≡ p.α`,
`p.k_l ≡ p.kₗ`), including for destructuring:

```@example quickstart
(; α, kₗ, kₛ) = p
α, kₗ, kₛ
```

Mutating a field does **not** recompute quantities derived from it. For example, `p.α`
influences `k_l`, `k_s`, and `F0` — changing `p.α` in place leaves those derived fields
stale. Always build a fresh parameter set with `compute_cg_parameters`/
`compute_gravity_parameters` keywords when changing a physical input:

```@example quickstart
p_new = compute_cg_parameters(T=80.0)   # consistently re-derives k_l, k_s, F0, ...
p_new.k_l, p_new.k_s
```

## Interactive theory and usage notebook

For a guided, runnable companion covering the physical setup, parameter construction,
pure-gravity and capillary–gravity workflows, and figure generation, see the
[`forced_interfacial_waves_usage.ipynb` notebook](https://github.com/yewalenikhil65/ForcedInterfacialWaves.jl/blob/main/notebooks/forced_interfacial_waves_usage.ipynb).

## Generating the paper figures

Plotting requires Plots.jl, kept out of the package's runtime dependencies. Use the
dedicated `scripts` environment:

```julia-repl
julia --project=scripts -e 'using Pkg; Pkg.instantiate()'
```

```bash
julia -t auto --project=scripts scripts/figure6_pure_gravity.jl
julia -t auto --project=scripts scripts/figure10_capillary_gravity.jl
```

## Contents

```@contents
Pages = [
    "theory/overview.md",
    "theory/steady_state.md",
    "theory/pure_gravity.md",
    "theory/capillary_gravity.md",
    "api.md",
]
Depth = 2
```
