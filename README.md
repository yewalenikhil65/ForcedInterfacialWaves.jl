# ForcedInterfacialWaves.jl

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://yewalenikhil65.github.io/ForcedInterfacialWaves.jl/)

Julia solver for the initial value problem (IVP) of pressure-forced interfacial waves
between two inviscid, incompressible fluids in relative uniform motion, following
[Kadari et al. (2026)](https://arxiv.org/abs/2605.12254). Two regimes are implemented:

- **Pure gravity** (`α = 0`) — analytical `T₀`–`T₄` decomposition, cross-checked against a
  direct numerical Cauchy principal value (CPV) evaluation (manuscript Figure 6).
- **Capillary–gravity** (`α > 0`) — a combined-integrand technique that cancels the two
  removable poles at `k_s`, `k_l` before quadrature (manuscript Figure 10).

This package is **not registered** in the Julia General registry; install it directly from
GitHub (see below).

## Installation

Requires Julia ≥ 1.9.

```julia
using Pkg
Pkg.add(url="https://github.com/yewalenikhil65/ForcedInterfacialWaves.jl.git")
```

Or clone and activate locally:

```bash
git clone https://github.com/yewalenikhil65/ForcedInterfacialWaves.jl.git
cd ForcedInterfacialWaves.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Quick start

```julia
using ForcedInterfacialWaves

# ─── Pure gravity (α = 0) ───
pg = compute_gravity_parameters()
t  = 1.0 / pg.t_c          # nondimensional time for t_dim = 1 s

sol = solve(ForcedGravityProblem(pg, -2.0, t))   # full IVP
sol.η            # total displacement
sol.η_steady     # steady part
sol.η_transient  # transient part

steady_sol = solve(ForcedGravityProblem(pg, -2.0); method=steady(rayleigh_dissipation=false))
steady_sol.η    # steady solution only

x_grid = make_gravity_xgrid(pg; Nx=2001)
prof = solve(ForcedGravityProblem(pg, x_grid, t))  # full spatial profile
prof.η

# ─── Capillary–gravity (α > 0) ───
p = compute_cg_parameters()
t_cg = 3.0 / p.t_c

sol_cg = solve(ForcedGCProblem(p, 3.0, t_cg))
sol_cg.η             # full IVP displacement
sol_cg.η_steady      # symmetric ηₛ from equation (4.5b)
sol_cg.η_transient   # η - ηₛ, corresponding to equations (4.5c–d)

# Explicitly request the asymmetric-cancellation/reference decomposition
sol_cg_asym = solve(
    ForcedGCProblem(p, 3.0, t_cg);
    method=IVP(asym_cancel=true),
)
sol_cg_asym.η_steady   # asymmetric classical radiation profile
sol_cg_asym.η_transient # η - sol_cg_asym.η_steady

# The same option works for a spatial profile
x_grid_cg = make_cg_xgrid(p; Nx=2001)
prof_cg = solve(ForcedGCProblem(p, x_grid_cg, t_cg))
prof_cg_asym = solve(
    ForcedGCProblem(p, x_grid_cg, t_cg);
    method=IVP(asym_cancel=true),
)
```

### Inspecting the pure-gravity `T₀`–`T₄` terms

The analytical decomposition behind `ForcedGravityProblem` is exposed directly. All four
transient/steady terms branch automatically on the side of the wavefront (`x` vs. `t`) — no
separate left/right functions are needed:

```julia
pg = compute_gravity_parameters()
t, x = 1.0 / pg.t_c, -2.0

T₀(x, pg)        # steady term
T₁(x, t, pg)     # closed-form transient term
T₂(x, t, pg)     # exponentially-damped oscillatory integral
T₃(x, t, pg)     # Fresnel-integral term (via FresnelIntegrals.jl, no quadrature)
T₄(x, t, pg)     # oscillatory integral
```

## Parameters: Unicode and ASCII, mutable, with a caveat

`CapillaryGravityParams` and `PureGravityParams` are **mutable** structs. The mathematical
fields listed below have both ASCII and Unicode names, and both forms work for reading,
writing, and destructuring:

```julia
p = compute_cg_parameters()

p.alpha == p.α        # true — ASCII/Unicode alias, both readable
p.α = 0.2              # writable through either name
p.alpha                # == 0.2

(; α, β, kₗ, kₛ, F₀) = p   # destructuring works with Unicode names too
```

| Struct | Unicode | ASCII |
|:-------|:--------|:------|
| `CapillaryGravityParams` | `α`, `β`, `ρᵣ`, `γ_ρ`, `kₗ`, `kₛ`, `F₀`, `ε` | `alpha`, `beta`, `rho_r`, `gamma_rho`, `k_l`, `k_s`, `F0`, `epsilon_pv` |
| `PureGravityParams` | `β`, `sqrtβ`, `ρᵣ`, `F₀`, `ε` | `beta`, `sqrt_beta`, `rho_r`, `F0`, `epsilon_cpv` |

Typing hints (Julia REPL / VS Code): `\alpha<TAB>` → α, `\beta<TAB>` → β, `\rho<TAB>` → ρ,
`\gamma<TAB>` → γ, `\varepsilon<TAB>` → ε, `\_l<TAB>` → ₗ, `\_s<TAB>` → ₛ, `\_r<TAB>` → ᵣ,
`\_0<TAB>` → ₀.

> [!WARNING]
> **Mutating a parameter does not recompute dependent quantities.** These structs store
> raw fields; setting `p.alpha = 0.2` after construction leaves `k_l`, `k_s`, `F0`, and every
> other derived field at their old values, since nothing re-runs the derivation. For any
> physically consistent change, build a *fresh* parameter set instead:
>
> ```julia
> p = compute_cg_parameters(T=80.0, F0_factor=0.02)   # new derived k_l, k_s, F0, ...
> pg = compute_gravity_parameters(rho_u=0.0005)
> ```
>
> Direct mutation is intended only for tolerances/quadrature knobs (`atol`, `rtol`,
> `k_max`, `epsilon_pv`, ...) that have no downstream dependents.

## Reproducing the manuscript figures

Figure generation depends on Plots.jl, which is **not** a runtime dependency of the
package (root `Project.toml` only requires CommonSolve, FresnelIntegrals, and QuadGK).
`scripts/Project.toml` provides a dedicated environment with `ForcedInterfacialWaves`
(via a local `[sources]` path pointing at the repository root) and `Plots`:

```bash
julia --project=scripts -e 'using Pkg; Pkg.instantiate()'   # one-time setup
julia -t auto --project=scripts scripts/figure6_pure_gravity.jl
julia -t auto --project=scripts scripts/figure10_capillary_gravity.jl
```

`-t auto` enables multithreaded quadrature for the spatial profiles. Output PNG/PDF files
are written to `output/figure6/` and `output/figure10/`.

## Documentation and tests

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl        # build the documentation locally

julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia -t auto --project=. test/runtests.jl   # run the test suite
```

Full docs (theory and API reference): <https://yewalenikhil65.github.io/ForcedInterfacialWaves.jl/>

For a guided, runnable theory and usage walkthrough (including Julia/MATLAB-oriented
figure reproduction), open [`notebooks/forced_interfacial_waves_usage.ipynb`](notebooks/forced_interfacial_waves_usage.ipynb).

## Dependencies

Runtime (root `Project.toml`):

| Package | Purpose |
|:--------|:--------|
| [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl) | Adaptive Gauss–Kronrod quadrature |
| [FresnelIntegrals.jl](https://github.com/kiranshila/FresnelIntegrals.jl) | Closed-form Fresnel cosine/sine integrals |
| [CommonSolve.jl](https://github.com/SciML/CommonSolve.jl) | `solve(prob)` interface convention |

Docs/scripts only (`docs/Project.toml`): [Plots.jl](https://github.com/JuliaPlots/Plots.jl)
(figure generation), [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl),
[LaTeXStrings.jl](https://github.com/stevengj/LaTeXStrings.jl), [LiveServer.jl](https://github.com/tlienart/LiveServer.jl).

## Repository layout

```text
ForcedInterfacialWaves.jl/
├── Project.toml
├── src/
│   ├── ForcedInterfacialWaves.jl
│   ├── parameters.jl
│   ├── solve.jl                  # ForcedGravityProblem, ForcedGCProblem, WaveSolution, solve
│   ├── pure_gravity.jl           # T₀–T₄, Fresnel terms, numerical CPV
│   └── capillary_gravity.jl      # combined integrand, partial integrals, profiles
├── scripts/
│   ├── Project.toml               # ForcedInterfacialWaves (local) + Plots
│   ├── figure6_pure_gravity.jl
│   └── figure10_capillary_gravity.jl
├── matlab/
│   ├── fig6_jfm_vinod.m           # authoritative MATLAB Figure 6 reference
│   └── Fig10_jfm_vinod.m          # authoritative MATLAB Figure 10 reference
├── docs/
│   ├── Project.toml
│   ├── make.jl
│   └── src/
│       ├── index.md
│       ├── theory/
│       │   ├── overview.md
│       │   ├── steady_state.md          # → steady_proof.md
│       │   ├── pure_gravity.md          # → pure_gravity_steady_proof.md, pure_gravity_transient_proof.md
│       │   └── capillary_gravity.md     # → capillary_gravity_asymmetric_cancellation.md, lamb_gx_equivalence.md
│       ├── steady_proof.md                          # hidden: contour-integration proof
│       ├── pure_gravity_steady_proof.md             # hidden: contour-integration proof
│       ├── pure_gravity_transient_proof.md          # hidden: contour-integration proof
│       ├── capillary_gravity_asymmetric_cancellation.md  # hidden: stationary-phase proof
│       ├── lamb_gx_equivalence.md                   # hidden: residue-theorem proof
│       └── api.md
├── notebooks/
│   └── forced_interfacial_waves_usage.ipynb          # guided theory/usage walkthrough
└── test/
    ├── runtests.jl
    ├── benchmark_test.jl
    └── validation_table.jl
```

The five "hidden" proof pages under `docs/src/` are linked from the corresponding
`theory/*.md` pages but omitted from the navigation sidebar (they are not listed in
`docs/make.jl`'s `pages`); they are still built and reachable by URL/link for readers who
want the full contour-integration derivations.

## Citation

If you use this package, please cite:

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

## Contributors

- Vinod Kumar Kadari
- Nikhil Yewale
- Prof. Y.S. Mayya
- Prof. Ratul Dasgupta
