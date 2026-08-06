# Implementation Notes — WavesDelta.jl

Complete context from the implementation session. Use as reference for the next session.

---

## Repository location

```
/home/nikhil/JFM_Vinod_supplementary/WavesDelta.jl/
```

Reference LaTeX files:
- `../jfm_matlab_codes_in_Latex_Vinod/jfm_matlab_codes.tex` — MATLAB listings + formulas
- `../supplementary_vinod_jfm.tex` — Paper supplementary (canonical notation)

---

## Paper citation

```bibtex
@article{kadari2026interfacial,
  title={Interfacial waves from pressure forcing: revisiting classical theories from an IVP perspective},
  author={Kadari, Vinod Kumar and Yewale, Nikhil and Farsoiya, Palas Kumar and Mayya, YS and Dasgupta, Ratul},
  journal={arXiv preprint arXiv:2605.12254},
  year={2026}
}
```

---

## Mathematical formulation (supplementary notation)

### Pure gravity (α = 0)

Pole at $k = \beta$. Analytical decomposition:
- $T_0$: steady (Laplace-type integral + closed sin term)
- $T_1^\pm$: closed-form $\mp\sin(\beta x)/(1+\rho_r)$
- $T_2^\pm$: exponentially-damped oscillatory integral (quadgk, Kmax=100)
- $T_3^\pm$: Fresnel integrals C(X), S(X) — **no quadrature** (O(1) via FresnelIntegrals.jl)
- $T_4^\pm$: oscillatory integral (quadgk, Kmax=100)

Direct numerical CPV: $G(k;x,t)$ integrand split around $k=\beta$

Left region uses $(0.5 - C)$, $(0.5 - S)$; right region uses $(0.5 + C)$, $(0.5 + S)$.

---

### Capillary–gravity (α > 0)

Nondimensionalisation: $l_c = U^2/g$, $t_c = U/g$, $F_c = \rho_l U^2 l_c$.

Parameters:
- $\alpha = g\tilde{T}/(\rho_l U^4)$ — surface tension parameter
- $\rho_r = \rho_u/\rho_l$ — density ratio
- $\beta = (1-\rho_r)/(1+\rho_r)$ — Atwood number
- $\gamma_\rho = 1/(1+\rho_r)$

Dispersion: $\chi(k) = \sqrt{\alpha k^3/(1+\rho_r) + \beta k}$

Poles: roots of $\alpha k^2 - (1+\rho_r)k + (1-\rho_r) = 0$ → $k_l$ (capillary), $k_s$ (gravity)

Formal solution (eq. 2.17): $\eta = \eta_s + \eta_{\mathrm{tr}}$

The MATLAB implementation combines all terms into one integrand $\mathbb{I}(k;x,t)$ before integration so poles cancel:

$$\mathbb{I} = \frac{2\cos(kx)}{\alpha(k-k_l)(k-k_s)} - \frac{(1+\rho_r)(k+\chi)}{1+\alpha k^2-\rho_r}\frac{\cos[k(t-x)-t\chi]}{\alpha(k-k_l)(k-k_s)} - \frac{(1+\rho_r)(k-\chi)}{1+\alpha k^2-\rho_r}\frac{\cos[k(t-x)+t\chi]}{\alpha(k-k_l)(k-k_s)}$$

$$\eta(x,t) = -\frac{F_0}{2\pi}\int_0^\infty \mathbb{I}\,dk$$

Split: $\int_0^{k_s-\varepsilon} + \int_{k_s+\varepsilon}^{k_l-\varepsilon} + \int_{k_l+\varepsilon}^\infty$

Steady solution: $G(x) = \frac{1}{k_l-k_s}\int_0^\infty[\cos(kx)/(k+k_s) - \cos(kx)/(k+k_l)]dk$

## Numerical parameters (identical to MATLAB)

| Parameter | CG value | PG value |
|:----------|:--------:|:--------:|
| ε (pole exclusion) | 1e-6 | 1e-6 |
| atol | 1e-10 | 1e-10 |
| rtol | 1e-8 | 1e-8 |
| k_max | Inf | — |
| K_max (analytical/CPV) | — | 100 |
| k_max_steady | Inf | — |
| front_band δ_x | — | 1.0 |

---

## Julia implementation details

### Dependencies

- `QuadGK.jl` — adaptive Gauss–Kronrod (quadgk, quadgk!)
- `FresnelIntegrals.jl` — fresnelc, fresnels (erf-based, O(1))
- `Plots.jl` + `LaTeXStrings.jl` — plotting (GR backend)
- `Documenter.jl` — documentation (KaTeX math engine)

### Source structure

```
src/WavesDelta.jl         — module, imports (quadgk, quadgk!, alloc_segbuf, fresnelc, fresnels)
src/parameters.jl         — CapillaryGravityParams, PureGravityParams structs + constructors
src/pure_gravity.jl       — T0–T4, Fresnel, CPV, profile
src/capillary_gravity.jl  — dispersion_chi, integrands, quadrature, profile algorithms
```

### Key design decisions

1. **Immutable parameter structs** — all physics + numerics in one struct, passed everywhere
2. **Individual integrands exposed** — `cg_integrand_steady`, `cg_integrand_I3`, `cg_integrand_I4` for validation
3. **Partial integrals exposed** — `cg_partial_integrals(x, t, p)` returns (I1, I2, I3)
4. **Immutable callable struct** for thread-safe scalar quadrature: `CGScalarIntegrand`
5. **Vector-valued quadgk!** for production profiles — factors integrand as $C(k,t)\cos(kx) + S(k,t)\sin(kx)$
6. **Trigonometric recurrence** on uniform grids — replaces sincos per point with 2 muls + 2 adds, reset every 64 points
7. **Maximum-component norm** — `cg_profile_norm(v) = maximum(abs, v)` so tolerances apply pointwise
8. **Threaded spatial chunks** — `Threads.@threads :static` over `cld(N, nthreads())` chunks

### Profile algorithms (`compute_cg_ivp_profile`)

| Method | Description | Use case |
|:-------|:------------|:---------|
| `:threaded_vector` (default) | One `quadgk!` per thread-chunk, factored integrand + recurrence | Production |
| `:vector` | Single `quadgk!` over full grid (serial) | Testing |
| `:threaded_scalar` | Independent `quadgk` per x-point, thread-parallel | Reference validation |

### Performance (8 threads, Nx=2000, k_max=Inf, atol=1e-10, rtol=1e-8)

- Threaded vector: **~9 s** per time index
- Threaded scalar: **~75 s** per time index  
- Speedup: **~8–10×**
- Pointwise difference vector vs scalar: **2.6e-10**

### @fastmath

Applied to `cg_combined_integrand` and `gravity_combined_integrand`. Safe because integrands are smooth away from excluded poles. Gives ~5-10% speedup on inner loop.

### QuadGK order=15

Used for:
- CG tail integral [k_l+ε, k_max] — highly oscillatory
- PG T4 integrals — oscillatory with quadratic phase
- PG CPV above-pole integral [β+ε, Kmax]

Default order=7 used for smooth/short intervals.

### QuadGK Inf handling

Uses coordinate transform $k = a + t/(1-t)$. Works because capillary dispersion $\chi \sim k^{3/2}$ gives super-linear phase growth → integrand effectively decays. But requires many subdivisions at tight tolerances. MATLAB's `integral` handles this better for oscillatory semi-infinite integrals.

**Practical workaround:** `compute_cg_parameters(k_max=50.0)` converges the integral to ~0.1% with ~100× fewer evaluations.

---

## Validation results

### Parameters: Julia = MATLAB to 15 digits

### Pure gravity T-terms (x=-2, t_dim=1s):
- T0: 15 digits (same quadrature target)
- T1: 15 digits (closed-form)
- T2: 11 digits (different adaptive paths)
- T3: 14 digits (both use fresnelc/fresnels — nearly identical)
- T4: 11 digits (different adaptive paths)
- η_analytical: 12 digits
- η_CPV: 10 digits
- |analytical − CPV| = 2.3250e-06 (identical in both)

### Integrand pointwise: 15 digits (1 ULP difference at k=2, x=3, t=110)

### Partial integrals (x=3, t=110):
- I1: 12+ digits agreement
- I2: 12 digits agreement  
- I3: ~5 digits (MATLAB hits interval limit on [k_l,∞); Julia resolves more)
---

## Documenter.jl setup

```
docs/
├── make.jl        — Documenter.HTML with KaTeX, pages list
├── serve.jl       — LiveServer at 0.0.0.0:8001
├── Project.toml   — deps: Documenter, LiveServer, Plots, LaTeXStrings, WavesDelta
└── src/
    ├── index.md          — Citation, problem, quick start
    ├── theory.md         — Full LaTeX math (supplementary notation)
    ├── julia_matlab.md   — Side-by-side code + outputs + plots
    ├── validation.md     — Numerical tables
    └── api.md            — @docs autodocs
```

Build: `julia --project=docs docs/make.jl`
Serve: `julia --project=docs docs/serve.jl` → http://<IP>:8001

Docstrings use `\\` escaping for LaTeX in triple-quoted strings:
```julia
"""
    my_func(x)

```math
\\eta = \\int_0^\\infty f(k)\\,dk
```
"""
```

---

## Scripts

| Script | Purpose |
|:-------|:--------|
| `scripts/figure6_pure_gravity.jl` | Figure 6 driver (Plots.jl) |
| `scripts/figure10_capillary_gravity.jl` | Figure 10 driver (Plots.jl) |
| `scripts/validate_integrals.jl` | Print all T-terms and partial integrals |
| `scripts/validate_output.jl` | Formatted output for docs comparison |
| `scripts/generate_comparison_plots.jl` | Overlay Julia + MATLAB CSV data |
| `matlab/validate_output.m` | MATLAB equivalent of validate_output.jl |
| `matlab/generate_comparison_data.m` | MATLAB CSV export for overlay plots |
| `matlab/fig6_jfm_vinod.m` | Extracted authoritative MATLAB Fig6 |
| `matlab/Fig10_jfm_vinod.m` | Extracted authoritative MATLAB Fig10 |

---

## Tests

```bash
julia -t auto --project=. test/runtests.jl
```

57 tests total:
- Fresnel values (8)
- PG parameters (4)
- PG T0 symmetry (2)
- PG analytical vs CPV left (3)
- PG analytical vs CPV right (2)
- PG individual terms finite (9)
- CG parameters (8)
- CG dispersion consistency (2)
- CG IVP→steady convergence (12)
- CG optimized profile equivalence (7): integrand, threaded_vector vs scalar, steady, compute_steady=false

---

## Known issues / TODO for next session

1. **MATLAB CSV comparison plots not yet generated** — need to run `matlab/generate_comparison_data.m` and copy CSVs to `docs/src/assets/`, then run `scripts/generate_comparison_plots.jl`
2. **Figure 6 at t_dim=1s** — entire visible domain is behind the front (wavefront at x≈36.7 is beyond plot range ±12.6). Consider using t_dim=0.07s for a plot showing both left and right regions.
3. **Remaining ASCII docstrings** — some internal/helper functions still have plain-text docs (profile integrands, chunk helpers). Not user-facing.
4. **`OscillatoryIntegralsODE.jl`** — investigated but incompatible (supports only constant-frequency harmonic/Bessel kernels; our phase is nonlinear in k).
5. **DelimitedFiles** added to Project.toml for CSV reading in comparison plot script.
6. **Plotting defaults** used in scripts:
```julia
plot_font = "Computer Modern"
default(fontfamily=plot_font, margin=6Plots.mm, linewidth=3,
        framestyle=:box, label=nothing, color="blue", grid=false,
        fg_legend=false, background_color_legend=false)
```
