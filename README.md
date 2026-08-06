# ForcedInterfacialWaves.jl

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://yewalenikhil65.github.io/ForcedInterfacialWaves.jl/)

Julia implementation and open validation of the pressure-forced interfacial-wave IVPs documented in `jfm_matlab_codes.tex`:

- **Pure gravity, α = 0 (Figure 6)**
- **Capillary–gravity, α > 0 (Figure 10)**

The package keeps scalar functions corresponding directly to the MATLAB formulas, while vectorised profile functions use Julia-specific optimizations without changing the default MATLAB numerical parameters.

## Installation

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Dependencies are QuadGK.jl, FresnelIntegrals.jl, Plots.jl, and Printf.

## MATLAB-equivalent defaults

| Case | Parameter | Value |
|---|---:|---:|
| Pure gravity | `epsilon_cpv` | `1e-6` |
| Pure gravity | CPV `atol`, `rtol` | `1e-10`, `1e-8` |
| Pure gravity | transformed `atol`, `rtol` | `1e-10`, `1e-8` |
| Pure gravity | `k_max_analytical`, `k_max_cpv` | `100`, `100` |
| Pure gravity | `front_band` | `1.0` |
| Capillary–gravity | `epsilon_pv` | `1e-6` |
| Capillary–gravity | `atol`, `rtol` | `1e-10`, `1e-8` |
| Capillary–gravity | `k_max`, `k_max_steady` | `Inf`, `Inf` |

## Pure-gravity API

```julia
pg = compute_gravity_parameters()
t = 1.0 / pg.t_c

η_cpv = gravity_numerical_cpv(-2.0, t, pg)
η, η_s, η_tr = gravity_analytical_left(-2.0, t, pg)

T0 = gravity_T0(-2.0, pg)
T1 = gravity_T1_left(-2.0, pg)
T2 = gravity_T2_left(-2.0, t, t + 2.0, pg)
T3 = gravity_T3_left(-2.0, t, t + 2.0, pg)
T4 = gravity_T4_left(-2.0, t, t + 2.0, pg)
```

`T3` uses FresnelIntegrals.jl's `fresnelc` and `fresnels`, which have the same normalization as MATLAB and avoid numerical quadrature of the Fresnel kernels.

## Capillary–gravity API

```julia
using ForcedInterfacialWaves

p = compute_cg_parameters()
t = 3.0 / p.t_c
x = make_cg_xgrid(p; Nx=2001)

# Vectorised profile computation
η_ivp = compute_cg_ivp_profile(x, t, p)
η_steady = compute_cg_steady_profile(x, p)

# Compatibility wrapper; skip steady work when it is not plotted
η_ivp, η_steady = compute_cg_profile(
    x, t, p; compute_steady=true)

# Scalar MATLAB-reference path
eta_at_x = ivp_surface_elevation(3.0, t, p)
I1, I2, I3 = cg_partial_integrals(3.0, t, p)
```

### Vectorised profile algorithm

The scalar MATLAB algorithm runs three adaptive quadratures independently at every spatial point. The optimised Julia profile instead:

1. factors the combined integrand as
   `I(k;x,t) = C(k,t)cos(kx) + S(k,t)sin(kx)`;
2. computes χ(k), transient amplitudes, and time phases once per quadrature node;
3. uses `QuadGK.quadgk!` for allocation-efficient vector-valued integration;
4. uses the maximum component norm so tolerances apply pointwise to the worst profile component;
5. divides the spatial grid into thread-local chunks;
6. generates `sin(kx)` and `cos(kx)` on each uniform chunk with a recurrence, reset every 64 points to control roundoff.

The default `method=:threaded_vector` is the recommended path. Two reference methods remain available:

```julia
compute_cg_ivp_profile(x, t, p; method=:vector)
compute_cg_ivp_profile(x, t, p; method=:threaded_scalar)
```

The steady profile is computed only when requested. Figure 10 requests it only for `time_index == 300`.

### Measured performance

On 8 Julia threads with MATLAB-exact `k_max=Inf`, `atol=1e-10`, and `rtol=1e-8`:

- 201-point first-time benchmark: **1.995 s** threaded-vector vs **19.410 s** threaded-scalar (**9.73×**)
- maximum pointwise difference: **2.581e-10**
- actual 2000-point IVP profile: **9.16 s** at `time_index=1`
- actual 2000-point IVP profile: **8.51 s** at `time_index=300`
- actual 2000-point steady profile: **10.25 s**

Run profile scripts with all available threads:

```bash
julia -t auto --project=. scripts/figure10_capillary_gravity.jl
```

## Reproducing Figures 6 and 10

```bash
julia -t auto --project=. scripts/figure6_pure_gravity.jl
julia -t auto --project=. scripts/figure10_capillary_gravity.jl
```

Plots use Plots.jl with the GR backend.

## Integral-level open validation

```bash
julia --project=. scripts/validate_integrals.jl
```

This prints:

- T₀ and T₁–T₄ separately on both sides of the pure-gravity front;
- transformed analytical sum vs direct numerical CPV;
- capillary–gravity integrands `I_s`, `I3`, `I4`, and their combined value;
- capillary–gravity partial integrals below, between, and above the two roots;
- Fresnel reference values.

See [`docs/src/validation.md`](docs/src/validation.md) for the equation-to-function correspondence and numerical tables.

## Tests

```bash
julia -t auto --project=. test/runtests.jl
```

The suite validates parameters, Fresnel normalization, individual terms, analytical-vs-CPV agreement, scalar-vs-vector profile equivalence, and large-time IVP behavior.

## Repository layout

```text
ForcedInterfacialWaves.jl/
├── Project.toml
├── src/
│   ├── ForcedInterfacialWaves.jl
│   ├── parameters.jl
│   ├── pure_gravity.jl
│   └── capillary_gravity.jl
├── scripts/
│   ├── figure6_pure_gravity.jl
│   ├── figure10_capillary_gravity.jl
│   └── validate_integrals.jl
├── docs/
│   ├── make.jl
│   └── src/
│       ├── index.md
│       ├── theory.md
│       ├── validation.md
│       └── api.md
└── test/
    └── runtests.jl
```

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
