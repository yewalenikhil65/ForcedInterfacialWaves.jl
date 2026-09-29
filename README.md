# ForcedInterfacialWaves

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://yewalenikhil65.github.io/ForcedInterfacialWaves.jl/)

Theory, derivations, and numerical validation for the initial value problem (IVP) of
pressure-forced interfacial waves between two inviscid, incompressible fluids in relative
uniform motion, following [Kadari et al. (2026)](https://arxiv.org/abs/2605.12254). Two
regimes are covered:

- **Pure gravity** (`α = 0`) — analytical `T₀`–`T₄` decomposition, cross-checked against a
  direct numerical Cauchy principal value (CPV) evaluation (manuscript Figure 5).
- **Capillary–gravity** (`α > 0`) — a combined-integrand technique that cancels the two
  removable poles at `k_s`, `k_l` before quadrature (manuscript Figures 6–9), validated
  against nonlinear Basilisk (Navier–Stokes/VOF) simulations.

The site is built with [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) from
self-contained Julia/MATLAB code blocks embedded directly in the theory pages — there is no
installable Julia package here. For a guided, runnable walkthrough of the equivalent
computations, see the notebooks described below.

## Building the documentation

Requires Julia ≥ 1.9.

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate()'   # one-time setup
julia --project=docs docs/make.jl                        # build into docs/build/
```

To preview with live-reload while editing:

```bash
julia --project=docs docs/serve.jl
```

Full docs: <https://yewalenikhil65.github.io/ForcedInterfacialWaves.jl/>

## Notebooks

[`notebooks/`](notebooks/) contains Jupyter notebooks and supporting data reproducing the
gravity–capillary wave comparisons presented in the paper:

- [`forced_interfacial_waves_usage.ipynb`](notebooks/forced_interfacial_waves_usage.ipynb) —
  guided usage template following the paper's mathematical organisation (steady-state
  Fourier representation, pure-gravity `T₀`–`T₄` IVP, capillary–gravity IVP via numerical
  quadrature).
- [`gc_comoving_ivp_readable.ipynb`](notebooks/gc_comoving_ivp_readable.ipynb) — solves the
  gravity–capillary co-moving IVP and compares it against a Basilisk VOF simulation
  (`basilisk_gc_ivp/`).

See [`notebooks/README.md`](notebooks/README.md) for details, data provenance, and
instructions to regenerate the Basilisk interface data.

## Repository layout

```text
ForcedInterfacialWaves.jl/
├── README.md
├── docs/
│   ├── Project.toml            # Documenter, LiveServer
│   ├── make.jl                 # site build
│   ├── serve.jl                # local live-reload preview
│   └── src/
│       ├── index.md
│       ├── theory/
│       │   ├── overview.md
│       │   ├── steady_state.md               # → steady_proof.md, steady_proof_withRayleigh.md
│       │   ├── pure_gravity.md                # → pure_gravity_steady_proof.md, pure_gravity_transient_proof.md,
│       │   │                                  #   pure_gravity_transient_asymptotic_proof.md
│       │   ├── capillary_gravity.md           # → capillary_gravity_asymmetric_cancellation.md,
│       │   │                                  #   capillary_gravity_rayleigh_dissipation.md, lamb_gx_equivalence.md
│       │   └── basilisk_capillary_gravity.md  # Basilisk CFD setup for Fig. 9
│       ├── validation/
│       │   └── cfd_conformal_mapping.md       # nonlinear-regime CFD vs. conformal-mapping IVP
│       ├── steady_proof.md                          # hidden: contour-integration proof
│       ├── steady_proof_withRayleigh.md             # hidden: contour-integration proof (Rayleigh dissipation)
│       ├── pure_gravity_steady_proof.md             # hidden: contour-integration proof
│       ├── pure_gravity_transient_proof.md          # hidden: contour-integration proof
│       ├── pure_gravity_transient_asymptotic_proof.md  # hidden: asymptotic analysis
│       ├── capillary_gravity_asymmetric_cancellation.md  # hidden: stationary-phase proof
│       ├── capillary_gravity_rayleigh_dissipation.md     # hidden: Rayleigh dissipation treatment
│       ├── lamb_gx_equivalence.md                   # hidden: residue-theorem proof
│       └── assets/                                  # figures, CSV reference data, animations
└── notebooks/
    ├── README.md
    ├── forced_interfacial_waves_usage.ipynb
    ├── gc_comoving_ivp_readable.ipynb
    ├── if_*.csv                             # interfacial wave data for the usage notebook
    └── basilisk_gc_ivp/                     # Basilisk simulation bundle (source, ICs, output)
```

The "hidden" proof pages under `docs/src/` are linked from the corresponding `theory/*.md`
pages but omitted from the navigation sidebar (they are not listed in `docs/make.jl`'s
`pages`); they are still built and reachable by URL/link for readers who want the full
derivations.

## Citation

If you use this material, please cite:

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
