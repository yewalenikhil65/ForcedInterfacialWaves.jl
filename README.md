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

## Citation

If you use this material, please cite:

> Kadari, V.K., Yewale, N., Farsoiya, P.K., Mayya, Y.S. & Dasgupta, R. (2026).
> Interfacial waves from pressure forcing: revisiting classical theories from an IVP perspective.
> *arXiv preprint* [arXiv:2605.12254](https://arxiv.org/abs/2605.12254).

```bibtex
@article{kadari2026interfacialwavespressureforcing,
      title={Interfacial waves from pressure forcing: role of the transients}, 
      author={Vinod Kumar Kadari and Nikhil Yewale and Palas Kumar Farsoiya and Y. S. Mayya and Ratul Dasgupta},
      year={2026},
      eprint={2605.12254},
      archivePrefix={arXiv},
      primaryClass={physics.flu-dyn},
      url={https://arxiv.org/abs/2605.12254}, 
}
```

## Contributors

- Vinod Kumar Kadari
- Nikhil Yewale
- Prof. Y.S. Mayya
- Prof. Ratul Dasgupta

  
## Building the documentation (contributors only)

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
