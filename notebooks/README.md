# Notebooks

This directory contains Jupyter notebooks and supporting data for
reproducing the gravity–capillary wave comparisons presented in the paper.

---

## Contents

```
notebooks/
├── README.md                              ← this file
├── forced_interfacial_waves_usage.ipynb   ← ForcedInterfacialWaves.jl usage guide
├── if_*.csv                               ← interfacial wave data for forced_interfacial_waves_usage.ipynb
├── gc_comoving_ivp_readable.ipynb         ← conformal IVP vs. Basilisk VOF comparison
└── basilisk_gc_ivp/                       ← Basilisk simulation bundle
    ├── run_gc_ivp.c                       ← Basilisk IVP solver (Step 4)
    ├── params.h                           ← physical/numerical parameters
    ├── FreeSurface.dat                    ← initial VOF interface geometry
    ├── velocity_interpolated_below.dat    ← initial liquid-cell velocity field
    └── interface_data/
        └── interface-0.dat ... interface-26.dat   ← Basilisk output_facets() snapshots
```

---

## 1. `forced_interfacial_waves_usage.ipynb`

A guided Julia usage template for the `ForcedInterfacialWaves.jl` package.
It follows the paper's mathematical organisation and demonstrates:

- The steady-state Fourier representation of pressure-forced interfacial waves.
- The pure-gravity (α = 0) IVP with the analytical T₀–T₄ decomposition.
- The capillary–gravity (α > 0) IVP with numerical quadrature.

The `if_*.csv` files in this directory provide precomputed interfacial
wave profiles used by `forced_interfacial_waves_usage.ipynb`.

---

## 2. `gc_comoving_ivp_readable.ipynb`

A self-contained, readable notebook that:

1. Solves the N = 1024 pure-gravity steady wave at steepness ε = 0.9 via
   Newton continuation.
2. Seeds the gravity–capillary co-moving IVP and integrates it in time
   using `OrdinaryDiffEq.jl` (ROK4a with Krylov linear solves).
3. Loads native Basilisk PLIC facet endpoints directly from
   `basilisk_gc_ivp/interface_data/interface-<index>.dat`.
4. Produces an animated GIF (`gc_ivp_comparison.gif`) comparing the
   conformal-mapping IVP surface (blue) with the Basilisk VOF interface
   (red scatter) at each matched time step.

### Physical parameters (CGS)

| Parameter     | Value                |
|---------------|----------------------|
| N             | 1024                 |
| ε (steepness) | 0.9                  |
| λ             | 3.194 cm             |
| c             | 24 cm/s              |
| g             | 981 cm/s²            |
| ρ_water       | 1.0 g/cm³            |
| σ             | 72 dyn/cm            |
| B (Bond)      | 0.00719              |
| Δt_save       | 0.01 s               |
| T_final       | 0.18 s (19 frames)   |

### Julia dependencies

```
FFTW, OrdinaryDiffEq, OrdinaryDiffEqRosenbrock, LinearSolve,
NonlinearSolve, Plots, LaTeXStrings, Printf, DelimitedFiles
```

---

## 3. `basilisk_gc_ivp/` — Basilisk simulation bundle

This directory contains everything needed to run (or re-run) the Basilisk
VOF simulation that produces the interface data compared in the notebook.

### Pre-supplied initial conditions

The two initial-condition files are already computed and included:

- `FreeSurface.dat` — dimensional surface profile for Basilisk `distance()`
  (generated upstream by `generate_ic.jl` in InterfacialWaves.jl).
- `velocity_interpolated_below.dat` — velocity (u, v) at each liquid cell
  centre, one CSV row per cell in Basilisk `foreach()` traversal order
  (generated upstream by `compute_velocity_below.jl`).

### How to regenerate `interface_data/`

If you need to re-run the Basilisk simulation from the provided initial
conditions:

```bash
cd basilisk_gc_ivp/

# Compile (requires Basilisk/qcc)
qcc -O2 -o run_gc_ivp run_gc_ivp.c -lm

# Or with MPI:
# CC99='mpicc -std=c99' qcc -D_MPI=1 -O2 -o run_gc_ivp run_gc_ivp.c -lm

# Run
./run_gc_ivp
# MPI: mpirun -np 4 ./run_gc_ivp
```

This reads `FreeSurface.dat` and `velocity_interpolated_below.dat`, then
writes snapshots every 0.01 s to:

```
interface_data/interface-0.dat
interface_data/interface-1.dat
...
```

Each file contains Basilisk `output_facets()` PLIC line segments:

```
x1 y1
x2 y2
                      ← blank line separates facets
x1 y1
x2 y2
...
```

### Simulation setup (from `run_gc_ivp.c` and `params.h`)

- Domain: 2λ × 2λ, periodic in x, free-slip top/bottom.
- Two-phase Navier–Stokes (inviscid: μ₁ = μ₂ = 0).
- Adaptive mesh: maxlevel = 11, minlevel = 7.
- ρ_air = ρ_water / 1000.

---

## Workflow summary

```
                   ┌─────────────────────────────────┐
                   │  Pre-computed initial conditions │
                   │  (FreeSurface.dat + velocity)    │
                   └────────────┬────────────────────┘
                                │
                                ▼
                   ┌─────────────────────────────────┐
                   │  Basilisk: run_gc_ivp.c          │
                   │  qcc + ./run_gc_ivp              │
                   └────────────┬────────────────────┘
                                │
                                ▼
                   ┌─────────────────────────────────┐
                   │  interface_data/interface-*.dat   │
                   │  (native PLIC facet endpoints)   │
                   └────────────┬────────────────────┘
                                │
                                ▼
              ┌─────────────────────────────────────────┐
              │  gc_comoving_ivp_readable.ipynb          │
              │                                         │
              │  1. Solve pure-gravity seed (Newton)    │
              │  2. Seed GC IVP, integrate (ROK4a)      │
              │  3. Load direct Basilisk facets          │
              │  4. Animate conformal vs. Basilisk       │
              │     → gc_ivp_comparison.gif              │
              └─────────────────────────────────────────┘
```
