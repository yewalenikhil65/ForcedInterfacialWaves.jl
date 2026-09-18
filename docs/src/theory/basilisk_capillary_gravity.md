# Basilisk CFD: Capillary-Gravity Waves from Pressure Forcing

This page describes the Basilisk[^1] simulation used to generate the nonlinear CFD reference data compared against the IVP theory in the [Capillary–Gravity](capillary_gravity.md) section (Fig. 10 of the manuscript). The simulation solves the two-phase incompressible Navier–Stokes equations with surface tension for a localized Lorentzian pressure forcing applied at the interface between two fluids in uniform horizontal motion. The simulation is dimensional in CGS units. The full source file is [`notebooks/capillary_gravity_forced.c`](https://github.com/yewalenikhil65/ForcedInterfacialWaves.jl/blob/main/notebooks/capillary_gravity_forced.c).

---

## Headers

```c
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "reduced.h"
#include "tension.h"
```

The Basilisk headers provide: the cell-centred Navier–Stokes solver (`centered.h`), the two-phase VOF interface tracking (`two-phase.h`), the reduced-gravity body-force formulation (`reduced.h`), and the continuum surface force (CSF) surface tension model (`tension.h`).

---

## Physical Parameters

```c
#define g_    981.0
#define rho_w 1.0
#define rho_a (rho_w/1000.0)
#define L     75.5996
```

The lower (heavier) fluid has density $\rho_l = 1\ \mathrm{g/cm^3}$ and the upper (lighter) fluid $\rho_u = 10^{-3}\ \mathrm{g/cm^3}$, giving a density ratio

```math
\rho_r = \frac{\rho_u}{\rho_l} = 10^{-3}.
```

The computational domain is a square of side $L = 75.5996\ \mathrm{cm}$, centred at the origin.

---

## Boundary Conditions

```c
u.t[top]    = dirichlet(26.7046);
u.n[top]    = dirichlet(0.);

u.t[bottom] = dirichlet(26.7046);
u.n[bottom] = dirichlet(0.);

u.n[left]   = dirichlet(26.7046);
u.t[left]   = dirichlet(0.);

u.n[right]  = dirichlet(26.7046);
u.t[right]  = dirichlet(0.);

p[right]  = dirichlet(0.);
pf[right] = dirichlet(0.);

p[left]   = neumann(0.);
pf[left]  = neumann(0.);

p[top]    = neumann(0.);
pf[top]   = neumann(0.);

p[bottom]  = neumann(0.);
pf[bottom] = neumann(0.);

f[left]   = neumann(0.);
f[right]  = neumann(0.);
f[top]    = neumann(0.);
f[bottom] = neumann(0.);
```

A uniform horizontal base velocity $U = 26.7046\ \mathrm{cm/s}$ is imposed as a Dirichlet condition on the tangential velocity component at all four boundaries. The normal velocity is zero everywhere. The pressure is fixed to zero at the right boundary; zero-normal-gradient (Neumann) conditions are applied at all other boundaries. The volume fraction $f$ satisfies Neumann conditions at all boundaries to prevent artificial flux at the domain edges.

---

## Main Function

```c
int main()
{
  origin (-L/2., -L/2.);

  G.y = -g_;

  rho1 = rho_w;
  rho2 = rho_a;

  f.sigma = 72.0;

  size (L);
  init_grid (256);

  run();
}
```

The domain origin is placed at the domain centre $(-L/2, -L/2)$. Gravity acts downward ($G_y = -981\ \mathrm{cm/s^2}$). The lower-fluid density is `rho1` $= \rho_l$, the upper-fluid density is `rho2` $= \rho_u$, and the surface tension coefficient is $T = 72\ \mathrm{dyn/cm}$. The simulation starts on a uniform $256 \times 256$ base grid and calls `run()` to advance through all registered events. The simulation is inviscid (viscosity terms are commented out).

---

## Event: Pressure Forcing (`acceleration`)

```c
event acceleration (i++)
{
  face vector ia = a;
  scalar pressure_phi[];

  double T_      = 72.0;
  double prefac_ = 0.01;

  foreach_face(y) {

    double b_ = Delta;

    pressure_phi[] =
      (prefac_*T_*b_)/(pi*(sq(b_) + sq(x)));

    ia.y[] += alpha.y[]/(fm.y[] + SEPS)*
      pressure_phi[]*(f[] - f[0,-1])/Delta;
  }
}
```

A Lorentzian pressure distribution is applied at the interface at every time step:

```math
p_e(x) = \frac{F_0}{\pi}\frac{b}{b^2 + x^2}, \qquad F_0 = 0.01\,T.
```

The half-width $b = \Delta$ equals the local grid-cell size. In the limit $b \to 0$ this approaches a Dirac delta forcing. The forcing is implemented as an additional acceleration on the vertical face-velocity component, weighted by the VOF gradient $(f[] - f[0,-1])/\Delta$ which is non-zero only at the interface — confining the forcing to the two-fluid interface.

---

## Event: Initial Condition (`init`)

```c
event init (t = 0)
{
  fraction (f, -y);
  boundary ({f});

  foreach() {
    u.x[] = 26.7046;
    u.y[] = 0.;
  }

  boundary ({u, p});
}
```

The interface is initialized as a flat horizontal surface at $y = 0$: the `fraction` function sets $f = 1$ for $y < 0$ (lower fluid) and $f = 0$ for $y > 0$ (upper fluid). Both fluids start with the uniform horizontal velocity $u_x = U = 26.7046\ \mathrm{cm/s}$ and zero vertical velocity $u_y = 0$.

---

## Event: Adaptive Mesh Refinement (`adapt`)

```c
double uemax = 0.001;

#define fErr   (1e-3)
#define KErr   (1e-6)
#define VelErr (uemax)

#define MAXLEVEL 13

scalar KAPPA[];

event adapt (i++)
{
  curvature (f, KAPPA);

  foreach() {
    if (KAPPA[] == nodata)
      KAPPA[] = 0.;

    if (x < -30. || x > 30.)
      KAPPA[] = 0.;
  }

  boundary ((scalar *) {KAPPA});

  adapt_wavelet ((scalar *) {f, u.x, u.y, KAPPA},
                 (double[]) {fErr, VelErr, VelErr, KErr},
                 MAXLEVEL, 5);
}
```

At every time step, the mesh is adapted using wavelet-based refinement on the volume fraction $f$ (threshold $10^{-3}$), both velocity components $u_x$, $u_y$ (threshold $10^{-3}$), and the interface curvature $\kappa$ (threshold $10^{-6}$). The maximum refinement level is **13** (minimum **5**), corresponding to an effective maximum resolution of $2^{13} = 8192$ cells per side. Curvature-based refinement is suppressed for $|x| > 30\ \mathrm{cm}$ (the buffer regions near the left and right boundaries), preventing unnecessary refinement far from the active interface region. Cells where the curvature is undefined (`nodata`) are also set to zero before adaptation.

---

## Event: Output (`dump_i1`)

```c
event dump_i1 (t = 0; t += 0.01; t <= 5.02)
{
  char name4[128];
  sprintf (name4, "./data/snapshot-%g", t);

  scalar pid[];

  foreach()
    pid[] = fmod(pid()*(npe() + 37), npe());

  boundary ({pid});

  dump (name4);
}
```

The complete simulation state is written to `./data/snapshot-<t>` every $0.01\ \mathrm{s}$, from $t = 0$ to $t = 5.02\ \mathrm{s}$. Each dump file contains the full field data (velocity, pressure, volume fraction, and mesh) and can be restarted or post-processed independently. The `pid[]` scalar records the MPI process ID for each cell, useful for visualizing the domain decomposition in parallel runs. The interface profiles at selected times are extracted from these dump files and stored as [`if_1.csv`](https://github.com/yewalenikhil65/ForcedInterfacialWaves.jl/blob/main/notebooks/if_1.csv), [`if_3.csv`](https://github.com/yewalenikhil65/ForcedInterfacialWaves.jl/blob/main/notebooks/if_3.csv), [`if_7.csv`](https://github.com/yewalenikhil65/ForcedInterfacialWaves.jl/blob/main/notebooks/if_7.csv), [`if_15.csv`](https://github.com/yewalenikhil65/ForcedInterfacialWaves.jl/blob/main/notebooks/if_15.csv), [`if_25.csv`](https://github.com/yewalenikhil65/ForcedInterfacialWaves.jl/blob/main/notebooks/if_25.csv), [`if_60.csv`](https://github.com/yewalenikhil65/ForcedInterfacialWaves.jl/blob/main/notebooks/if_60.csv), [`if_145.csv`](https://github.com/yewalenikhil65/ForcedInterfacialWaves.jl/blob/main/notebooks/if_145.csv), and [`if_300.csv`](https://github.com/yewalenikhil65/ForcedInterfacialWaves.jl/blob/main/notebooks/if_300.csv), used in the [Capillary–Gravity](capillary_gravity.md#comparison-with-nonlinear-simulations-fig-10) comparison.

---

## Running the Simulation

```bash
# Compile (serial)
qcc -O2 -o capillary_gravity_forced capillary_gravity_forced.c -lm
mkdir -p data
./capillary_gravity_forced

# Compile and run in parallel (MPI)
qcc -O2 -D_MPI=1 -o capillary_gravity_forced capillary_gravity_forced.c -lm
mpirun -np <N> ./capillary_gravity_forced
```

---

## References

[^1]: Popinet, S., & collaborators. (2013–2026). *Basilisk: Free software for solving partial differential equations on adaptive Cartesian meshes*. [http://basilisk.fr](http://basilisk.fr)
