/**
 * run_gc_ivp.c
 *
 * Step 4 of the Basilisk initialisation pipeline.
 *
 * Initialises the VOF fraction from FreeSurface.dat and the velocity
 * field from velocity_interpolated_below.dat, then runs the gravity–
 * capillary IVP with adaptive mesh refinement.
 *
 * The velocity file is a flat CSV (u,v per line) matching the
 * foreach() traversal order for cells with f[] >= 0.5.
 *
 * Physical setup (CGS):
 *   - g = 981 cm/s², ρ_water = 1.0 g/cm³, ρ_air = 0.001 g/cm³
 *   - σ = 72 dyn/cm (surface tension)
 *   - Inviscid (μ = 0)
 *   - Domain: 2λ × 2λ, periodic in x, free-slip top/bottom
 *   - Adaptive grid: maxlevel=11, minlevel=7
 *
 * Compile:
 *   qcc -O2 -o run_gc_ivp run_gc_ivp.c -lm
 *   (MPI: CC99='mpicc -std=c99' qcc -D_MPI=1 -O2 -o run_gc_ivp run_gc_ivp.c -lm)
 *
 * Follows: TravellingViscousGC/basilisk_sim/run_StokesWave.c
 */

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "reduced.h"
#include "tension.h"
#include "distance.h"

#include "params.h"

#define L       (2.0 * WAVELENGTH)

double t_final = 50.0;
double saveAt  = 0.01;
double uemax   = 0.01;
int    maxlevel = 11;

/* free-slip BCs */
u.t[top]    = neumann(0.);
u.n[top]    = dirichlet(0.);
u.t[bottom] = neumann(0.);
u.n[bottom] = dirichlet(0.);

void put_velocity(FILE *fpt);

int main()
{
    origin(-L / 2., -L / 2.);
    periodic(right);
    G.y = -G_ACCEL;
    rho1 = RHO_W;
    rho2 = RHO_W / 1000.0;
    f.sigma = SIGMA_VAL;
    mu1 = 0.;
    mu2 = 0.;
    size(L);
    init_grid(RESOLUTION);
    system("mkdir -p dumpfile interface_data");
    run();
}

event init(i = 0)
{
    /* --- Set VOF fraction from FreeSurface.dat --- */
    FILE *ptr1 = fopen("FreeSurface.dat", "rb");
    if (ptr1 == NULL) {
        fprintf(ferr, "ERROR: FreeSurface.dat not found\n");
        exit(1);
    }
    coord *shapedata = input_xy(ptr1);
    fclose(ptr1);

    scalar d[];
    distance(d, shapedata);

    vertex scalar phi[];
    foreach_vertex() {
        double xx = (d[] + d[-1] + d[0, -1] + d[-1, -1]) / 4.0;
        phi[] = xx;
    }
    fractions(phi, f);

    /* --- Read and apply velocity field --- */
    FILE *fpt = fopen("velocity_interpolated_below.dat", "r");
    if (fpt == NULL) {
        fprintf(ferr, "ERROR: velocity_interpolated_below.dat not found\n");
        exit(1);
    }
    put_velocity(fpt);
    fclose(fpt);

    boundary((scalar *){u});
}

event logfile2(i++)
{
    if (pid() == 0)
        fprintf(ferr, "t:%g \t dt:%g \n", t, dt);
}

int pt = 0;
event snapshot(t = 0; t += saveAt; t <= t_final)
{
    pt = round(t / saveAt);
    char name[80];
    sprintf(name, "./dumpfile/dump-%d", pt);
    dump(name);
    char name_1[80];
    sprintf(name_1, "./interface_data/interface-%d.dat", pt);
    FILE *ptr_1 = fopen(name_1, "w");
    output_facets(f, ptr_1);
    fclose(ptr_1);
}

event adapt(i++)
{
    scalar KAPPA[];
    curvature(f, KAPPA);
    boundary((scalar *){KAPPA});
    adapt_wavelet({KAPPA, f, u}, (double[]){1e-4, 0.001, uemax, uemax, uemax},
                  maxlevel, 7);
}

/**
 * Read velocity from flat CSV (u,v per line) in foreach() traversal order.
 * Only reads for cells with f[] >= 0.5 (water phase).
 * This MUST match the same foreach() order used by dump_coordinates.c.
 */
void put_velocity(FILE *fpt)
{
    char str[250];
    if (fpt == NULL) {
        printf("Can't open velocity file\n");
        return;
    }
    int count = 0;
    foreach() {
        if (f[] >= 0.5) {
            if (fgets(str, 250, fpt)) {
                char *ptr = strtok(str, ",");
                int column = 1;
                while (ptr != NULL) {
                    if (column == 1)
                        u.x[] = atof(ptr);
                    else if (column == 2)
                        u.y[] = atof(ptr);
                    else {
                        printf("ERROR: unexpected column in velocity file\n");
                        fflush(stdout);
                        exit(1);
                    }
                    column++;
                    ptr = strtok(NULL, ",");
                }
                count++;
            }
        }
    }
    fprintf(ferr, "put_velocity: loaded %d cell velocities\n", count);
}
