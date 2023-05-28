#include <gkgeom.h>

void
psi(double t, const double *xn, double *fout, void *ctx)
{
  double R = xn[0], Z = xn[1];
  fout[0] = (R-2)*(R-2) + Z*Z/4;
}

int
main(void)
{
  struct gkgeom_inp inp = {
    .name = "gkgeom",

    .polyOrder = 2,
    .lower = { 0.1, -4.0 },
    .upper = { 6.0, 4.0 },
    .cells = { 64, 128 },

    .use_proj_on_basis = false,

    .psi = psi,
    .ctx = 0
  };

  gkgeom_app *app = gkgeom_app_new(&inp);

  double psi0 = 2.0, Z0 = -2.0;

  long nloop = 1;
  int nroots;
  double R[2], dR[2];
  struct timespec tm = gkyl_wall_clock();
  for (long i=0; i<nloop; ++i)
    nroots = gkgeom_app_R_psiz(app, psi0, Z0, 2, R, dR);
  
  printf("Took %g sec to run %ld roots (polyOrder %d on [%d X %d] grid\n",
    gkyl_time_diff_now_sec(tm), nloop,
    inp.polyOrder, inp.cells[0], inp.cells[1]
  );
  printf("psi = %g. Z = %g. Number of roots = %d\n", psi0, Z0, nroots);

  for (int i=0; i<nroots; ++i) {
    printf("R(%g,%g) : %g\n", psi0, Z0, R[i]);
    printf("dR/dZ(%g,%g) : %g\n", psi0, Z0, dR[i]);

    double psiE[1];
    psi(0.0, (double[]) { R[i], Z0 }, psiE, 0);
    printf("  psi(%g,%g) = %g\n", R[i], Z0, psiE[0]);
    printf("  dR/dZ(%g,%g) = %g\n", R[i], Z0, -0.25*Z0/(R[i]-2));
  }

  gkgeom_app_release(app);
  
  return 0;
}
