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

    .polyOrder = 1,
    .lower = {0.1, -4.0},
    .upper = {6.0, 4.0},
    .cells = { 64, 128 },

    .psi = psi,
    .ctx = 0
  };

  gkgeom_app *app = gkgeom_app_new(&inp);

  double psi0 = 4.0, Z0 = -3.0;

  long nloop = 1;
  double R[2], dR[2];
  struct timespec tm = gkyl_wall_clock();
  for (long i=0; i<nloop; ++i)
    gkgeom_app_R_psiz(app, psi0, Z0, R, dR);
  printf("Took %g sec to run %ld roots\n", gkyl_time_diff_now_sec(tm), nloop);
  
  printf("R(%g,%g) : %g, %g\n", psi0, Z0, R[0], R[1]);
  // check solution
  double psiL[1], psiR[1];
  psi(0.0, (double[]) { R[0], Z0 }, psiL, 0);
  psi(0.0, (double[]) { R[1], Z0 }, psiR, 0);

  printf("psi(%g,%g) = %g\n", R[0], Z0, psiL[0]);
  printf("psi(%g,%g) = %g\n", R[1], Z0, psiR[0]);

  gkgeom_app_release(app);
  
  return 0;
}
