#include <gkgeom.h>

struct solovev_ctx {
  double B0, R0, k, q0, Ztop;
};

static inline double sq(double x) { return x*x; }

void
psi(double t, const double *xn, double *fout, void *ctx)
{
  struct solovev_ctx *s = ctx;
  double B0 = s->B0, R0 = s->R0, k = s->k, q0 = s->q0;
  double R = xn[0], Z = xn[1];
  fout[0] = B0*k/(2*sq(R0)*q0)*(sq(R)*sq(Z)/sq(k) + sq(sq(R) - sq(R0))/4);
}

int
main(void)
{
  struct solovev_ctx ctx = {
    .B0 = 0.55, .R0 = 0.85, .k = 2, .q0 = 2, .Ztop = 1.5
  };

  double psi_sep = ctx.B0*ctx.k*sq(ctx.R0)/(8*ctx.q0);
  printf("psi_sep = %lg\n", psi_sep);
  
  struct gkgeom_inp inp = {
    .name = "solovev",

    .polyOrder = 2,
    .lower = { 0.0, -1.5 },
    .upper = { 1.5, 1.5 },
    .cells = { 64, 128 },

    .psi = psi,
    .ctx = &ctx,
  };

  gkgeom_app *app = gkgeom_app_new(&inp);

  //double psi0 = 0.05, Z0 = 0.1;
  double psi0 = 0.1, Z0 = 0.1;

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
    psi(0.0, (double[]) { R[i], Z0 }, psiE, &ctx);
    printf("  psi(%g,%g) = %g\n", R[i], Z0, psiE[0]);
  }

  double Ipsi = gkgeom_app_integrate_psi_contour(app, inp.lower[1], inp.upper[1], psi0);
  printf("Length of countour with psi=%g is %lg\n", psi0, Ipsi);
  printf(" Error: %lg percent \n", (Ipsi - 3.343782219297654)/4.052162043995597*100);
  //printf(" Error: %lg percent \n", (Ipsi - 4.052162043995597)/4.052162043995597*100);

  gkgeom_app_release(app);
  
  return 0;
}
