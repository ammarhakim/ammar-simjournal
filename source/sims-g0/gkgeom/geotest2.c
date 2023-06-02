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

void
basic_root_test(struct gkgeom_inp *inp, const struct gkgeom_app *app)
{
  /* //double psi0 = 0.05, Z0 = 0.1; */
  /* double psi0 = 0.1, Z0 = 0.1; */

  /* long nloop = 1; */
  /* int nroots; */
  /* double R[2], dR[2]; */
  /* struct timespec tm = gkyl_wall_clock(); */
  /* for (long i=0; i<nloop; ++i) */
  /*   nroots = gkgeom_app_R_psiz(app, psi0, Z0, 2, R, dR); */
  
  /* printf("Took %g sec to run %ld roots (polyOrder %d on [%d X %d] grid\n", */
  /*   gkyl_time_diff_now_sec(tm), nloop, */
  /*   inp.polyOrder, inp.cells[0], inp.cells[1] */
  /* ); */
  /* printf("psi = %g. Z = %g. Number of roots = %d\n", psi0, Z0, nroots); */

  /* for (int i=0; i<nroots; ++i) { */
  /*   printf("R(%g,%g) : %g\n", psi0, Z0, R[i]); */
  /*   printf("dR/dZ(%g,%g) : %g\n", psi0, Z0, dR[i]); */

  /*   double psiE[1]; */
  /*   psi(0.0, (double[]) { R[i], Z0 }, psiE, &ctx); */
  /*   printf("  psi(%g,%g) = %g\n", R[i], Z0, psiE[0]); */
  /* } */

  /* double Ipsi = gkgeom_app_integrate_psi_contour(app, inp.upper[1], */
  /*   inp.lower[1], inp.upper[1], psi0); */
  /* printf("Length of countour with psi=%g is %lg\n", psi0, Ipsi); */
  /* printf(" Error: %lg percent \n", (Ipsi - 3.343782219297654)/4.052162043995597*100); */
  /* //printf(" Error: %lg percent \n", (Ipsi - 4.052162043995597)/4.052162043995597*100);   */
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

  // Computational grid: theta X psi X alpha (only 2D for now)
  double lower[] = { -M_PI/2, 0.05 };
  double upper[] = { M_PI/2, 0.1 };
  int cells[] = { 16, 10 };

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 2, lower, upper, cells);

  struct gkyl_range clocal, clocal_ext;
  int nghost[GKYL_MAX_CDIM] = { 0, 0 };
  gkyl_create_grid_ranges(&cgrid, nghost, &clocal_ext, &clocal);

  int cpoly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 2, 1);

  struct gkyl_array *mapc2p = gkyl_array_new(GKYL_DOUBLE, 2*cbasis.num_basis, clocal_ext.volume);
  gkgeom_app_calcgeom(app, &cgrid, cpoly_order, mapc2p);

  long ncalls = gkgeom_app_nroots(app);
  printf("Total number of R(psi,Z) calls: %ld\n", ncalls);

  char fileNm[1024]; // buffer for file name
  do {
    const char *fmt = "%s_mapc2p.gkyl";
    snprintf(fileNm, sizeof fileNm, fmt, inp.name);  
    gkyl_grid_sub_array_write(&cgrid, &clocal, mapc2p, fileNm);
  } while (0);
  
  gkyl_array_release(mapc2p);
  gkgeom_app_release(app);
  
  return 0;
}
