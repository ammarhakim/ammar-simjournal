#include <gkgeom.h>

void
psi(double t, const double *xn, double *fout, void *ctx)
{
  double R = xn[0], Z = xn[1];
  fout[0] = (R-2)*(R-2) + Z*Z/4;
}

void
basic_root_test(struct gkgeom_inp *inp, const struct gkgeom_app *app)
{
  double psi0 = 10.0, Z0 = 0.0;
  //double psi0 = 1.0, Z0 = 0.0;
  
  long nloop = 1;
  int nroots;
  double R[2], dR[2];
  struct timespec tm = gkyl_wall_clock();
  for (long i=0; i<nloop; ++i)
    nroots = gkgeom_app_R_psiz(app, psi0, Z0, 2, R, dR);
  
  printf("Took %g sec to run %ld roots (polyOrder %d on [%d X %d] grid\n",
    gkyl_time_diff_now_sec(tm), nloop,
    inp->polyOrder, inp->cells[0], inp->cells[1]
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

  double Ipsi = gkgeom_app_integrate_psi_contour(app, inp->upper[1],
    inp->lower[1], inp->upper[1], psi0);
  printf("Length of countour with psi=%lg is %lg\n", psi0, Ipsi);
  printf(" Error: %lg\n", Ipsi - 8.174873682157976);
}

int
main(void)
{
  struct gkgeom_inp inp = {
    .name = "gkgeom",

    .polyOrder = 2,
    .lower = { 0.5, -4.0 },
    .upper = { 6.0, 4.0 },
    .cells = { 64, 128 },

    .use_proj_on_basis = false,

    .psi = psi,
    .ctx = 0
  };

  gkgeom_app *app = gkgeom_app_new(&inp);

  //basic_root_test(&inp, app);

  // Computational grid: theta X psi X alpha (only 2D for now)
  double lower[] = { -M_PI/2, 7.0 };
  double upper[] = { M_PI/2, 12.0 };
  int cells[] = { 16, 8 };

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
