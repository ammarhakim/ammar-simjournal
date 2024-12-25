#include <gkgeom.h>

struct solovev_ctx {
  double R0, psi_prefactor;
};

static inline double sq(double x) { return x*x; }
static inline double cub(double x) { return x*x*x; }
static inline double qad(double x) { return x*x*x*x; }
static inline double pen(double x) { return x*x*x*x*x; }
static inline double hex(double x) { return x*x*x*x*x*x; }

void
psi(double t, const double *xn, double *fout, void *ctx)
{
  struct solovev_ctx *s = ctx;
  double R0 = s->R0, psi_prefactor = s->psi_prefactor;
  double R = xn[0], Z = xn[1];
  double x = R/R0, y = Z/R0;

  fout[0] = psi_prefactor*(0.00373804283369699*hex(x)*log(x) - 0.00574955335438162*hex(x) - 0.0448565140043639*qad(x)*sq(y)*log(x) + 0.0503044260840946*qad(x)*sq(y) + 0.017623348727471*qad(x)*log(x) + 0.0956643504553683*qad(x) + 0.0299043426695759*sq(x)*qad(y)*log(x) - 0.0160920841654771*sq(x)*qad(y) - 0.0704933949098842*sq(x)*sq(y)*log(x) + 0.0644725519961135*sq(x)*sq(y) - 7.00898484784405e-5*sq(x)*log(x) - 0.303766642191745*sq(x) - 0.00199362284463839*hex(y) + 0.0117488991516474*qad(y) + 7.00898484784405e-5*sq(y) + 0.0145368720253975);
}

void
test_arc(struct gkgeom_inp *inp, gkgeom_app *app)
{
  double psi_curr = 0.001;

  double zmin = inp->lower[1], zmax = inp->upper[1];
  double zmid = 0.5*(zmin+zmax);
  double zqua = 0.5*(zmin+zmid);

  double arcL = gkgeom_app_integrate_psi_contour(app, 6.0, zmin, zmax, psi_curr);
  double arcL2 = gkgeom_app_integrate_psi_contour(app, 6.0, zmin, zmid , psi_curr);
  double arcL3 = gkgeom_app_integrate_psi_contour(app, 6.0, zmid, zmax, psi_curr);

  printf("arcL = %.10lg\n", arcL);  
  printf("arcL2 = %.10lg\n", arcL2);
  printf("arcL3 = %.10lg\n", arcL3);
  printf("arcL2 + arcL3 = %lg (%lg)\n", arcL2 + arcL3, arcL2+arcL3-arcL);

  int nr;
  double R[2], dR[2];
  nr = gkgeom_app_R_psiz(app, psi_curr, 0.0-1e-14, 2, R, dR);
  printf("nr: %d. R = %lg. dR = %lg\n", nr, R[1], dR[1]);

  nr = gkgeom_app_R_psiz(app, psi_curr, 0.0+1e-14, 2, R, dR);
  printf("nr: %d. R = %lg. dR = %lg\n", nr, R[1], dR[1]);  
}

int
main(void)
{
  struct solovev_ctx ctx = {
    .R0 = 2.5, .psi_prefactor = 1.0
  };

  // psi_max = 0.0075

  struct gkgeom_inp inp = {
    .name = "cerfon_dn",

    .polyOrder = 2,

    .lower = { 0.01, -6.0 },
    .upper = { 6.0, 6.0 },

    .cells = { 64, 128 },

    .psi = psi,
    .ctx = &ctx,

    .use_proj_on_basis = false
  };

  gkgeom_app *app = gkgeom_app_new(&inp);

  //test_arc(&inp, app);
  
  int npsi = 10;
  double psi_min = 0.0001, psi_max = 1.2;
  double dpsi = (psi_max-psi_min)/npsi;
  
  // Computational grid: theta X psi X alpha (only 2D for now)
  double lower[] = { -M_PI/2, psi_min };
  double upper[] = { M_PI/2, psi_max };
  int cells[] = { 16, npsi };

  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, 2, lower, upper, cells);

  struct gkyl_range clocal, clocal_ext;
  int nghost[GKYL_MAX_CDIM] = { 0, 0 };
  gkyl_create_grid_ranges(&cgrid, nghost, &clocal_ext, &clocal);

  int cpoly_order = 2;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 2, cpoly_order);

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
