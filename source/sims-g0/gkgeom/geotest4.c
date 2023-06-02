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
    .cells = { 64, 129 },

    .psi = psi,
    .ctx = &ctx,
  };

  gkgeom_app *app = gkgeom_app_new(&inp);

  int npsi = 10;
  double psi_min = 0.001, psi_max = 1.2;
  double dpsi = (psi_max-psi_min)/npsi;
  
  // Computational grid: theta X psi X alpha (only 2D for now)
  double lower[] = { -M_PI/2, psi_min };
  double upper[] = { M_PI/2, psi_min+dpsi };
  int cells[] = { 8, 1 };

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
