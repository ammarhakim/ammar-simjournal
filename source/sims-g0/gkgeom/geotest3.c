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

  fout[0] = psi_prefactor*(0.119946404324151*hex(x)*log(x) - 0.191087723964864*hex(x) - 1.43935685188981*qad(x)*sq(y)*log(x) + 1.69332066595762*qad(x)*sq(y) + 0.4138062761496*qad(x)*y*log(x) - 0.484565034576353*qad(x)*y + 0.66221050132831*qad(x)*log(x) - 0.0701095967287596*qad(x) + 0.959571234593207*sq(x)*qad(y)*log(x) - 0.569130557125708*sq(x)*qad(y) - 0.551741701532799*sq(x)*cub(y)*log(x) + 0.232280436618871*sq(x)*cub(y) - 2.64884200531324*sq(x)*sq(y)*log(x) - 1.12869311706989*sq(x)*sq(y) + 0.712672354321272*sq(x)*y*log(x) + 0.400592495501521*sq(x)*y + 0.319254744035335*sq(x)*log(x) + 0.144471037299054*sq(x) - 0.0639714156395472*hex(y) + 0.0551741701532799*pen(y) + 0.441473667552207*qad(y) - 0.237557451440424*cub(y) - 0.396754744035335*sq(y) + 0.0795862418377478*y + 0.0831914907404484);
}

int
main(void)
{
  struct solovev_ctx ctx = {
    .R0 = 1.0, .psi_prefactor = 1.0
  };

  // psi_max = 0.0075

  struct gkgeom_inp inp = {
    .name = "cerfon",

    .polyOrder = 2,
    .lower = { 0.5, -1.0 },
    .upper = { 1.5, 0.75 },
    .cells = { 64, 128 },

    .psi = psi,
    .ctx = &ctx,
  };

  gkgeom_app *app = gkgeom_app_new(&inp);

  //double psi0 = 0.05, Z0 = 0.1;
  double psi0 = 0.005, Z0 = 0.1;

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

  double Ipsi = gkgeom_app_integrate_psi_contour(app, inp.lower[1], 0.0, psi0);
  printf("Length of countour with psi=%g is %lg\n", psi0, Ipsi);
  //printf(" Error: %lg percent \n", (Ipsi - 3.343782219297654)/4.052162043995597*100);
  //printf(" Error: %lg percent \n", (Ipsi - 4.052162043995597)/4.052162043995597*100);

  gkgeom_app_release(app);
  
  return 0;
}
