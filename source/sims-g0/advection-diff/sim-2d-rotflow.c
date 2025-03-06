#include <adiff.h>
#include <math.h>

static inline double sq(double x) { return x*x; }

// velocity field
static void
velocity(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  // solid-body rotation about (0.5, 0.5)
  double ux = -y+0.5;
  double uy = x-0.5;
  fout[0] = ux; fout[1] = uy;
}

// function to set initial conditions
static void
init(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double r0 = 0.2;
  double x0 = 0.25, y0 = 0.5;
  double r = fmin(sqrt(sq(x-x0) + sq(y-y0)), r0)/r0;
  fout[0] = 0.25*(1 + cos(GKYL_PI*r));
}

int
main(void)
{
  struct adiff_app_inp app_inp = {
    .name = "sim-2d-rotflow",

    .ndim = 2,
    .cells = { 64, 64 },
    .lower = { 0.0, 0.0 },
    .upper = { 1.0, 1.0 },

    .nframe = 16,
    .tend = 2.0*GKYL_PI,
    .cfl_frac = 0.9,

    .scheme = SCHEME_C2_C2,

    .init = init,
    .velocity = velocity,

    .alpha = 1e-5,
  };

  adiff_app *app = adiff_app_new(&app_inp);
  adiff_app_run(app);
  adiff_app_release(app);
  
  return 0;
}
