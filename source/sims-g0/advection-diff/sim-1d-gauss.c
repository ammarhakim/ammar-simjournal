#include <adiff.h>

// velocity field
static void
velocity(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1.0;
}

// function to set initial conditions
static void
init(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = exp(-50*x*x);
}

int
main(int argc, char **argv)
{
  struct adiff_app_inp app_inp = {
    .name = "sim-1d-gauss",

    .ndim = 1,
    .cells = { 32 },
    .lower = { -1.0 },
    .upper = { 1.0 },

    .nframe = 5,
    .tend = 2.0,
    .cfl_frac = 0.9,

    .scheme = SCHEME_C2_C2,

    .init = init,
    .velocity = velocity,

    .alpha = 0.0
  };

  adiff_app *app = adiff_app_new(&app_inp);
  adiff_app_run(app);
  adiff_app_release(app);
  
  return 0;
}
