#include <adiff.h>

// velocity field
void
velocity(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double ux = 1.0, uy = 0.0;

  fout[0] = ux; fout[1] = uy;
}

// function to set initial conditions
void
init(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = exp(-50*x*x);
}

int
main(void)
{
  struct adiff_app_inp app_inp = {
    .name = "sim-1d-const-adiff",
    .cells = { 32, 32 },
    .lower = { -1.0, -1.0 },
    .upper = { 1.0, 1.0 },

    .nframe = 5,
    .tend = 2.0,
    .cfl_frac = 0.9,

    .init = init,
    .velocity = velocity,
  };

  adiff_app *app = adiff_app_new(&app_inp);
  adiff_app_run(app);
  adiff_app_release(app);
  
  return 0;
}
