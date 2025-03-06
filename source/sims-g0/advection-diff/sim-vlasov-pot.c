#include <adiff.h>
#include <math.h>

static inline double sq(double x) { return x*x; }

// velocity field
static void
velocity(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0], v = xn[1];
  // characteristics for Vlasov with a fixed potential of
  // phi(x) = cos(x)
  double ux = v; // characteristic vel in conf space
  double Ex = sin(x); // Ex = -grad phi
  fout[0] = ux; fout[1] = Ex;
}

// function to set initial conditions
static void
init(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0], v = xn[1];
  // Mawellian with vth = 1
  fout[0] = 1/sqrt(2*GKYL_PI)*exp(-v*v/2);
}

int
main(void)
{
  struct adiff_app_inp app_inp = {
    .name = "sim-vlasov-pot",

    .ndim = 2,
    .cells = { 64, 64 },
    // [0, 2 pi] x [-6vth, 6vth]
    .lower = { 0.0, -6.0 },
    .upper = { 2*GKYL_PI, 6.0 },

    .nframe = 16,
    .tend = 3.0,
    .cfl_frac = 0.9,

    .scheme = SCHEME_C2_C2,

    .init = init,
    .velocity = velocity,

    .alpha = 0.0,
  };

  adiff_app *app = adiff_app_new(&app_inp);
  adiff_app_run(app);
  adiff_app_release(app);
  
  return 0;
}
