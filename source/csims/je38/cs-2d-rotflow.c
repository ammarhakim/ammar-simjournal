#include <adiff.h>
#include <math.h>
#include <string.h>

#include <gkyl_const.h>

static inline double sq(double x) { return x*x; }

static void
velocity(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  // solid-body rotation about (0.5, 0.5)
  double ux = -y+0.5;
  double uy = x-0.5;
  fout[0] = ux; fout[1] = uy;
}

static void
init(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double r0 = 0.2;
  double x0 = 0.25, y0 = 0.5;
  double r = fmin(sqrt(sq(x-x0) + sq(y-y0)), r0)/r0;
  fout[0] = 0.25*(1 + cos(GKYL_PI*r));
}

struct sim_params {
  const char *name;
  int nx;
  double cfl_frac;
  double ux;
  double alpha;
  enum adiff_advection_scheme a_scheme;
  enum adiff_diffusion_scheme d_scheme;
};

static void
run_sim(struct sim_params params)
{
  struct adiff_app_inp app_inp = {
    .ndim = 2,
    .cells = { 64, 64 },
    .lower = { 0.0, 0.0 },
    .upper = { 1.0, 1.0 },

    .nframe = 1,
    .tend = 2.0*GKYL_PI,
    .cfl_frac = 0.9,

    .a_scheme = params.a_scheme,
    .d_scheme = params.d_scheme,

    .init = init,
    .velocity = velocity,

    .alpha = params.alpha
  };
  strcpy(app_inp.name, params.name);

  adiff_app *app = adiff_app_new(&app_inp);
  adiff_app_run(app);
  adiff_app_release(app);
}

int
main(void)
{

  // upwind scheme scans
  run_sim( (struct sim_params) {
      .name = "cs-2d-rotflow-c",
      .alpha = 0.0,
      .a_scheme = ADV_SCHEME_C2,
      .d_scheme = DIF_SCHEME_C2,
    }
  );

  // upwind scheme scans
  run_sim( (struct sim_params) {
      .name = "cs-2d-rotflow-u",
      .alpha = 0.0,
      .a_scheme = ADV_SCHEME_U3,
      .d_scheme = DIF_SCHEME_C2,
    }
  );  
  
  return 0;
}
