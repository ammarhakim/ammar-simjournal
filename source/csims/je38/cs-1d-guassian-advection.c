#include <adiff.h>
#include <math.h>
#include <string.h>

#include <gkyl_const.h>

static inline double sq(double x) { return x*x; }

struct velocity_ctx { double ux; };

static void
velocity(double t, const double *xn, double *fout, void *ctx)
{
  struct velocity_ctx *vctx = ctx;
  fout[0] = vctx->ux;
}

static void
init(double t, const double *xn, double *fout, void *ctx)
{
  double x = xn[0];
  double sig = 0.05;
  fout[0] = exp(-x*x/sq(sig));
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
  struct velocity_ctx vctx = { .ux = params.ux };
  
  struct adiff_app_inp app_inp = {
    .ndim = 1,
    .cells = { params.nx },
    .lower = { -0.5 },
    .upper = { 0.5 },

    .nframe = 1,
    .tend = 1.0,
    .cfl_frac = params.cfl_frac,

    .a_scheme = params.a_scheme,
    .d_scheme = params.d_scheme,

    .init = init,

    .velocity_ctx = &vctx,
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
      .name = "cs-1d-gaussian-advection-nx-16",
      .nx = 16,
      .cfl_frac = 1.0-2, // to reduce time-step error
      .ux = 1.0,
      .alpha = 0.0,      
      .a_scheme = ADV_SCHEME_U3,
      .d_scheme = DIF_SCHEME_C2,
    }
  );

  run_sim( (struct sim_params) {
      .name = "cs-1d-gaussian-advection-nx-32",
      .nx = 32,
      .cfl_frac = 1.0-2, // to reduce time-step error
      .ux = 1.0,
      .alpha = 0.0,      
      .a_scheme = ADV_SCHEME_U3,
      .d_scheme = DIF_SCHEME_C2,
    }
  );

  run_sim( (struct sim_params) {
      .name = "cs-1d-gaussian-advection-nx-64",
      .nx = 64,
      .cfl_frac = 1.0-2, // to reduce time-step error
      .ux = 1.0,
      .alpha = 0.0,      
      .a_scheme = ADV_SCHEME_U3,
      .d_scheme = DIF_SCHEME_C2,
    }
  );

  run_sim( (struct sim_params) {
      .name = "cs-1d-gaussian-advection-nx-128",
      .nx = 128,
      .cfl_frac = 1.0-2, // to reduce time-step error
      .ux = 1.0,
      .alpha = 0.0,      
      .a_scheme = ADV_SCHEME_U3,
      .d_scheme = DIF_SCHEME_C2,
    }
  );  

  // central scheme scans  
  run_sim( (struct sim_params) {
      .name = "cs-1d-gaussian-advection-nx-c32",
      .nx = 32,
      .cfl_frac = 1.0-2, // to reduce time-step error
      .ux = 1.0,
      .alpha = 0.0,      
      .a_scheme = ADV_SCHEME_C2,
      .d_scheme = DIF_SCHEME_C2,
    }
  );

  run_sim( (struct sim_params) {
      .name = "cs-1d-gaussian-advection-nx-c64",
      .nx = 64,
      .cfl_frac = 1.0-2, // to reduce time-step error
      .ux = 1.0,
      .alpha = 0.0,      
      .a_scheme = ADV_SCHEME_C2,
      .d_scheme = DIF_SCHEME_C2,
    }
  );

  run_sim( (struct sim_params) {
      .name = "cs-1d-gaussian-advection-nx-c128",
      .nx = 128,
      .cfl_frac = 1.0-2, // to reduce time-step error
      .ux = 1.0,
      .alpha = 0.0,      
      .a_scheme = ADV_SCHEME_C2,
      .d_scheme = DIF_SCHEME_C2,
    }
  );  

  // central scheme L2 scans
  run_sim( (struct sim_params) {
      .name = "cs-1d-gaussian-advection-nx-l2-c64-a",
      .nx = 64,
      .cfl_frac = 0.9,
      .ux = 1.0,
      .alpha = 0.0,      
      .a_scheme = ADV_SCHEME_C2,
      .d_scheme = DIF_SCHEME_C2,
    }
  );

  run_sim( (struct sim_params) {
      .name = "cs-1d-gaussian-advection-nx-l2-c64-b",
      .nx = 64,
      .cfl_frac = 0.9/2.0,
      .ux = 1.0,
      .alpha = 0.0,      
      .a_scheme = ADV_SCHEME_C2,
      .d_scheme = DIF_SCHEME_C2,
    }
  );  

  run_sim( (struct sim_params) {
      .name = "cs-1d-gaussian-advection-nx-l2-c64-c",
      .nx = 64,
      .cfl_frac = 0.9/4.0,
      .ux = 1.0,
      .alpha = 0.0,      
      .a_scheme = ADV_SCHEME_C2,
      .d_scheme = DIF_SCHEME_C2,
    }
  );  
  
  return 0;
}
