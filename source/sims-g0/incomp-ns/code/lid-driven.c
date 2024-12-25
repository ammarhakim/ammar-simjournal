#include <incompns.h>

struct eval_ctx {
  
};

void
rho_init(double t, const double *xn, double *fout, void *ctx)
{
  struct eval_ctx *incomp_ctx = ctx;  
  double x = xn[0], y = xn[1];

  double rho = 1.0;
  if (x < 0.5 && y < 0.5)
    rho = 0.5;
  if (x < 0.5 && y > 0.5)
    rho = 0.25;
  if (x > 0.5 && y < 0.5)
    rho = 0.75;
  if (x > 0.5 && y > 0.5)
    rho = 1.0;

  fout[0] = rho;
}

void
vel_init(double t, const double *xn, double *fout, void *ctx)
{
  struct eval_ctx *incomp_ctx = ctx;  
  double x = xn[0], y = xn[1];

  fout[0] = 0.0;
  fout[1] = 0.0;
}

int
main(void)
{
  struct eval_ctx ctx = {
  };
  
  struct incompns_inp inp = {
    .name = "lid-driven",
    .ndim = 2,

    .lower = { 0.0, 0.0 },
    .upper = { 1.0, 1.0 },
    .cells = { 64, 64 },

    .cfl_frac = 0.9,

    .rho_init = rho_init,
    .vel_init = vel_init,
    .ctx = &ctx,

    .bcx = { GKYL_SPECIES_NO_SLIP, GKYL_SPECIES_NO_SLIP },
    .bcy = { GKYL_SPECIES_NO_SLIP, GKYL_SPECIES_NO_SLIP },

    .bcy_vel = { 0.0, 1.0 }, // drive upper wall
  };

  struct incomp_run_params run_param = {
    .nframe = 1,
    .tstart = 0.0,
    .tend = 1.0    
  };

  incompns_app *app = incompns_app_new(&inp);
  incompns_run(app, run_param);
  incompns_app_release(app);
  
  return 0;
}
