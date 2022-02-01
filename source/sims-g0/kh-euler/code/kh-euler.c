#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <rxi_ini.h>
#include <gkylzero.h>

struct euler_ctx {
    int cells[2]; // grid resolution
    double tend; // time to run sim
    int nframe; // number of frames to write
    double gas_gamma; // gas constant
    enum gkyl_wave_limiter limiter; // limiter to use
    pcg32_random_t rng; // RNG for use in IC
};

// Names of limiters to enum mapping
static const struct { char const *const nm; enum gkyl_wave_limiter lm; } lm_names[] = {
  { "min-mod", GKYL_MIN_MOD },
  { "superbee", GKYL_SUPERBEE },
  { "van-leer", GKYL_VAN_LEER },
  { "monotonized-centered", GKYL_MONOTONIZED_CENTERED },
  { "beam-warming", GKYL_BEAM_WARMING }
};

// Return limiter code given then name of the limiter
static enum gkyl_wave_limiter
get_limiter(const char *nm)
{
  int n = sizeof(lm_names)/sizeof(lm_names[0]);
  for (int i=0; i<n; ++i)
    if (strcmp(lm_names[i].nm, nm) == 0)
      return lm_names[i].lm;
  return GKYL_MONOTONIZED_CENTERED;
}

void
evalEulerInit(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  struct euler_ctx *app = ctx;
  double gas_gamma = app->gas_gamma;

  double x = xn[0], y = xn[1];

  double rho = 1.0, vx = 0.5, pr = 2.5;
  if (fabs(y)<0.25) {
    rho = 2.0;
    vx = -0.5;
  }

  vx = vx + 0.01*2*(0.5*gkyl_pcg32_rand_double(&app->rng)-1);
  double vy = 0.01*2*(0.5*gkyl_pcg32_rand_double(&app->rng)-1);
  
  fout[0] = rho;
  fout[1] = rho*vx; fout[2] = rho*vy; fout[3] = 0.0;
  fout[4] = 0.5*rho*(vx*vx+vy*vy) + pr/(gas_gamma-1);
}

struct euler_ctx
create_ctx(rxi_ini_t *inp)
{
  struct euler_ctx ctx = {
    .gas_gamma = 1.4,
    .nframe = 1,
    .limiter = GKYL_MONOTONIZED_CENTERED,
    .rng = gkyl_pcg32_init(true)
  };

  int read_failed = 0;
  rxi_ini_sget(inp, "params", "gas_gamma", "%lg", &ctx.gas_gamma);
  rxi_ini_sget(inp, "params", "nframe", "%d", &ctx.nframe);

  if (!rxi_ini_sget(inp, "params", "cells", "%d, %d", &ctx.cells[0], &ctx.cells[1] )) {
    fprintf(stderr, "Must provide 'cells = NX, NY'!\n");
    read_failed = 1;
  }
  if (!rxi_ini_sget(inp, "params", "tend", "%lg", &ctx.tend)) {
    fprintf(stderr, "Must provide 'tend'!\n");
    read_failed = 1;
  }

  const char *lm_nm = 0;
  if ( (lm_nm = rxi_ini_get(inp, "params", "limiter")) )
    ctx.limiter = get_limiter(lm_nm);

  if (read_failed) {
    fprintf(stderr, "... aborting!\n");
    exit(1);
  }
  
  return ctx;
}

void
write_data(struct gkyl_tm_trigger *iot, const gkyl_moment_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr))
    gkyl_moment_app_write(app, tcurr, iot->curr-1);
}

int
main(int argc, char **argv)
{
  const char *inp_name = argc>1 ? argv[1] : "kh-euler.ini";
  rxi_ini_t *inp = rxi_ini_load(inp_name);

  if (0 == inp) {
    fprintf(stderr, "Unable to open input file %s!\n", inp_name);
    exit(1);
  }

  struct euler_ctx ctx = create_ctx(inp); // context for init functions  

  // equation object
  struct gkyl_wv_eqn *euler = gkyl_wv_euler_new(ctx.gas_gamma);

  struct gkyl_moment_species fluid = {
    .name = "euler",

    .equation = euler,
    .limiter = ctx.limiter,
    
    .evolve = 1,
    .ctx = &ctx,
    .init = evalEulerInit,
  };

  // VM app
  struct gkyl_moment app_inp = {

    .ndim = 2,
    .lower = { -0.5, -0.5 },
    .upper = { 0.5, 0.5 }, 
    .cells = { ctx.cells[0], ctx.cells[1] },

    .num_species = 1,
    .species = { fluid },

    .num_periodic_dir = 2,
    .periodic_dirs = { 0, 1 },    
  };

  // construct sim name based on input file name
  char const *const inp_last_slash = strrchr(inp_name, '/');
  char const *const inp_no_slash = inp_last_slash ? inp_last_slash+1 : inp_name;
  strncpy(app_inp.name, inp_no_slash, strcspn(inp_no_slash, ".ini"));

  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = ctx.tend/ctx.nframe };

  // create app object
  gkyl_moment_app *app = gkyl_moment_app_new(app_inp);

  // start, end
  double tcurr = 0.0, tend = ctx.tend;

  // initialize simulation
  gkyl_moment_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);

  // compute estimate of maximum stable time-step
  double dt = gkyl_moment_app_max_dt(app);

  long step = 1;
  while (tcurr < tend) {
    printf("Taking time-step %ld at t = %g ...", step, tcurr);
    struct gkyl_update_status status = gkyl_moment_update(app, dt);
    printf(" dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      printf("** Update method failed! Aborting simulation ....\n");
      break;
    }
    tcurr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, tcurr);

    step += 1;
  }
  gkyl_moment_app_stat_write(app);
  
  struct gkyl_moment_stat stat = gkyl_moment_app_stat(app);

  // simulation complete, free resources
  gkyl_wv_eqn_release(euler);
  gkyl_moment_app_release(app);
  rxi_ini_free(inp);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of failed time-steps %ld\n", stat.nfail);
  printf("Species updates took %g secs\n", stat.species_tm);
  printf("Field updates took %g secs\n", stat.field_tm);
  printf("Total updates took %g secs\n", stat.total_tm);
  
  return 0;
}
