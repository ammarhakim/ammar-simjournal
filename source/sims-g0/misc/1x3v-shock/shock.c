#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_vlasov.h>
#include <rt_arg_parse.h>
#include <rxi_ini.h>

struct shock_inp {
  int cells[4]; // grid dimensions
  double lower[4], upper[4]; // grid lower/upper
  double tend; // end-time
  int nframe; // number of frames to write
  double mass; // mass of species
  double n; // number density of jet
  double u; // velocity of jet
  double vth; // thermal speed of jet
};

static inline double sq(double x) { return x*x; }

inline double
maxwellian(double n, double vx, double vy, double vz, double u, double vth)
{
  double v2 = (vx - u)*(vx - u) + vy*vy + vz*vz;
  return n/pow(sqrt(2*M_PI*vth*vth),3)*exp(-v2/(2*vth*vth));
}

void
evalDistFunc(double t, const double* GKYL_RESTRICT xn,
  double* GKYL_RESTRICT fout, void *ctx)
{
  struct shock_inp *app = ctx;
  double x = xn[0], vx = xn[1], vy = xn[2], vz = xn[3];
  double n = app->n;
  double u = app->u;
  double vth = app->vth;
  if (x<0)
    fout[0] = maxwellian(n, vx, vy, vz, u, vth);
  else
    fout[0] = maxwellian(n, vx, vy, vz, -u, vth);
}

void
evalNu(double t, const double* GKYL_RESTRICT xn, double* GKYL_RESTRICT fout,
  void* ctx)
{
  struct free_stream_ctx *app = ctx;
  //double x = xn[0], v = xn[1];
  fout[0] = 15e3/5e-4;
}

struct shock_inp
create_shock_inp(rxi_ini_t *inp)
{
  struct shock_inp sinp;

  int read_failed = 0;
  // read from input file
  if (!rxi_ini_sget(inp, "conf-grid", "cells", "%d", &sinp.cells[0])) {
    fprintf(stderr, "Must provide 'cells' in section '[conf-grid]'!\n");
    read_failed = 1;
  }
  if (!rxi_ini_sget(inp, "conf-grid", "lower", "%lg", &sinp.lower[0])) {
    fprintf(stderr, "Must provide 'lower' in section '[conf-grid]'!\n");
    read_failed = 1;
  }
  if (!rxi_ini_sget(inp, "conf-grid", "upper", "%lg", &sinp.upper[0])) {
    fprintf(stderr, "Must provide 'upper' in section '[conf-grid]'!\n");
    read_failed = 1;
  }  
  if (!rxi_ini_sget(inp, "conf-grid", "tend", "%lg", &sinp.tend)) {
    fprintf(stderr, "Must provide 'tend' in section '[conf-grid]'!\n");
    read_failed = 1;
  }
  if (!rxi_ini_sget(inp, "conf-grid", "nframe", "%d", &sinp.nframe)) {
    fprintf(stderr, "Must provide 'nframe' in section '[conf-grid]'!\n");
    read_failed = 1;
  }  

  if (!rxi_ini_sget(inp, "neutrals", "mass", "%lg", &sinp.mass)) {
    fprintf(stderr, "Must provide 'mass' in section '[neutrals]'!\n");
    read_failed = 1;
  }
  if (!rxi_ini_sget(inp, "neutrals", "n", "%lg", &sinp.n)) {
    fprintf(stderr, "Must provide 'n' in section '[neutrals]'!\n");
    read_failed = 1;
  }
  if (!rxi_ini_sget(inp, "neutrals", "u", "%lg", &sinp.u)) {
    fprintf(stderr, "Must provide 'u' in section '[neutrals]'!\n");
    read_failed = 1;
  }
  if (!rxi_ini_sget(inp, "neutrals", "vth", "%lg", &sinp.vth)) {
    fprintf(stderr, "Must provide 'vth' in section '[neutrals]'!\n");
    read_failed = 1;
  }
  if (!rxi_ini_sget(inp, "neutrals", "cells", "%d, %d, %d", &sinp.cells[1], &sinp.cells[2], &sinp.cells[3])) {
    fprintf(stderr, "Must provide 'cells' in section '[neutrals]'!\n");
    read_failed = 1;
  }
  if (!rxi_ini_sget(inp, "neutrals", "lower", "%lg, %lg, %lg", &sinp.lower[1], &sinp.lower[2], &sinp.lower[3])) {
    fprintf(stderr, "Must provide 'lower' in section '[neutrals]'!\n");
    read_failed = 1;
  }
  if (!rxi_ini_sget(inp, "neutrals", "upper", "%lg, %lg, %lg", &sinp.upper[1], &sinp.upper[2], &sinp.upper[3])) {
    fprintf(stderr, "Must provide 'upper' in section '[neutrals]'!\n");
    read_failed = 1;
  }

  if (read_failed) {
    fprintf(stderr, "... aborting!\n");
    exit(1);
  }  

  return sinp;
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_vlasov_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr)) {
    gkyl_vlasov_app_write(app, tcurr, iot->curr-1);
    gkyl_vlasov_app_calc_mom(app);
    gkyl_vlasov_app_write_mom(app, tcurr, iot->curr-1);
  }
}

int
main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }

  if (strcmp(app_args.file_name, APP_ARGS_DEFAULT_FILE_NAME) == 0)
    strcpy(app_args.file_name, "shock.ini");

  struct shock_inp sinp;
  
  rxi_ini_t *inp = rxi_ini_load(app_args.file_name);
  if (inp)
    sinp = create_shock_inp(inp);
  else
    exit(1);
  
  // electrons
  struct gkyl_vlasov_species neut = {
    .name = "neut",
    .charge = 0.0, .mass = sinp.mass,
    
    .lower = { sinp.lower[1], sinp.lower[2], sinp.lower[3] },
    .upper = { sinp.upper[1], sinp.upper[2], sinp.upper[3] },
    .cells = { sinp.cells[1], sinp.cells[2], sinp.cells[3] },

    .evolve = 1,
    .ctx = &sinp,
    .init = evalDistFunc,
    .nu = evalNu,
    .collision_id = GKYL_LBO_COLLISIONS,
    .num_diag_moments = 3,
    .diag_moments = { "M0", "M1i", "M2" },
  };

  // VM app
  struct gkyl_vm vm = {
    .name = "shock_g0",

    .cdim = 1, .vdim = 3,
    .lower = { sinp.lower[0] },
    .upper = { sinp.upper[0] },
    .cells = { sinp.cells[0] },
    .poly_order = 2,
    .basis_type = app_args.basis_type,

    .num_periodic_dir = 0,
    .periodic_dirs = {  },

    .num_species = 1,
    .species = { neut },
    .skip_field = true,

    .use_gpu = app_args.use_gpu,
  };

  // create app object
  gkyl_vlasov_app *app = gkyl_vlasov_app_new(vm);

  // start, end and initial time-step
  double tcurr = 0.0, tend = sinp.tend;
  double dt = tend-tcurr;

  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = tend/sinp.nframe };

  // initialize simulation
  gkyl_vlasov_app_apply_ic(app, tcurr);
  write_data(&io_trig, app, tcurr);

  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps)) {
    printf("Taking time-step at t = %g ...", tcurr);
    struct gkyl_update_status status = gkyl_vlasov_update(app, dt);
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

  gkyl_vlasov_app_write(app, tcurr, 1);
  gkyl_vlasov_app_calc_mom(app); gkyl_vlasov_app_write_mom(app, tcurr, 1);
  gkyl_vlasov_app_stat_write(app);

  // fetch simulation statistics
  struct gkyl_vlasov_stat stat = gkyl_vlasov_app_stat(app);

  // simulation complete, free app
  rxi_ini_free(inp);
  gkyl_vlasov_app_release(app);

  printf("\n");
  printf("Number of update calls %ld\n", stat.nup);
  printf("Number of forward-Euler calls %ld\n", stat.nfeuler);
  printf("Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    printf("Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    printf("Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }  
  printf("Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  printf("Species RHS calc took %g secs\n", stat.species_rhs_tm);
  printf("Species collisions took %g secs\n", stat.species_coll_tm);
  printf("Updates took %g secs\n", stat.total_tm);
  
  return 0;
}
