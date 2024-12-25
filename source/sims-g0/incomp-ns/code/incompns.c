#include <incompns.h>
#include <mpack.h>

#include <limits.h>

struct gkyl_msgpack_data*
incompns_array_meta_new(int nframe, double tm)
{
  struct gkyl_msgpack_map_elem elist[] = {
    { .key = "frame", .elem_type = GKYL_MP_INT, .ival = nframe },
    { .key = "time", .elem_type = GKYL_MP_DOUBLE, .dval = tm }
  };
  
  return gkyl_msgpack_create(2, elist);
}

// Incompressible NS App
struct incompns_app {
  char name[128];
  int ndim;
  double tcurr;
  double cfl;

  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext;
  struct gkyl_range global, global_ext;

  struct gkyl_rect_decomp *decomp;
  struct gkyl_comm *comm;

  double nu, kappa;

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions

  // non-periodic BCs
  enum gkyl_species_bc_type bcx[2], bcy[2], bcz[2];
  // for driven simulations, set velocity of fluid at boundary
  double bcx_vel[2], bcy_vel[2], bcz_vel[2];  

  struct gkyl_array *rho, *rho0;
  struct gkyl_array *vel, *vel0;

  evalf_t rho_init, vel_init;
  void *ctx;
};

struct incompns_app *
incompns_app_new(struct incompns_inp *inp)
{
  struct incompns_app *app = gkyl_malloc(sizeof *app);
  strcpy(app->name, inp->name);
  int ndim = app->ndim = inp->ndim;

  app->tcurr = 0.0;
  double cfl_frac = inp->cfl_frac > 0 ? inp->cfl_frac : 1.0;

  gkyl_rect_grid_init(&app->grid, inp->ndim, inp->lower, inp->upper, inp->cells);

  int nghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&app->grid, nghost, &app->global_ext, &app->global);

  int cuts[3] = { 1, 1, 1 };
  app->decomp = gkyl_rect_decomp_new_from_cuts(ndim, cuts, &app->global);
    
  app->comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = app->decomp,
      .sync_corners = true,
      }
  );
    
  // global and local ranges are same, and so just copy them.
  memcpy(&app->local, &app->global, sizeof(struct gkyl_range));
  memcpy(&app->local_ext, &app->global_ext, sizeof(struct gkyl_range));  

  app->num_periodic_dir = inp->num_periodic_dir;
  for (int d=0; d<ndim; ++d)
    app->periodic_dirs[d] = inp->periodic_dirs[d];

  for (int d=0; d<2; ++d) {
    app->bcx_vel[d] = inp->bcx_vel[d];
    app->bcy_vel[d] = inp->bcy_vel[d];
    app->bcz_vel[d] = inp->bcz_vel[d];
  }

  app->nu = inp->nu;
  app->kappa = inp->kappa;

  app->rho_init = inp->rho_init;
  app->vel_init = inp->vel_init;
  app->ctx = inp->ctx;

  app->rho = gkyl_array_new(GKYL_DOUBLE, 1, app->local_ext.volume);
  app->rho0 = gkyl_array_new(GKYL_DOUBLE, 1, app->local_ext.volume);

  app->vel = gkyl_array_new(GKYL_DOUBLE, ndim, app->local_ext.volume);
  app->vel0 = gkyl_array_new(GKYL_DOUBLE, ndim, app->local_ext.volume);  

  return app;
}

void
incompns_app_apply_ic(incompns_app *app)
{
  // initialize density
  struct gkyl_offset_descr offset_rho[] = {
    { 0.0, 0.0, 0.0 },
  };
  gkyl_eval_offset_fd *ev_rho = gkyl_eval_offset_fd_new(
    &(struct gkyl_eval_offset_fd_inp) {
      .grid = &app->grid,
      .num_ret_vals = 1,
      .eval = app->rho_init,
      .ctx = app->ctx,
      .offsets = offset_rho
    }
  );
  gkyl_eval_offset_fd_advance(ev_rho, 0.0, &app->local, app->rho);
  gkyl_eval_offset_fd_release(ev_rho);

  // initialize velocity
  struct gkyl_offset_descr offset_vel[] = {
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0 },
  };
  gkyl_eval_offset_fd *ev_vel = gkyl_eval_offset_fd_new(
    &(struct gkyl_eval_offset_fd_inp) {
      .grid = &app->grid,
      .num_ret_vals = app->ndim,
      .eval = app->vel_init,
      .ctx = app->ctx,
      .offsets = offset_vel
    }
  );
  gkyl_eval_offset_fd_advance(ev_vel, 0.0, &app->local, app->vel);
  gkyl_eval_offset_fd_release(ev_vel);
}


// apply boundary condition to density
static void
apply_rho_bc(incompns_app *app, struct gkyl_array *rho)
{
  int num_periodic_dir = app->num_periodic_dir, ndim = app->ndim, is_non_periodic[3] = {1, 1, 1};
  
  for (int d=0; d<num_periodic_dir; ++d)
    is_non_periodic[app->periodic_dirs[d]] = 0;

  for (int d=0; d<ndim; ++d)
    if (is_non_periodic[d]) {

    }

  // sync interior ghost cells
  gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, rho);
  // sync periodic ghost cells
  gkyl_comm_array_per_sync(app->comm, &app->local, &app->local_ext, num_periodic_dir,
    app->periodic_dirs, rho);
}

// apply boundary condition to velocity
static void
apply_vel_bc(incompns_app *app, struct gkyl_array *vel)
{
  int num_periodic_dir = app->num_periodic_dir, ndim = app->ndim, is_non_periodic[3] = {1, 1, 1};
  
  for (int d=0; d<num_periodic_dir; ++d)
    is_non_periodic[app->periodic_dirs[d]] = 0;

  for (int d=0; d<ndim; ++d)
    if (is_non_periodic[d]) {

    }

  // sync interior ghost cells
  gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, vel);
  // sync periodic ghost cells
  gkyl_comm_array_per_sync(app->comm, &app->local, &app->local_ext, num_periodic_dir,
    app->periodic_dirs, vel);  
}

void
incompns_write(incompns_app *app, double tm, int frame)
{
  char fileNm[1024]; // buffer for file name

  struct gkyl_msgpack_data *mdata = incompns_array_meta_new(frame, tm);
  
  do {
    const char *fmt = "%s_rho_%d.gkyl";
    snprintf(fileNm, sizeof fileNm, fmt, app->name, frame);
    gkyl_grid_sub_array_write(&app->grid, &app->local, mdata, app->rho, fileNm);
  } while (0);

  do {
    const char *fmt = "%s_vel_%d.gkyl";
    snprintf(fileNm, sizeof fileNm, fmt, app->name, frame);
    gkyl_grid_sub_array_write(&app->grid, &app->local, mdata, app->vel, fileNm);
  } while (0);

  gkyl_msgpack_data_release(mdata);
}

struct gkyl_update_status
incompns_update(incompns_app *app, double dt)
{
  struct gkyl_update_status status = {
    .dt_actual = dt,
    .dt_suggested = dt,
    .success = true
  };

  return status;
}

static void
write_data(struct gkyl_tm_trigger* iot, incompns_app* app, double t_curr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    incompns_write(app, t_curr, frame);
  }
}

bool
incompns_run(incompns_app *app, struct incomp_run_params run_param)
{
  double t_curr = run_param.tstart, t_end = run_param.tend;
  int num_frames = run_param.nframe;
  int frame_curr = 0;
  
  struct gkyl_tm_trigger io_trig = { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };

  int max_steps = INT_MAX;
  if (run_param.use_single_step)
    max_steps = run_param.max_step;
  
  incompns_app_apply_ic(app);
  incompns_write(app, 0.0, frame_curr);

  double dt = t_end-t_curr;
  
  long step = 1;
  while ((t_curr < t_end) && (step <= max_steps)) {
    fprintf(stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    
    struct gkyl_update_status status = incompns_update(app, dt);
    fprintf(stdout, " dt = %g\n", status.dt_actual);
    
    if (!status.success) {
      fprintf(stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }

    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    write_data(&io_trig, app, t_curr);

    step += 1;
  }
  write_data(&io_trig, app, t_curr);

  return true;
}

void
incompns_app_release(incompns_app *app)
{
  gkyl_comm_release(app->comm);
  gkyl_rect_decomp_release(app->decomp);
  
  gkyl_array_release(app->rho);
  gkyl_array_release(app->rho0);

  gkyl_array_release(app->vel);
  gkyl_array_release(app->vel0);
  
  gkyl_free(app);
}
